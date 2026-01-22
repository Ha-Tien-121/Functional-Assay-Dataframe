import pandas as pd
import pathlib
import re
import ast
import os
import difflib

# Cache for loaded file columns to avoid re-reading large CSVs
# Key: file_path, Value: list of columns
COLUMN_CACHE = {}

def get_file_columns(file_path):
    """
    Reads the file and attempts to find the correct header row.
    Handles .csv, .csv.gz, .xlsx, .xls
    Scans first few rows to find the one with the most columns/strings.
    """
    if file_path in COLUMN_CACHE:
        return COLUMN_CACHE[file_path]
        
    if not os.path.exists(file_path):
        print(f"  [Error] File not found: {file_path}")
        return []

    try:
        # Determine extension
        ext = pathlib.Path(file_path).suffix.lower()
        df = None
        
        # Read first few rows
        if ext == '.gz':
            df = pd.read_csv(file_path, nrows=5, header=None)
        elif ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, nrows=5, header=None)
        elif ext in ['.tsv']:
             df = pd.read_csv(file_path, sep='\t', nrows=5, header=None)
        else:
            # Assume CSV
            try:
                df = pd.read_csv(file_path, nrows=5, header=None)
            except:
                # Fallback for weird encodings or separators
                try:
                    df = pd.read_csv(file_path, nrows=5, header=None, encoding='latin1')
                except:
                    pass

        if df is None or df.empty:
            return []

        # Heuristic to find header row
        best_header_idx = 0
        max_cols = 0
        
        for i in range(len(df)):
            row = df.iloc[i]
            # Count non-null string-like values
            # meaningful column names usually don't contain NaN and are strings
            non_nulls = row.count()
            
            # Simple heuristic: row with most non-nulls is likely header
            # Tie-breaker: prefer row with more string values
            if non_nulls > max_cols:
                max_cols = non_nulls
                best_header_idx = i
            elif non_nulls == max_cols:
                # Check stringiness
                current_strs = sum(1 for x in df.iloc[best_header_idx] if isinstance(x, str))
                new_strs = sum(1 for x in row if isinstance(x, str))
                if new_strs > current_strs:
                    best_header_idx = i
        
        # Extract header
        header_row = df.iloc[best_header_idx]
        cols = [str(c).strip() for c in header_row if pd.notna(c)]
        
        # If the file has a "Gene" column or similar, it's a strong indicator
        # But we just return what we found
        
        COLUMN_CACHE[file_path] = cols
        return cols
    except Exception as e:
        print(f"  [Error] Failed to read columns from {file_path}: {e}")
        return []

def parse_build_script(script_path):
    """
    Parses a build_[gene]_dataframe.py script to extract:
    1. File variables (e.g. FINDLAY_FILE -> 'path/to/file')
    2. Explicit column renames from load_* functions.
    """
    with open(script_path, 'r') as f:
        content = f.read()
    
    # 1. Extract file variables and their values
    # Regex to find XXX_FILE = INPUT_DIR / "..." or similar
    # We'll just look for assignments where the value contains a string
    file_vars = {}
    
    # Simple regex to capture VAR = ...
    # We will then eval the RHS roughly
    
    # Regex: VAR_NAME = INPUT_DIR / "string"
    # or VAR_NAME = pathlib.Path("string")
    # We'll rely on the specific pattern in these files: INPUT_DIR / "path"
    
    # Pattern: VAR = INPUT_DIR / "path"
    pattern = r'^\s*([A-Z0-9_]+_FILE(?:_\d+)?)\s*=\s*INPUT_DIR\s*/\s*["\']([^"\']+)["\']'
    
    for line in content.splitlines():
        match = re.search(pattern, line)
        if match:
            var_name = match.group(1)
            rel_path = match.group(2)
            # Resolve to absolute or relative to CWD
            # Assuming INPUT_DIR is CWD
            file_vars[var_name] = rel_path.strip()

    # 2. Parse AST for renames (same as before)
    renames = {} # target_col -> (source_col, dataset_name)
    
    try:
        tree = ast.parse(content)
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                func_name = node.name
                if func_name.startswith('load_'):
                    dataset_guess = func_name.replace('load_', '').capitalize()
                    for subnode in ast.walk(node):
                        if isinstance(subnode, ast.Call) and isinstance(subnode.func, ast.Attribute) and subnode.func.attr == 'rename':
                            for keyword in subnode.keywords:
                                if keyword.arg == 'columns' and isinstance(keyword.value, ast.Dict):
                                    keys = []
                                    values = []
                                    for k in keyword.value.keys:
                                        if isinstance(k, ast.Constant): keys.append(k.value)
                                        elif isinstance(k, ast.Str): keys.append(k.s)
                                        else: keys.append(None)
                                    for v in keyword.value.values:
                                        if isinstance(v, ast.Constant): values.append(v.value)
                                        elif isinstance(v, ast.Str): values.append(v.s)
                                        else: values.append(None)
                                    for k, v in zip(keys, values):
                                        if k and v:
                                            renames[v] = (k, dataset_guess)
    except Exception as e:
        print(f"Error parsing AST for {script_path}: {e}")

    return file_vars, renames

def infer_dataset_name(col_name, gene, file_vars):
    """
    Infers the dataset name from the column name and available file variables.
    """
    if col_name.upper().startswith(gene.upper() + '_'):
        rest = col_name[len(gene) + 1:]
    else:
        rest = col_name
        
    sorted_vars = sorted(list(file_vars.keys()), key=len, reverse=True)
    
    for var in sorted_vars:
        # Strip _FILE suffix for matching against column name
        var_base = var.replace('_FILE', '')
        # Handle LEE_FILE_1 -> LEE
        var_base = re.sub(r'_\d+$', '', var_base)
        
        if var_base in rest.upper():
             parts = var_base.split('_')
             formatted = "_".join([p.capitalize() if not p.isdigit() else p for p in parts])
             
             year_match = re.search(r'20\d\d', rest)
             if year_match and year_match.group(0) not in formatted:
                 formatted += f"_{year_match.group(0)}"
                 
             return formatted
             
    match = re.match(r'^([A-Za-z]+)_(\d{4})_', rest)
    if match:
        return f"{match.group(1)}_{match.group(2)}"
        
    match = re.match(r'^([A-Za-z]+)_', rest)
    if match:
        return match.group(1)

    return "?"

def find_best_column_match(target_col, dataset_name, gene, file_vars):
    """
    Finds the best matching column in the files associated with the dataset_name.
    """
    # 1. Identify relevant files from file_vars
    # Dataset name might be "Findlay" or "Bouwman_2013"
    # file_vars keys are like "FINDLAY_FILE", "BOUWMAN_2013_FILE"
    
    candidate_files = []
    
    # Normalize dataset name for matching
    # "Bouwman_2013" -> "BOUWMAN_2013"
    # "Lee_2010" -> "LEE" (to catch LEE_FILE_1, etc.)
    
    ds_upper = dataset_name.upper()
    
    # Heuristic mapping from Dataset Name back to FILE keys
    for var_key, file_path in file_vars.items():
        # var_key: FINDLAY_FILE
        # ds_upper: FINDLAY
        
        # Check if ds_upper contains the var base or vice versa
        # Strip _FILE
        var_base = var_key.replace('_FILE', '')
        # Strip numeric suffix like _1
        var_base_root = re.sub(r'_\d+$', '', var_base)
        
        # Match checks:
        # 1. Exact match (BOUWMAN_2013 == BOUWMAN_2013)
        # 2. Root match (LEE == LEE_2010 root?) -> No, LEE_2010 contains LEE
        
        if var_base == ds_upper:
            candidate_files.append(file_path)
        elif var_base_root in ds_upper and len(var_base_root) > 3: # Avoid short matches
            candidate_files.append(file_path)
        elif ds_upper.replace('_', '') in var_base.replace('_', ''):
            candidate_files.append(file_path)
            
    if not candidate_files:
        # Fallback: try finding any file with the first part of dataset name
        # e.g. "Adamovich_Hdr" -> find "ADAMOVICH"
        prefix = ds_upper.split('_')[0]
        for var_key, file_path in file_vars.items():
            if prefix in var_key:
                candidate_files.append(file_path)
    
    candidate_files = list(set(candidate_files))
    
    if not candidate_files:
        return "?"

    # 2. Search in candidates
    best_match = None
    best_score = 0
    
    # Prepare target for matching
    # Remove Gene prefix
    clean_target = target_col
    if clean_target.upper().startswith(gene.upper() + '_'):
        clean_target = clean_target[len(gene) + 1:]
        
    # Remove Dataset prefix if present
    # dataset_name might be "Bouwman_2020" or "Findlay"
    # clean_target might be "Bouwman_2020_nGFP+b" or "Findlay_auth_func_score"
    
    # Try to strip dataset name
    # 1. Exact start match
    if clean_target.lower().startswith(dataset_name.lower() + '_'):
        clean_target = clean_target[len(dataset_name) + 1:]
    
    # 2. Try stripping parts (e.g. "Bouwman" from "Bouwman_2020") if the year is redundant or missing
    ds_parts = dataset_name.split('_')
    # If the first part is in the target, strip it and maybe the year
    if clean_target.lower().startswith(ds_parts[0].lower()):
        # This is a bit aggressive but often correct for this schema
        # But we must be careful not to strip "Score" if dataset was named "Score" (unlikely)
        # Let's try to match as much of the dataset name as possible
        pass 
    
    # Clean target for better matching (replace underscores with spaces for fuzzy match?)
    # "HDR_score" vs "HDR Score" -> "HDR score" vs "HDR Score"
    clean_target_spaced = clean_target.replace('_', ' ')
    
    for fpath in candidate_files:
        file_cols = get_file_columns(fpath)
        if not file_cols:
            continue
            
        # Matching Logic
        
        # 1. Exact match (case insensitive)
        for col in file_cols:
            if col.lower() == clean_target.lower():
                return col
            if col.lower() == target_col.lower():
                return col
            # Match with spaces vs underscores
            if col.replace(' ', '_').lower() == clean_target.lower():
                return col
                
        # 2. Suffix match
        # Target: auth_func_score (after strip)
        # Col: auth_func_score
        for col in file_cols:
            if clean_target.lower().endswith(col.lower()) and len(col) > 4:
                # Check if it's too generic
                if col.lower() in ['score', 'class', 'category', 'function', 'comment']:
                     pass
                else:
                    return col
                    
        # 3. Fuzzy match
        # Match clean_target (without dataset prefix) against cols
        matches = difflib.get_close_matches(clean_target, file_cols, n=1, cutoff=0.6)
        if matches:
            return matches[0]
            
        # Try spaced version
        matches = difflib.get_close_matches(clean_target_spaced, file_cols, n=1, cutoff=0.6)
        if matches:
            return matches[0]

        # 4. Keyword search (Best Guess)
        # If target has "score" and col has "score"
        target_lower = clean_target.lower()
        
        # Extract key terms from target
        keywords = []
        if 'score' in target_lower: keywords.append('score')
        if 'class' in target_lower: keywords.append('class')
        if 'category' in target_lower: keywords.append('category')
        if 'p_val' in target_lower or 'pval' in target_lower: keywords.append('pval')
        if 'func' in target_lower: keywords.append('func')
        
        # If we have keywords, look for col with most overlapping keywords
        best_kw_col = None
        max_kw_count = 0
        
        for col in file_cols:
            col_lower = col.lower()
            count = sum(1 for kw in keywords if kw in col_lower)
            if count > max_kw_count:
                max_kw_count = count
                best_kw_col = col
        
        if max_kw_count > 0:
            return best_kw_col

    return "?"

def main():
    csv_path = "Functional Assay Mapping - Sheet1.csv"
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    print(f"Reading {csv_path}...")
    df = pd.read_csv(csv_path)
    
    # Ensure columns exist
    if 'Mapped to ' not in df.columns:
        df['Mapped to '] = ""
    if 'Found in this dataset' not in df.columns:
        df['Found in this dataset'] = ""
        
    genes = df['Gene'].unique()
    
    filled_count = 0
    skipped_count = 0
    heuristics_count = 0
    unresolved_rows = []
    
    for gene in genes:
        if gene == 'All Genes': 
            continue
            
        script_name = f"build_{gene.lower()}_dataframe.py"
        if not os.path.exists(script_name):
            print(f"Warning: Script {script_name} not found. Skipping {gene}.")
            continue
            
        print(f"Parsing {script_name} for {gene}...")
        file_vars, renames = parse_build_script(script_name)
        
        # Process rows for this gene
        gene_mask = df['Gene'] == gene
        indices = df[gene_mask].index
        
        for idx in indices:
            row = df.loc[idx]
            col_name = str(row['Column Name'])
            
            # Check current status
            curr_mapped = str(row['Mapped to '])
            curr_found = str(row['Found in this dataset'])
            
            # Logic Phase 1: Ensure 'Found in this dataset' is populated (from previous step)
            dataset_name = curr_found
            mapped_col = curr_mapped
            
            needs_dataset = dataset_name in ['', 'nan', '?']
            needs_mapping = mapped_col in ['', 'nan', '?']
            
            if not needs_dataset and not needs_mapping:
                skipped_count += 1
                continue

            # 1. Infer Dataset if needed
            if needs_dataset:
                if col_name in renames:
                    _, dataset_name = renames[col_name]
                else:
                    dataset_name = infer_dataset_name(col_name, gene, file_vars)
                
                df.at[idx, 'Found in this dataset'] = dataset_name

            # 2. Find Mapped Column if needed
            if needs_mapping:
                # Try explicit rename first
                if col_name in renames:
                    mapped_col, _ = renames[col_name]
                else:
                    # Heuristic Search
                    if dataset_name != "?":
                        print(f"  Searching match for {col_name} in {dataset_name}...")
                        mapped_col = find_best_column_match(col_name, dataset_name, gene, file_vars)
                        if mapped_col != "?":
                            print(f"    -> Found: {mapped_col}")
                            heuristics_count += 1
                        else:
                            print(f"    -> No match found.")
                            unresolved_rows.append(f"{gene} - {col_name} ({dataset_name})")
                    else:
                        mapped_col = "?"
                
                df.at[idx, 'Mapped to '] = mapped_col
                filled_count += 1
            
    # Save
    print(f"Saving updated CSV to {csv_path}...")
    df.to_csv(csv_path, index=False)
    print(f"\nSummary:")
    print(f"  Filled/Updated rows: {filled_count}")
    print(f"  Successful Heuristic Matches: {heuristics_count}")
    print(f"  Skipped (already done): {skipped_count}")
    print(f"  Unresolved: {len(unresolved_rows)}")
    
    if unresolved_rows:
        print("\nUnresolved Rows (Sample):")
        for r in unresolved_rows[:10]:
            print(f"  {r}")

if __name__ == "__main__":
    main()

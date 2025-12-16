import re
import pandas as pd
import pathlib
import numpy as np
from variant_helpers import parse_hgvsp

# Centralized Join Key Map
JOIN_KEY_COLUMN_MAP = {
    # BRCA1
    "FINDLAY_FILE": ["hgvs_p", "variant", "mutant", "HGVSp", "HGVSp.", "join_key"], 
    "ADAMOVICH_HDR_FILE": ["hgvs_p", "variant", "mutant"],
    "ADAMOVICH_CISPLATIN_FILE": ["hgvs_p", "variant", "mutant"],
    "FERNANDES_FILE": ["variant", "hgvs_p", "HGVSp", "p_variant"],
    "BOUWMAN_2013_FILE": ["Base_Variant", "variant", "hgvs_p"],
    "BOUWMAN_2020_FILE": ["Protein change", "variant"],
    "CALECA_FILE": ["Protein Change", "variant"],
    "GOU_FILE": ["Amino_acid_change", "HGVS", "hgvs_p"],
    "FAYER_FILE": ["p_variant", "HGVS", "HGVSp.", "hgvs_p", "variant"],
    "BASSI_FILE": ["Amino_acid_change", "variant"],
    "LANGERUD_FILE": ["HGVSc.", "Variant"],
    "LEE_FILE_1": ["Variant", "hgvs_p"],
    "LEE_FILE_2": ["Variant", "hgvs_p"],
    "LEE_FILE_3": ["Variant", "hgvs_p"],
    "STARITA_FILE": ["variant", "hgvs_p"],
    "HART_FILE": ["AminoAcid", "variant", "hgvs_p"],
    
    # BRCA2
    "RICHARDSON_FILE": ["hgvs_p", "variant", "Variant"],
    "HU_FILE": ["hgvs_p", "variant"],
    "HUANG_FILE": ["hgvs_p", "variant", "Amino acid change (p.)"],
    "IKEGAMI_FILE": ["hgvs_p", "variant", "Variant"],
    "BISWAS_FILE": ["hgvs_p", "variant", "Variant"],
    "MESMAN_FILE": ["hgvs_p", "variant", "Amino acid"],
    "GUIDUGLI_FILE": ["hgvs_p", "variant", "HGVS_Amino_Acid"],
    "SAHU_2023_FILE": ["hgvs_p", "variant"],
    "SAHU_2025_FILE": ["hgvs_p", "variant"],
    
    # MSH2
    "JIA_FILE": ["hgvs_p", "variant"],
    "OLLODART_FILE": ["hgvs_p", "variant"],
    "BOUVET_FILE": ["hgvs_p", "variant"],
    
    # PTEN
    "MATREYEK_FILE": ["hgvs_p", "variant"],
    "MIGHELL_FILE": ["hgvs_p", "variant"],
    
    # TP53
    "FUNK_FILE": ["hgvs_p", "variant"],
    "KOTLER_FILE": ["hgvs_p", "variant"],
    "KAWAGUCHI_FILE": ["hgvs_p", "variant"],
    "KATO_FILE": ["hgvs_p", "variant"],
    "GIACOMELLI_FILE": ["hgvs_p", "variant"],
    "TP53_FUNCTIONAL_WORKSHEET": ["Protein change", "variant"],
    
    # VHL
    "BUCKLEY_FILE": ["hgvs_p", "variant"]
}

# Tokens to identify the specific header row for each file (if not row 0)
HEADER_TOKENS_MAP = {
    "BOUWMAN_2020_FILE": ["Protein change", "Cisplatin"],
    "HART_FILE": ["AminoAcid", "Functional Score"], 
    "FAYER_FILE": ["p_variant"],
    "RICHARDSON_FILE": ["HDR Score", "Variant"], 
    "HUANG_FILE": ["Amino acid change (p.)"], 
    "IKEGAMI_FILE": ["Variant"], # "Olaparib" is on next line
    "MESMAN_FILE": ["Complementation", "Nucleotide"], 
    "GUIDUGLI_FILE": ["HGVS_Amino_Acid"], 
}

def normalize_join_key_value(val):
    """
    Robustly parses a value into a standardized protein change key (e.g. 'K80Q', 'C218*', 'P283=').
    Handles 3-letter codes, synonymous, stops, parenthesis, etc.
    """
    if pd.isna(val):
        return pd.NA
        
    val_str = str(val).strip()
    if not val_str or val_str.lower() in ['nan', 'none', '']:
        return pd.NA

    # Handle parenthesis like p.(Ser4Phe) -> p.Ser4Phe
    if '(' in val_str and ')' in val_str:
        val_str = val_str.replace('(', '').replace(')', '')
        
    # 1. Try parse_hgvsp first (shared helper)
    # Ensure p. prefix if missing but looks like protein change
    if not val_str.startswith('p.') and not val_str.startswith('c.'):
        # Heuristic: if it starts with an AA code
        if re.match(r'^([A-Z][a-z]{2}|[A-Z])\d+', val_str):
             val_str = 'p.' + val_str
             
    _, _, _, key = parse_hgvsp(val_str)
    if key:
        return key
        
    # 2. Manual fallbacks
    clean = val_str.replace('p.', '')
    
    # Handle synonyms with '='
    if '=' in clean:
        match = re.match(r'^([A-Z])\d+=$', clean)
        if match:
            return clean 
            
    # Simple 1-letter format K80Q
    if re.match(r'^[A-Z]\d+([A-Z]|\*|=)$', clean):
        return clean
        
    return pd.NA

def build_join_key(df, dataset_name, join_key_cols=None):
    """
    Constructs a standardized join_key for the DataFrame.
    """
    # 1. Try explicit AA columns
    if {'aa_pos', 'aa_ref', 'aa_alt'}.issubset(df.columns):
        if not df[['aa_pos', 'aa_ref', 'aa_alt']].isna().all().all():
            print(f"  Building join_key from aa_ref/pos/alt columns...")
            return df.apply(
                lambda row: f"{row['aa_ref']}{int(row['aa_pos']) if pd.notna(row['aa_pos']) else ''}{row['aa_alt']}" 
                if pd.notna(row['aa_ref']) and pd.notna(row['aa_pos']) and pd.notna(row['aa_alt']) else pd.NA, 
                axis=1
            )
    
    # 2. Try map-based columns
    candidate_cols = []
    if join_key_cols:
        candidate_cols = join_key_cols
    elif dataset_name in JOIN_KEY_COLUMN_MAP:
        candidate_cols = JOIN_KEY_COLUMN_MAP[dataset_name]
    
    target_col = None
    for col in candidate_cols:
        matches = [c for c in df.columns if c.strip().lower() == col.lower()]
        if matches:
            target_col = matches[0]
            break
            
    if target_col:
        print(f"  Building join_key from '{target_col}'...")
        
        # Specific override for Langerud HGVSc (cDNA) matching
        if dataset_name == "LANGERUD_FILE" and (target_col.lower().startswith("hgvs") or target_col == "Variant"):
             return df[target_col].astype(str).str.strip().replace({"nan": pd.NA})

        return df[target_col].apply(normalize_join_key_value)
    
    # 3. Fail
    print(f"  Warning: No suitable column found for join_key in {dataset_name}. Candidates: {candidate_cols}. Columns: {list(df.columns)}")
    return pd.Series([pd.NA] * len(df), index=df.index)

def find_header_row(path, required_tokens, encoding='utf-8'):
    """
    Scans the first N rows to find a header row containing required tokens.
    Returns 0-based row index or None.
    """
    try:
        # Read first 50 lines
        with open(path, 'r', encoding=encoding, errors='replace') as f:
            lines = [f.readline() for _ in range(50)]
            
        for i, line in enumerate(lines):
            # Check if line contains a significant number of tokens
            # ALL tokens must be present to be strict
            if all(token.lower() in line.lower() for token in required_tokens):
                return i
    except Exception:
        pass
    return None

def post_process_dataset(df, name):
    """
    Applies dataset-specific column renaming or fixes after loading.
    """
    if name == "BISWAS_FILE":
        # Rename columns to match mapping sheet
        rename_map = {
            "PIF_HAT_DS": "PIF[HAT+DS]a",
            "Functional_Classification": "Functional classification"
        }
        df = df.rename(columns=rename_map)
        
    if name == "IKEGAMI_FILE":
        # Rename columns to "Column X" (A=0, B=1, etc.) to match mapping sheet
        # Mapping uses Column F, G, K, L, P, Q, U, V
        # F=5, G=6, K=10, L=11, P=15, Q=16, U=20, V=21
        cols_len = len(df.columns)
        for idx in [5, 6, 10, 11, 15, 16, 20, 21]:
            if idx < cols_len:
                col_letter = chr(65 + idx) # 0=A, 5=F
                col_name = f"Column {col_letter}"
                # Create alias
                df[col_name] = df.iloc[:, idx]
                
    return df

def load_generic_dataset(filepath, dataset_key, file_type=None, strict=False, **kwargs):
    path = pathlib.Path(filepath)
    name = dataset_key
    print(f"Loading {name} from {path}...")
    
    if not path.exists():
        # Fallback for FAYER
        if "FAYER" in name:
            print(f"  File not found at {path}. Searching for alternative...")
            parent = path.parent
            matches = list(parent.glob("*Fayer*Table*.csv")) + list(parent.glob("*Fayer*Table*.xlsx"))
            if matches:
                print(f"  Found alternative: {matches[0]}")
                path = matches[0]
            else:
                print(f"  Error: File not found: {path} and no alternatives found.")
                if strict: raise FileNotFoundError(f"File not found: {path}")
                return pd.DataFrame(columns=["join_key"])
        else:
            print(f"  Error: File not found: {path}")
            if strict: raise FileNotFoundError(f"File not found: {path}")
            return pd.DataFrame(columns=["join_key"])

    # Dataset specific overrides for reading
    read_kwargs = kwargs.copy()
    
    # Header detection
    if name in HEADER_TOKENS_MAP:
        header_row = find_header_row(path, HEADER_TOKENS_MAP[name])
        if header_row is not None:
            print(f"  Detected header at row {header_row}")
            read_kwargs['header'] = header_row
            
    # Fernandes/Caleca defaults
    if name == "CALECA_FILE":
        read_kwargs.setdefault("engine", "python")
        read_kwargs.setdefault("on_bad_lines", "skip")
        read_kwargs.setdefault("index_col", False) 
    
    if name == "FERNANDES_FILE":
        # Maybe latin1?
        pass

    try:
        # Auto-detect file type
        detected_type = file_type
        if detected_type is None:
            if path.suffix.lower() in ['.xlsx', '.xls']:
                detected_type = 'excel'
            elif path.suffix.lower() == '.tsv':
                detected_type = 'tsv'
            else:
                detected_type = 'csv'

        def do_read(f_path, f_type, **kws):
            if f_type == 'excel':
                return pd.read_excel(f_path, **kws)
            elif f_type == 'tsv':
                return pd.read_csv(f_path, sep='\t', **kws)
            else:
                if 'low_memory' not in kws and 'engine' not in kws: 
                    kws['low_memory'] = False
                return pd.read_csv(f_path, **kws)

        try:
            df = do_read(path, detected_type, **read_kwargs)
        except Exception as e:
            # Fallback retry strategies
            print(f"  Read failed ({e}), retrying with latin1 encoding...")
            read_kwargs['encoding'] = 'latin1'
            if detected_type == 'csv':
                read_kwargs['engine'] = 'python' # More robust parser
            df = do_read(path, detected_type, **read_kwargs)
            
        # Strip whitespace from column names
        df.columns = df.columns.astype(str).str.strip()
        
        # Build join key
        df['join_key'] = build_join_key(df, dataset_key)
        
        # Validation
        initial_len = len(df)
        
        # Keep rows with valid keys
        df = df.dropna(subset=['join_key'])
        dropped_null = initial_len - len(df)
        
        if df.duplicated(subset=['join_key']).any():
            dupes = df.duplicated(subset=['join_key']).sum()
            df = df.drop_duplicates(subset=['join_key'])
            print(f"  Dropping {dupes} duplicate keys.")
            
        final_len = len(df)
        print(f"  Loaded {final_len} variants (dropped {dropped_null} null keys).")
        
        # Post-processing (renaming, aliasing)
        df = post_process_dataset(df, dataset_key)
        
        return df
        
    except Exception as e:
        print(f"Error loading {name}: {e}")
        if strict: raise e
        return pd.DataFrame(columns=["join_key"])

def apply_functional_mappings(master_df, mapping_df, dataset_dfs, key_col="join_key"):
    print("Applying functional assay mappings...")
    
    DATASET_ALIASES = {
        "Fayer": "FAYER_FILE",
        "Fayer_2021": "FAYER_FILE",
        "FAYER": "FAYER_FILE",
        "Pillar": "PILLAR_FILE",
        "PILLAR": "PILLAR_FILE",
        "Cravat": "CRAVAT_FILE",
        "CRAVAT": "CRAVAT_FILE",
        "Findlay": "FINDLAY_FILE",
        "Adamovich_Hdr": "ADAMOVICH_HDR_FILE",
        "Adamovich_Cisplatin": "ADAMOVICH_CISPLATIN_FILE",
        "Fernandes": "FERNANDES_FILE",
        "Bouwman_2013": "BOUWMAN_2013_FILE",
        "Bouwman_2020": "BOUWMAN_2020_FILE",
        "Caleca": "CALECA_FILE",
        "Caleca_2019": "CALECA_FILE",
        "Gou": "GOU_FILE",
        "Gou_2023": "GOU_FILE",
        "Bassi": "BASSI_FILE",
        "Bassi_2023": "BASSI_FILE",
        "Langerud": "LANGERUD_FILE",
        "Langerud_2018": "LANGERUD_FILE",
        "Lee_2010_1": "LEE_FILE_1",
        "Lee_2010_2": "LEE_FILE_2",
        "Lee_2010_3": "LEE_FILE_3",
        "Starita": "STARITA_FILE",
        "Hart_2018": "HART_FILE",
        "Richardson_2021": "RICHARDSON_FILE",
        "Hu_2024": "HU_FILE",
        "Huang_2025": "HUANG_FILE",
        "Ikegami_2020": "IKEGAMI_FILE",
        "Hart_2021": "HART_FILE",
        "Biswas_2020": "BISWAS_FILE",
        "Mesman_2021": "MESMAN_FILE",
        "Guidugli_2018": "GUIDUGLI_FILE",
        "Sahu_2023": "SAHU_2023_FILE",
        "Sahu_2025": "SAHU_2025_FILE",
        "Jia": "JIA_FILE",
        "Ollodart": "OLLODART_FILE",
        "Bouvet": "BOUVET_FILE",
        "Matreyek_2018": "MATREYEK_FILE",
        "Mighell_2018": "MIGHELL_FILE",
        "Funk_2025": "FUNK_FILE",
        "Kotler": "KOTLER_FILE",
        "Kawaguchi": "KAWAGUCHI_FILE",
        "Kato": "KATO_FILE",
        "Giacomelli": "GIACOMELLI_FILE",
        "Buckley_2016": "BUCKLEY_FILE"
    }
    
    for _, row in mapping_df.iterrows():
        dataset_raw = str(row["Found in this dataset"]).strip()
        dataset = DATASET_ALIASES.get(dataset_raw, dataset_raw)
        
        src_col = str(row["Mapped to "]).strip()
        dest_col = str(row["Column Name"]).strip()

        if dataset in ["?", "nan", ""] or src_col in ["?", "nan", ""]:
            continue

        df = dataset_dfs.get(dataset)
        
        if df is None:
            # print(f"  Warning: Dataset '{dataset}' (from '{dataset_raw}') not found.")
            continue
            
        if src_col not in df.columns:
            print(f"  Skipping {dataset}.{src_col} -> {dest_col}: missing column")
            # print(f"    Available columns (first 20): {list(df.columns)[:20]}")
            continue

        if dest_col in master_df.columns:
            if master_df[dest_col].notna().any():
                continue
            else:
                master_df = master_df.drop(columns=[dest_col])

        print(f"  Mapping {dataset}.{src_col} -> {dest_col}")
        
        # Support dataset-specific join keys
        current_key_col = key_col
        left_key = key_col
        right_key = key_col # In the loaded df, it is always normalized to 'join_key' column name
        
        if dataset == "LANGERUD_FILE":
             left_key = "HGVSc." # Use HGVSc. column in master_df
             
        if right_key not in df.columns:
            print(f"    Skipping {dataset}: '{right_key}' missing.")
            continue
            
        temp_df = df[[right_key, src_col]].copy()
        temp_df = temp_df.rename(columns={src_col: dest_col})
        
        # Special rename for Langerud to avoid clashing with 'join_key' column in master_df
        # and to merge on the correct column name 'HGVSc.'
        if dataset == "LANGERUD_FILE":
             temp_df = temp_df.rename(columns={right_key: "HGVSc."})
             right_key = "HGVSc."
        
        master_df = master_df.merge(temp_df, left_on=left_key, right_on=right_key, how='left')
        
    return master_df

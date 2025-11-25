import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "MSH2_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "cravat_msh2_cleaned.tsv"
JIA_FILE = INPUT_DIR / "Jia MSH2 2021.xlsx"
OLLODART_FILE = INPUT_DIR / "Ollodart2020.xlsx"
BOUVET_FILE = INPUT_DIR / "MSH2_Bouvet.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Transcript", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "MSH2 Jia auth_func_score", "Jia_auth_reported_functional_class", 
    "MSH2_Ollodart_auth_func_score", "Ollodart_auth_reported_functional_class", 
    "MSH2_Bouvet_auth_func_score", "Bouvet_auth_reported_functional_class", 
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

def load_mave_metadata(filepath):
    """
    Loads MAVE curation metadata to get Interval definitions for MSH2 (Jia dataset).
    Returns a dictionary of interval values.
    """
    print("Loading MAVE metadata...")
    try:
        df = pd.read_csv(filepath, encoding='latin1')
        # Filter for MSH2 and Jia (assuming Jia is primary source for intervals as it appears first in schema)
        # Looking for Dataset_tag containing 'Jia'
        msh2_jia = df[(df['Gene (HGNC symbol)'] == 'MSH2') & (df['Dataset_tag'].str.contains('Jia', na=False))]
        
        if msh2_jia.empty:
            print("Warning: No MSH2 Jia dataset found in MAVE curation.")
            return {}
        
        row = msh2_jia.iloc[0]
        
        metadata = {}
        for i in range(1, 4): # Intervals 1, 2, 3
            metadata[f'Interval {i} name'] = row.get(f'Interval {i} name')
            metadata[f'Interval {i} range'] = row.get(f'Interval {i} range')
            metadata[f'Interval {i} MaveDB class'] = row.get(f'Interval {i} MaveDB class')
            
        # Also get transcript info if available
        metadata['Transcript'] = row.get('Ensembl_transcript_ID') or row.get('Ref_seq_transcript_ID')
        metadata['HGNC ID'] = row.get('HGNC Gene ID')
            
        return metadata
    except Exception as e:
        print(f"Error loading MAVE metadata: {e}")
        return {}

def load_cravat(filepath):
    """
    Loads and cleans CRAVAT output.
    Returns DataFrame with normalized columns and 'join_key' (HGVSp short).
    """
    print(f"Loading CRAVAT file from {filepath}...")
    
    # Columns to keep (mapping)
    usecols = [
        'Variant_Annotation_Gene', 
        'Variant_Annotation_Transcript',
        'Variant_Annotation_cDNA_change', 
        'Variant_Annotation_Protein_Change',
        'Variant_Annotation_Chrom', 
        'Variant_Annotation_Position',
        'Variant_Annotation_End_Position', 
        'Variant_Annotation_Ref_Base', 
        'Variant_Annotation_Alt_Base',
        'gnomAD_Global_AF',
        'ClinVar_Clinical_Significance', 
        'ClinVar_Review_Status',
        'SpliceAI_Acceptor_Gain_Score', 
        'SpliceAI_Acceptor_Loss_Score', 
        'SpliceAI_Donor_Gain_Score', 
        'SpliceAI_Donor_Loss_Score'
    ]
    
    # Check actual columns in file to avoid errors
    df_iter = pd.read_csv(filepath, sep='\t', chunksize=1)
    actual_cols = list(next(df_iter).columns)
    valid_cols = [c for c in usecols if c in actual_cols]
    
    # Read full file (memory permitting, 58k rows is fine)
    df = pd.read_csv(filepath, sep='\t', usecols=valid_cols, low_memory=False)
    
    # Filter for MSH2 gene if needed (though file name implies MSH2)
    if 'Variant_Annotation_Gene' in df.columns:
        df = df[df['Variant_Annotation_Gene'] == 'MSH2'].copy()
        
    # Rename columns
    rename_map = {
        'Variant_Annotation_Gene': 'Gene',
        'Variant_Annotation_Transcript': 'Transcript_Cravat',
        'Variant_Annotation_cDNA_change': 'HGVSc.',
        'Variant_Annotation_Protein_Change': 'HGVSp.',
        'Variant_Annotation_Chrom': 'Chrom',
        'Variant_Annotation_Position': 'hg38_start',
        'Variant_Annotation_End_Position': 'hg38_end',
        'Variant_Annotation_Ref_Base': 'ref_allele',
        'Variant_Annotation_Alt_Base': 'alt_allele',
        'gnomAD_Global_AF': 'gnomad_MAF',
        'ClinVar_Clinical_Significance': 'clinvar_sig_2025',
        'ClinVar_Review_Status': 'clinvar_star_2025_raw',
        'SpliceAI_Acceptor_Gain_Score': 'spliceAI_DS_AG',
        'SpliceAI_Acceptor_Loss_Score': 'spliceAI_DS_AL',
        'SpliceAI_Donor_Gain_Score': 'spliceAI_DS_DG',
        'SpliceAI_Donor_Loss_Score': 'spliceAI_DS_DL'
    }
    df = df.rename(columns=rename_map)
    
    # Map ClinVar stars
    def map_stars(status):
        if pd.isna(status): return np.nan
        status = str(status).lower()
        if 'practice guideline' in status: return 4
        if 'expert panel' in status: return 3
        if 'multiple submitters' in status and 'conflicts' not in status: return 2
        if 'single submitter' in status or 'conflicting' in status: return 1
        if 'no assertion' in status: return 0
        return 0 # Default/Uncertain
        
    if 'clinvar_star_2025_raw' in df.columns:
        df['clinvar_star_2025'] = df['clinvar_star_2025_raw'].apply(map_stars)
    
    # Generate Join Key and AA fields
    print("Parsing HGVSp in CRAVAT data...")
    parsed = df['HGVSp.'].apply(parse_hgvsp)
    
    # parsed is a series of tuples. Expand.
    df['aa_pos'] = parsed.apply(lambda x: x[0])
    df['aa_ref'] = parsed.apply(lambda x: x[1])
    df['aa_alt'] = parsed.apply(lambda x: x[2])
    df['join_key'] = parsed.apply(lambda x: x[3])
    
    # Create a backup join_key from HGVSc if HGVSp is missing (optional, but functional sets are usually protein)
    
    # Handle duplicates in Cravat: 
    # Do NOT drop duplicates on join_key (Protein change) because we want to keep distinct genomic variants (HGVSc).
    # Only drop if strict duplicate rows exist (which read_csv might handle, but we can check).
    if df.duplicated().any():
       df = df.drop_duplicates()

    return df

def load_jia(filepath):
    print(f"Loading Jia data from {filepath}...")
    try:
        df = pd.read_excel(filepath)
        # Variant column seems to be e.g. M1A
        df['join_key'] = df['Variant'].astype(str).str.strip()
        
        # Deduplicate
        if df['join_key'].duplicated().any():
            print(f"Warning: duplicates found in Jia join_key. Dropping.")
            df = df.drop_duplicates(subset=['join_key'])
        
        # Rename columns
        df = df.rename(columns={
            'LOF score': 'MSH2 Jia auth_func_score',
            # Choosing 'Consensus classifiaction' as primary, or 'Clinical classification'? 
            # User prompt says: "Jia_auth_reported_functional_class"
            'Consensus classifiaction': 'Jia_auth_reported_functional_class' 
        })
        
        # Keep only necessary columns for merge
        cols = ['join_key', 'MSH2 Jia auth_func_score', 'Jia_auth_reported_functional_class']
        # Check existence
        cols = [c for c in cols if c in df.columns]
        
        return df[cols]
    except Exception as e:
        print(f"Error loading Jia: {e}")
        return pd.DataFrame()

def load_ollodart(filepath):
    print(f"Loading Ollodart data from {filepath}...")
    try:
        df = pd.read_excel(filepath)
        # Human Msh2 Genotypeb e.g. P27L
        df['join_key'] = df['Human Msh2 Genotypeb'].astype(str).str.strip()
        
        if df['join_key'].duplicated().any():
            print(f"Warning: duplicates found in Ollodart join_key. Dropping.")
            df = df.drop_duplicates(subset=['join_key'])
        
        # Mapping Scoreg to score.
        # Mapping Sigi to class (assuming 'ns' = Normal/neutral, others might be pathogenic?)
        # Let's verify 'Sigi' values later or just copy it.
        # User wants 'auth_reported_functional_class'.
        
        df = df.rename(columns={
            'Scoreg': 'MSH2_Ollodart_auth_func_score',
            'Sigi': 'Ollodart_auth_reported_functional_class' 
        })
        
        cols = ['join_key', 'MSH2_Ollodart_auth_func_score', 'Ollodart_auth_reported_functional_class']
        cols = [c for c in cols if c in df.columns]
        
        return df[cols]
    except Exception as e:
        print(f"Error loading Ollodart: {e}")
        return pd.DataFrame()

def load_bouvet(filepath):
    print(f"Loading Bouvet data from {filepath}...")
    try:
        df = pd.read_csv(filepath)
        # Protein e.g. L85P
        df['join_key'] = df['Protein'].astype(str).str.strip()
        
        if df['join_key'].duplicated().any():
            print(f"Warning: duplicates found in Bouvet join_key. Dropping.")
            df = df.drop_duplicates(subset=['join_key'])
        
        # MT assay result: "77.12% Potentially damaging"
        # We need to extract the score (77.12) and maybe the class if "Potentially damaging" is the class.
        # Or use 'French Classification' / 'InSiGHT Classification'.
        # User wants 'Bouvet_auth_reported_functional_class'. French or InSiGHT seems more "classification" than the raw text.
        # I will use 'French Classification' as it had '3*' in the sample.
        
        # Extract score number from 'MT assay result'
        def extract_score(val):
            if pd.isna(val): return np.nan
            m = re.match(r'([\d\.]+)', str(val))
            if m:
                return float(m.group(1))
            return np.nan
            
        df['MSH2_Bouvet_auth_func_score'] = df['MT assay result'].apply(extract_score)
        
        df = df.rename(columns={
            'French Classification': 'Bouvet_auth_reported_functional_class'
        })
        
        cols = ['join_key', 'MSH2_Bouvet_auth_func_score', 'Bouvet_auth_reported_functional_class']
        cols = [c for c in cols if c in df.columns]
        
        return df[cols]
    except Exception as e:
        print(f"Error loading Bouvet: {e}")
        return pd.DataFrame()

def main():
    print("Starting pipeline...")
    
    # 1. Load MAVE Metadata (Intervals)
    mave_meta = load_mave_metadata(MAVE_FILE)
    
    # 2. Load Cravat (Backbone)
    cravat_df = load_cravat(CRAVAT_FILE)
    print(f"Cravat variants loaded: {len(cravat_df)}")
    
    # 3. Load Functional Datasets
    jia_df = load_jia(JIA_FILE)
    ollodart_df = load_ollodart(OLLODART_FILE)
    bouvet_df = load_bouvet(BOUVET_FILE)
    
    print(f"Jia variants: {len(jia_df)}")
    print(f"Ollodart variants: {len(ollodart_df)}")
    print(f"Bouvet variants: {len(bouvet_df)}")
    
    # 4. Merge
    # Left join onto Cravat. 
    # Note: If cravat doesn't have a variant that is in functional data, it will be lost.
    # Assumption: Cravat file (58k variants) is comprehensive for MSH2.
    
    merged = cravat_df.merge(jia_df, on='join_key', how='left')
    merged = merged.merge(ollodart_df, on='join_key', how='left')
    merged = merged.merge(bouvet_df, on='join_key', how='left')
    
    # 5. Fill Metadata Columns
    # HGNC ID, Transcript from MAVE if available, else from Cravat
    merged['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:7325') # MSH2 HGNC ID
    if 'Transcript_Cravat' in merged.columns:
        # Use Cravat transcript for now, but if MAVE specifies one, we might want to ensure we're using that one or just label it.
        # The target schema just says 'Transcript'. 
        merged['Transcript'] = merged['Transcript_Cravat']
    else:
        merged['Transcript'] = mave_meta.get('Transcript', np.nan)

    # Fill Interval columns (constant values from MAVE metadata)
    for i in range(1, 4):
        merged[f'Interval {i} name'] = mave_meta.get(f'Interval {i} name', np.nan)
        merged[f'Interval {i} range'] = mave_meta.get(f'Interval {i} range', np.nan)
        merged[f'Interval {i} MaveDB class'] = mave_meta.get(f'Interval {i} MaveDB class', np.nan)
        
    # clinvar_date_last_reviewed_2025 -> NaN for now (not in Cravat header)
    merged['clinvar_date_last_reviewed_2025'] = np.nan
    
    # 6. Final Column Selection & Ordering
    # Create missing columns with NaNs
    for col in TARGET_COLUMNS:
        if col not in merged.columns:
            merged[col] = np.nan
            
    final_df = merged[TARGET_COLUMNS]
    
    # 7. Write Output
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    # Optional: Parquet
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

if __name__ == "__main__":
    main()


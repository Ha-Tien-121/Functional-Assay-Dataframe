import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "MSH2_master_dataframe_updated.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "MSH2/MSH2_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "MSH2/MSH2_pillar_data.csv.gz"
JIA_FILE = INPUT_DIR / "MSH2/Jia MSH2 2021.xlsx"
OLLODART_FILE = INPUT_DIR / "MSH2/Ollodart2020.xlsx"
BOUVET_FILE = INPUT_DIR / "MSH2/MSH2_Bouvet_new.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
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
    Now includes ALL columns from Cravat file.
    """
    print(f"Loading CRAVAT file from {filepath}...")
    
    # Load FULL file
    # The file extension is .csv.gz, so we should assume comma separator.
    # The previous code used sep='\t' which caused it to read the entire header as one column.
    df = pd.read_csv(filepath, sep=',', low_memory=False)
    
    # Filter for MSH2 gene if needed (though file name implies MSH2)
    if 'hugo' in df.columns:
        df = df[df['hugo'] == 'MSH2'].copy()
        
    # Rename key columns for consistent merging/schema
    # NOTE: The keys here must match EXACTLY with the CSV header
    rename_map = {
        'hugo': 'Gene',
        'transcript': 'Ensembl_transcript_ID', 
        'extra_vcf_info.HGVS_TRANSCRIPT': 'Ref_seq_transcript_ID', 
        'cchange': 'HGVSc.',
        'achange': 'HGVSp.',
        'chrom': 'Chrom',
        'pos': 'hg38_start',
        'gposend': 'hg38_end',
        'ref_base': 'ref_allele',
        'alt_base': 'alt_allele',
        'gnomad.af': 'gnomad_MAF',
        'clinvar.sig': 'clinvar_sig_2025',
        'clinvar.rev_stat': 'clinvar_star_2025_raw',
        'spliceai.ds_ag': 'spliceAI_DS_AG',
        'spliceai.ds_al': 'spliceAI_DS_AL',
        'spliceai.ds_dg': 'spliceAI_DS_DG',
        'spliceai.ds_dl': 'spliceAI_DS_DL'
    }

    # Strip any potential whitespace from column names before renaming
    df.columns = df.columns.str.strip()
    
    # Perform renaming
    df = df.rename(columns=rename_map)
    
    # Ensure HGVSp. exists
    if 'HGVSp.' not in df.columns:
        print("Error: HGVSp. column missing after rename. Available columns:")
        print(df.columns.tolist())
        # Return empty DF to avoid crashing later, or raise error?
        # Let's raise to fail fast
        raise KeyError("HGVSp. column missing")
            
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
    
    # Handle duplicates in Cravat: 
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
        return pd.DataFrame(columns=['join_key', 'MSH2 Jia auth_func_score', 'Jia_auth_reported_functional_class'])

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
        return pd.DataFrame(columns=['join_key', 'MSH2_Ollodart_auth_func_score', 'Ollodart_auth_reported_functional_class'])

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
        # Extract score and class
        extracted = df['MT assay result'].astype(str).str.extract(r'([\d\.]+)%\s*(.*)')
        
        df['MSH2_Bouvet_auth_func_score'] = pd.to_numeric(extracted[0], errors='coerce')
        df['Bouvet_auth_reported_functional_class'] = extracted[1]
        
        cols = ['join_key', 'MSH2_Bouvet_auth_func_score', 'Bouvet_auth_reported_functional_class']
        cols = [c for c in cols if c in df.columns]
        
        return df[cols]
    except Exception as e:
        print(f"Error loading Bouvet: {e}")
        return pd.DataFrame(columns=['join_key', 'MSH2_Bouvet_auth_func_score', 'Bouvet_auth_reported_functional_class'])

def load_pillar(filepath):
    """
    Loads MSH2 pillar data to get ClinVar 2025 date/stars.
    Returns DataFrame with 'join_key', 'clinvar_date_last_reviewed_2025', 'clinvar_star_2025'.
    (Transcripts now come from Cravat map)
    """
    print(f"Loading Pillar file from {filepath}...")
    try:
        # Load only necessary columns
        usecols = [
            'hgvs_p', 'aa_pos', 'aa_ref', 'aa_alt',
            'clinvar_star_2025', 'clinvar_date_last_reviewed_2025'
        ]
        df = pd.read_csv(filepath, usecols=lambda x: x in usecols, low_memory=False)
        
        # Map ClinVar stars in Pillar data
        def map_stars_pillar(status):
            if pd.isna(status): return np.nan
            status = str(status).lower()
            if 'practice guideline' in status: return 4
            if 'expert panel' in status: return 3
            if 'multiple submitters' in status and 'conflicts' not in status: return 2
            if 'single submitter' in status or 'conflicting' in status: return 1
            if 'no assertion' in status: return 0
            return 0 # Default/Uncertain

        if 'clinvar_star_2025' in df.columns:
            df['clinvar_star_2025'] = df['clinvar_star_2025'].apply(map_stars_pillar)

        # Method 1: Use aa columns if available to construct join_key
        if {'aa_pos', 'aa_ref', 'aa_alt'}.issubset(df.columns):
             df['join_key'] = df.apply(
                 lambda row: f"{row['aa_ref']}{int(row['aa_pos']) if pd.notna(row['aa_pos']) else ''}{row['aa_alt']}" 
                 if pd.notna(row['aa_ref']) and pd.notna(row['aa_pos']) and pd.notna(row['aa_alt']) else np.nan, 
                 axis=1
             )
        
        # Method 2: Fallback to HGVSp
        mask = df.get('join_key', pd.Series([np.nan]*len(df))).isna()
        if mask.any() and 'hgvs_p' in df.columns:
            parsed = df.loc[mask, 'hgvs_p'].apply(parse_hgvsp)
            df.loc[mask, 'join_key'] = parsed.apply(lambda x: x[3])
            
        # Select final columns to merge
        out_cols = [
            'join_key', 
            'clinvar_star_2025', 
            'clinvar_date_last_reviewed_2025'
        ]
        
        df = df[out_cols].dropna(subset=['join_key'])
        
        if df.duplicated(subset=['join_key']).any():
            # If multiple entries, prefer one with date?
            # Sort by date descending (NaNs last) to pick most reviewed/recent
            df = df.sort_values('clinvar_date_last_reviewed_2025', ascending=False)
            df = df.drop_duplicates(subset=['join_key'])
            
        return df
    except Exception as e:
        print(f"Error loading Pillar: {e}")
        return pd.DataFrame(columns=['join_key', 'clinvar_star_2025', 'clinvar_date_last_reviewed_2025'])

def main():
    print("Starting pipeline...")
    
    # 1. Load MAVE Metadata (Intervals)
    mave_meta = load_mave_metadata(MAVE_FILE)
    
    # 2. Load Cravat (Backbone) - NOW INCLUDING ALL COLUMNS
    cravat_df = load_cravat(CRAVAT_FILE)
    print(f"Cravat variants loaded: {len(cravat_df)}")
    
    # 2b. Load Pillar for ClinVar and Transcript updates
    pillar_df = load_pillar(PILLAR_FILE)
    print(f"Pillar variants loaded: {len(pillar_df)}")

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
    
    # Merge Pillar data for ClinVar fields
    # We want to overwrite or fill the clinvar columns from Cravat with Pillar data
    merged = merged.merge(pillar_df, on='join_key', how='left', suffixes=('', '_pillar'))
    
    # Update/Override ClinVar columns
    if 'clinvar_star_2025_pillar' in merged.columns:
        # Use Pillar star if available, else keep existing
        merged['clinvar_star_2025'] = merged['clinvar_star_2025_pillar'].fillna(merged['clinvar_star_2025'])
    
    if 'clinvar_date_last_reviewed_2025_pillar' in merged.columns:
        merged['clinvar_date_last_reviewed_2025'] = merged['clinvar_date_last_reviewed_2025_pillar']
        
    # Transcripts are now mapped in Cravat load step, so no need to map from Pillar here.
    
    # 5. Fill Metadata Columns
    # HGNC ID, Transcript from MAVE if available, else from Cravat
    merged['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:7325') # MSH2 HGNC ID
    
    # Fallback for Transcript IDs if missing in Cravat?
    # Mave metadata might have them
    if 'Ensembl_transcript_ID' not in merged.columns:
        merged['Ensembl_transcript_ID'] = mave_meta.get('Ensembl_transcript_ID', np.nan)
    else:
        merged['Ensembl_transcript_ID'] = merged['Ensembl_transcript_ID'].fillna(mave_meta.get('Ensembl_transcript_ID', np.nan))
        
    if 'Ref_seq_transcript_ID' not in merged.columns:
        merged['Ref_seq_transcript_ID'] = mave_meta.get('Ref_seq_transcript_ID', np.nan)
    else:
        merged['Ref_seq_transcript_ID'] = merged['Ref_seq_transcript_ID'].fillna(mave_meta.get('Ref_seq_transcript_ID', np.nan))

    # Fill Interval columns (constant values from MAVE metadata)
    for i in range(1, 4):
        merged[f'Interval {i} name'] = mave_meta.get(f'Interval {i} name', np.nan)
        merged[f'Interval {i} range'] = mave_meta.get(f'Interval {i} range', np.nan)
        merged[f'Interval {i} MaveDB class'] = mave_meta.get(f'Interval {i} MaveDB class', np.nan)
        
    # 6. Final Column Selection & Ordering
    # We want ALL columns from Cravat + our Target Columns.
    # First, ensure TARGET_COLUMNS exist
    for col in TARGET_COLUMNS:
        if col not in merged.columns:
            merged[col] = np.nan
            
    # Construct final column list:
    # 1. TARGET_COLUMNS (in specific order)
    # 2. All other columns from merged that are NOT in TARGET_COLUMNS (Cravat extra cols)
    
    other_cols = [c for c in merged.columns if c not in TARGET_COLUMNS and not c.endswith('_pillar')]
    final_cols = TARGET_COLUMNS + other_cols
    
    final_df = merged[final_cols]
    
    # 7. Write Output
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    # Optional: Parquet
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

if __name__ == "__main__":
    main()

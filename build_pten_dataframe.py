import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp
from dataset_loader import load_generic_dataset, apply_functional_mappings

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "PTEN_master_dataframe.csv"
MAPPING_FILE = INPUT_DIR / "Functional Assay Mapping - Sheet1.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "PTEN/PTEN_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "PTEN/PTEN_pillar_data.csv.gz"
FAYER_FILE = INPUT_DIR / "PTEN/Fayer et al data.xlsx - Table_S12.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "PTEN_Fayer_2021_Activity_score", "PTEN_Fayer_2021_Activity_class", "PTEN_Fayer_2021_Abundance_score", "PTEN_Fayer_2021_Abundance_class",
    "PTEN_Matreyek_2018_func_score",
    "PTEN_Mighell_2018_func_score",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

def load_mave_metadata(filepath):
    print("Loading MAVE metadata...")
    try:
        df = pd.read_csv(filepath, encoding='latin1')
        gene_data = df[df['Gene (HGNC symbol)'] == 'PTEN']
        
        if gene_data.empty:
            print("Warning: No PTEN dataset found in MAVE curation.")
            return {}
        
        row = gene_data.iloc[0]
        metadata = {}
        for i in range(1, 4):
            metadata[f'Interval {i} name'] = row.get(f'Interval {i} name')
            metadata[f'Interval {i} range'] = row.get(f'Interval {i} range')
            metadata[f'Interval {i} MaveDB class'] = row.get(f'Interval {i} MaveDB class')
            
        metadata['Transcript'] = row.get('Ensembl_transcript_ID') or row.get('Ref_seq_transcript_ID')
        metadata['HGNC ID'] = row.get('HGNC Gene ID')
        return metadata
    except Exception as e:
        print(f"Error loading MAVE metadata: {e}")
        return {}

def load_cravat(filepath):
    print(f"Loading CRAVAT file from {filepath}...")
    df = pd.read_csv(filepath, sep=',', low_memory=False)
    
    if 'hugo' in df.columns:
        df = df[df['hugo'] == 'PTEN'].copy()
        
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
        'clinvar.date_last_reviewed': 'clinvar_date_last_reviewed_2025',
        'spliceai.ds_ag': 'spliceAI_DS_AG',
        'spliceai.ds_al': 'spliceAI_DS_AL',
        'spliceai.ds_dg': 'spliceAI_DS_DG',
        'spliceai.ds_dl': 'spliceAI_DS_DL'
    }

    df.columns = df.columns.str.strip()
    df = df.rename(columns=rename_map)
    
    if 'HGVSp.' not in df.columns:
        raise KeyError("HGVSp. column missing")
            
    def map_stars(status):
        if pd.isna(status): return np.nan
        status = str(status).lower()
        if 'practice guideline' in status: return 4
        if 'expert panel' in status: return 3
        if 'multiple submitters' in status and 'conflicts' not in status: return 2
        if 'single submitter' in status or 'conflicting' in status: return 1
        if 'no assertion' in status: return 0
        return 0 
        
    if 'clinvar_star_2025_raw' in df.columns:
        df['clinvar_star_2025'] = df['clinvar_star_2025_raw'].apply(map_stars)

    print("Parsing HGVSp in CRAVAT data...")
    parsed = df['HGVSp.'].apply(parse_hgvsp)
    
    df['aa_pos'] = parsed.apply(lambda x: x[0])
    df['aa_ref'] = parsed.apply(lambda x: x[1])
    df['aa_alt'] = parsed.apply(lambda x: x[2])
    df['join_key'] = parsed.apply(lambda x: x[3])
    
    if df.duplicated().any():
       df = df.drop_duplicates()

    return df

def load_pillar(filepath, usecols: list[str] | None = None):
    print(f"Loading Pillar file from {filepath}...")
    try:
        df = pd.read_csv(filepath, low_memory=False, usecols=usecols)
        
        def map_stars_pillar(status):
            if pd.isna(status): return np.nan
            status = str(status).lower()
            if 'practice guideline' in status: return 4
            if 'expert panel' in status: return 3
            if 'multiple submitters' in status and 'conflicts' not in status: return 2
            if 'single submitter' in status or 'conflicting' in status: return 1
            if 'no assertion' in status: return 0
            return 0 

        if 'clinvar_star_2025' in df.columns:
            df['clinvar_star_2025'] = df['clinvar_star_2025'].apply(map_stars_pillar)

        if {'aa_pos', 'aa_ref', 'aa_alt'}.issubset(df.columns):
             df['join_key'] = df.apply(
                 lambda row: f"{row['aa_ref']}{int(row['aa_pos']) if pd.notna(row['aa_pos']) else ''}{row['aa_alt']}" 
                 if pd.notna(row['aa_ref']) and pd.notna(row['aa_pos']) and pd.notna(row['aa_alt']) else np.nan, 
                 axis=1
             )
        elif 'hgvs_p' in df.columns:
            parsed = df['hgvs_p'].apply(parse_hgvsp)
            df['join_key'] = parsed.apply(lambda x: x[3])
            
        if 'join_key' in df.columns:
            df = df.dropna(subset=['join_key'])
            if df.duplicated(subset=['join_key']).any():
                df = df.drop_duplicates(subset=['join_key'])
                
        return df
    except Exception as e:
        print(f"Error loading Pillar: {e}")
        return pd.DataFrame()

def main():
    print("Starting PTEN pipeline...")
    
    if not MAPPING_FILE.exists():
        print(f"Error: Mapping file not found at {MAPPING_FILE}")
        return
        
    full_mapping = pd.read_csv(MAPPING_FILE)
    gene_mapping = full_mapping[full_mapping['Gene'] == 'PTEN'].copy()
    
    mave_meta = load_mave_metadata(MAVE_FILE)
    cravat_df = load_cravat(CRAVAT_FILE)
    print(f"Cravat variants loaded: {len(cravat_df)}")
    
    # Identify columns needed from Pillar based on mapping
    pillar_base_cols = {'hgvs_p', 'aa_pos', 'aa_ref', 'aa_alt', 'clinvar_star_2025'}
    pillar_map_rows = gene_mapping[gene_mapping['Found in this dataset'] == 'PILLAR_FILE']
    
    pillar_extra_cols = set()
    for _, row in pillar_map_rows.iterrows():
        if pd.notna(row['Mapped to ']):
             pillar_extra_cols.add(str(row['Mapped to ']).strip())
             
    pillar_usecols = list(pillar_base_cols.union(pillar_extra_cols))
    pillar_df = load_pillar(PILLAR_FILE, usecols=pillar_usecols)
    
    # Load Fayer using generic loader
    fayer_df = load_generic_dataset(FAYER_FILE, "FAYER_FILE", header=1)
    # Fix specific typo in Fayer data
    if 'Abudance_class' in fayer_df.columns:
        fayer_df['Abundance_class'] = fayer_df['Abudance_class']

    dataset_dfs = {
        "CRAVAT_FILE": cravat_df,
        "PILLAR_FILE": pillar_df,
        "FAYER_FILE": fayer_df
    }

    master_df = cravat_df.copy()
    master_df = master_df.merge(pillar_df, on='join_key', how='left', suffixes=('', '_pillar'))
    master_df = apply_functional_mappings(master_df, gene_mapping, dataset_dfs)
    
    # Fix typo from mapping sheet if present (Cleanup step)
    if 'PTEN_Fayer_2021_Abudance_class' in master_df.columns:
        # If the mapping sheet generated the typo column, rename it to the correct one
        if 'PTEN_Fayer_2021_Abundance_class' not in master_df.columns:
             master_df.rename(columns={'PTEN_Fayer_2021_Abudance_class': 'PTEN_Fayer_2021_Abundance_class'}, inplace=True)
        else:
             # Both exist? Copy data if needed and drop typo
             master_df['PTEN_Fayer_2021_Abundance_class'] = master_df['PTEN_Fayer_2021_Abundance_class'].fillna(master_df['PTEN_Fayer_2021_Abudance_class'])
             master_df.drop(columns=['PTEN_Fayer_2021_Abudance_class'], inplace=True)

    master_df['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:9588')
    
    if 'Ensembl_transcript_ID' not in master_df.columns:
        master_df['Ensembl_transcript_ID'] = mave_meta.get('Ensembl_transcript_ID', np.nan)
    else:
        master_df['Ensembl_transcript_ID'] = master_df['Ensembl_transcript_ID'].fillna(mave_meta.get('Ensembl_transcript_ID', np.nan))
        
    if 'Ref_seq_transcript_ID' not in master_df.columns:
        master_df['Ref_seq_transcript_ID'] = mave_meta.get('Ref_seq_transcript_ID', np.nan)
    else:
        master_df['Ref_seq_transcript_ID'] = master_df['Ref_seq_transcript_ID'].fillna(mave_meta.get('Ref_seq_transcript_ID', np.nan))

    for i in range(1, 4):
        master_df[f'Interval {i} name'] = mave_meta.get(f'Interval {i} name', np.nan)
        master_df[f'Interval {i} range'] = mave_meta.get(f'Interval {i} range', np.nan)
        master_df[f'Interval {i} MaveDB class'] = mave_meta.get(f'Interval {i} MaveDB class', np.nan)
        
    for col in TARGET_COLUMNS:
        if col not in master_df.columns:
            master_df[col] = np.nan
            
    other_cols = [c for c in master_df.columns if c not in TARGET_COLUMNS and not c.endswith('_pillar')]
    final_cols = TARGET_COLUMNS + other_cols
    
    final_df = master_df[final_cols]
    
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

    # 9. Write Unmapped Variants
    UNMAPPED_FILE = INPUT_DIR / "PTEN_variants_not_mapped.csv"
    unmapped_mask = final_df['aa_pos'].isna() & final_df['aa_ref'].isna() & final_df['aa_alt'].isna()
    unmapped_df = final_df[unmapped_mask].copy()
    
    if unmapped_df.duplicated().any():
        unmapped_df = unmapped_df.drop_duplicates()
        
    print(f"Writing {len(unmapped_df)} unmapped variants to {UNMAPPED_FILE}...")
    unmapped_df.to_csv(UNMAPPED_FILE, index=False)

    # 10. Write Variants Not Merged
    NOT_MERGED_FILE = INPUT_DIR / "PTEN_variants_not_merged.csv"
    print("Computing variants not merged...")

    master_keys = final_df["join_key"].dropna().drop_duplicates()

    records = []
    for dataset_name, df in dataset_dfs.items():
        if "join_key" not in df.columns:
            continue
        
        tmp = df[["join_key"]].dropna().drop_duplicates().copy()
        tmp["dataset"] = dataset_name
        records.append(tmp)

    if records:
        all_source_keys = pd.concat(records, ignore_index=True)
    else:
        all_source_keys = pd.DataFrame(columns=["join_key", "dataset"])

    if not all_source_keys.empty:
        datasets_per_variant = (
            all_source_keys
            .groupby("join_key")["dataset"]
            .apply(lambda x: ";".join(sorted(x.unique())))
            .reset_index(name="datasets_present")
        )
    else:
        datasets_per_variant = pd.DataFrame(columns=["join_key", "datasets_present"])

    not_merged_mask = ~datasets_per_variant["join_key"].isin(master_keys)
    variants_not_merged = datasets_per_variant[not_merged_mask].copy()

    print("Enriching unmerged variants with source data...")
    unmerged_full = variants_not_merged
    
    for name, df in dataset_dfs.items():
        if 'join_key' in df.columns:
            unmerged_full = unmerged_full.merge(df, on='join_key', how='left', suffixes=('', f'_{name}'))

    unmerged_full['Gene'] = 'PTEN'
    unmerged_full['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:9588')
    
    if 'Ensembl_transcript_ID' not in unmerged_full.columns:
        unmerged_full['Ensembl_transcript_ID'] = np.nan
    unmerged_full['Ensembl_transcript_ID'] = unmerged_full['Ensembl_transcript_ID'].fillna(mave_meta.get('Ensembl_transcript_ID', np.nan))
    
    if 'Ref_seq_transcript_ID' not in unmerged_full.columns:
        unmerged_full['Ref_seq_transcript_ID'] = np.nan
    unmerged_full['Ref_seq_transcript_ID'] = unmerged_full['Ref_seq_transcript_ID'].fillna(mave_meta.get('Ref_seq_transcript_ID', np.nan))

    for i in range(1, 4):
        unmerged_full[f'Interval {i} name'] = mave_meta.get(f'Interval {i} name', np.nan)
        unmerged_full[f'Interval {i} range'] = mave_meta.get(f'Interval {i} range', np.nan)
        unmerged_full[f'Interval {i} MaveDB class'] = mave_meta.get(f'Interval {i} MaveDB class', np.nan)

    def parse_key_to_aa(key):
        match = re.match(r'^([A-Z])(\d+)([A-Z\*])$', str(key))
        if match:
            return match.groups()
        return (np.nan, np.nan, np.nan)

    if 'aa_pos' not in unmerged_full.columns or unmerged_full['aa_pos'].isna().all():
        parsed_aa = unmerged_full['join_key'].apply(parse_key_to_aa)
        unmerged_full['aa_ref'] = parsed_aa.apply(lambda x: x[0])
        unmerged_full['aa_pos'] = parsed_aa.apply(lambda x: x[1])
        unmerged_full['aa_alt'] = parsed_aa.apply(lambda x: x[2])
    
    if 'HGVSp.' not in unmerged_full.columns or unmerged_full['HGVSp.'].isna().all():
         unmerged_full['HGVSp.'] = unmerged_full.apply(
             lambda r: f"p.{r['aa_ref']}{r['aa_pos']}{r['aa_alt']}" 
             if pd.notna(r['aa_pos']) and pd.notna(r['aa_ref']) else np.nan,
             axis=1
         )

    for col in final_df.columns:
        if col not in unmerged_full.columns:
            unmerged_full[col] = np.nan
            
    final_cols_ordered = ['datasets_present'] + list(final_df.columns)
    variants_not_merged_enriched = unmerged_full[final_cols_ordered].copy()

    print(f"Writing {len(variants_not_merged_enriched)} not-merged variants to {NOT_MERGED_FILE}...")
    variants_not_merged_enriched.to_csv(NOT_MERGED_FILE, index=False)

    # Verification Logs
    print("\n--- Dataset Loading Verification ---")
    print(f"{'Dataset':<25} | {'Loaded':<8} | {'Overlap w/ Cravat':<18}")
    print("-" * 55)
    
    # Check Pillar overlap
    pillar_ov = cravat_df.join_key.isin(pillar_df.join_key).mean()
    print(f"{'PILLAR_FILE':<25} | {len(pillar_df):<8} | {pillar_ov:.2%}")
    
    for name, df in dataset_dfs.items():
        if name in ["CRAVAT_FILE", "PILLAR_FILE"]: continue
        if 'join_key' in df.columns:
            ov = cravat_df.join_key.isin(df.join_key).mean()
            print(f"{name:<25} | {len(df):<8} | {ov:.2%}")
        else:
            print(f"{name:<25} | {len(df):<8} | N/A (no key)")

if __name__ == "__main__":
    main()

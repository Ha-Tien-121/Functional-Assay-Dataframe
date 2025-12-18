import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp
from dataset_loader import load_generic_dataset, apply_functional_mappings

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "BRCA1_master_dataframe.csv"
MAPPING_FILE = INPUT_DIR / "Functional Assay Mapping - Sheet1.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "BRCA1/BRCA1_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
FINDLAY_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_HDR_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_CISPLATIN_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
FERNANDES_FILE = INPUT_DIR / "BRCA1/supp_RA118.005274_139996_2_supp_282032_pm8gkd.xlsx - Supp. Table 3.csv"
BOUWMAN_2013_FILE = INPUT_DIR / "BRCA1/supplementary_table_s3_mean_BRCA1_Bouwman_2013_PMID23867111.csv"
BOUWMAN_2020_FILE = INPUT_DIR / "BRCA1/Supplementary Table S2_BRCA1_Bouwman_2020_PMID32546644.xlsx - Sheet1.csv"
CALECA_FILE = INPUT_DIR / "BRCA1/table3_Caleca_BRCA1_BRCA2_2019_PMID30696104.csv"
GOU_FILE = INPUT_DIR / "BRCA1/table1_Guo_BRCA1_2023_PMID37731132.csv"
# FAYER_FILE = INPUT_DIR / "BRCA1/Fayer et al data.xlsx - Table_S7.csv"
BASSI_FILE = INPUT_DIR / "BRCA1/table1_BRCA1_Bassi_2023_PMID37085799.csv"
LANGERUD_FILE = INPUT_DIR / "BRCA1/table2_BRCA1_Langerud_2018_PMID30458859.csv"
LEE_FILE_1 = INPUT_DIR / "BRCA1/Extracted_Supplement_BRCA1_Lee_2010_PMID20516115.csv"
LEE_FILE_2 = INPUT_DIR / "BRCA1/Extracted_Table2_BRCA1_Lee_2010_PMID20516115.csv.csv"
LEE_FILE_3 = INPUT_DIR / "BRCA1/Extracted_Fig3_BRCA1_Lee_2010_PMID20516115.csv"
STARITA_FILE = INPUT_DIR / "BRCA1/TableS7 2_Starita_BRCA1_2018_PMID30219179.tsv"
HART_FILE = INPUT_DIR / "BRCA1/Table S2.xls - Table S2.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns (Must match mapping sheet 'Column Name')
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "BRCA1_Findlay_auth_func_score", "BRCA1_Findlay_reported_functional_class", 
    "BRCA1_Adamovich_HDR_auth_func_score", "BRCA1_Adamovich_HDR_auth_reported_functional_class", 
    "BRCA1_Adamovich_Cisplatin_auth_func_score", "BRCA1_Adamovich_Cisplatin_auth_reported_functional_class", 
    "BRCA1_Fernandes_eta", "BRCA1_Fernandes_fClass_Category",
    "BRCA1_Bouwman_2013_RSE_mean", "BRCA1_Bouwman_2013_IC50_mean", "BRCA1_Bouwman_2013_Pval_mean", "BRCA1_Bouwman_2013_Classification",
    "BRCA1_Bouwman_2020_Cisplatin_nIC50b", "BRCA1_Bouwman_2020_Cisplatin_prediction", "BRCA1_Bouwman_2020_Olaparib_nIC50b", "BRCA1_Bouwman_2020_Olaparib_prediction", "BRCA1_Bouwman_2020_nGFP+b", "BRCA1_Bouwman_2020_DR-GFP_prediction",
    "BRCA1_Caleca_2019_BARD1_Binding", "BRCA1_Caleca_2019_UbCH5a_Binding",
    "BRCA1_Gou_2023_Relative_HR_activity", "BRCA1_Gou_2023_HR_function",
    # "BRCA1_Fayer_2021_function_score", "BRCA1_Fayer_2021_functional_class",
    "BRCA1_Bassi_2023_HDR_assay_His_BRCA1_HeLa_DR_GFP_percent", "BRCA1_Bassi_2023_TA_assay_DBD_BRCT_HEK293_percent", "BRCA1_Bassi_2023_Results_in_functional_assays_this_study",
    "BRCA1_Langerud_HEK293T_TA_Activity_Percent","BRCA1_Langerud_MDA_MB_231_TA_Activity_Percent", "BRCA1_Langerud_2018_Risk_Category",
    "BRCA1_Lee_2010_Structural_Stability_Table_SD", "BRCA1_Lee_2010_Binding_Activity_Table_SD", "BRCA1_Lee_2010_Binding_Specificity_Table","BRCA1_Lee_2010_Transcription_Activity_BarGraph", "BRCA1_Lee_2010_Protease_Sensitivity", "BRCA1_Lee_2010_Binding_Activity", "BRCA1_Lee_2010_Binding_Specificity", "BRCA1_Lee_2010_Transcription_Activity", "BRCA1_Lee_2010_Functional_Effect", 
    "BRCA1_Starita_HDR_predit",
    "BRCA1_Hart_2018_Functional_score", "BRCA1_Hart_2018_Functional_classification",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

def load_mave_metadata(filepath):
    """
    Loads MAVE curation metadata.
    """
    print("Loading MAVE metadata...")
    try:
        df = pd.read_csv(filepath, encoding='latin1')
        # Filter for BRCA1
        gene_data = df[df['Gene (HGNC symbol)'] == 'BRCA1']
        
        if gene_data.empty:
            print("Warning: No BRCA1 dataset found in MAVE curation.")
            return {}
        
        # Take the first available row for intervals (assuming consistent across gene entries or specific one needed)
        # For now, just take the first one
        row = gene_data.iloc[0]
        
        metadata = {}
        for i in range(1, 4): # Intervals 1, 2, 3
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
    """
    Loads and cleans CRAVAT output.
    """
    print(f"Loading CRAVAT file from {filepath}...")
    
    # Load FULL file
    df = pd.read_csv(filepath, sep=',', low_memory=False)
    
    if 'hugo' in df.columns:
        df = df[df['hugo'] == 'BRCA1'].copy()
        
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

def load_pillar(filepath):
    """
    Loads Pillar data (similar to MSH2 logic).
    """
    print(f"Loading Pillar file from {filepath}...")
    try:
        usecols = ['hgvs_p', 'aa_pos', 'aa_ref', 'aa_alt', 'clinvar_star_2025']
        # Load all columns first to avoid error if some missing
        df = pd.read_csv(filepath, low_memory=False)
        
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

        # Join Key
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
    print("Starting BRCA1 pipeline...")
    
    # 0. Load Mapping Configuration
    if not MAPPING_FILE.exists():
        print(f"Error: Mapping file not found at {MAPPING_FILE}")
        return
        
    full_mapping = pd.read_csv(MAPPING_FILE)
    gene_mapping = full_mapping[full_mapping['Gene'] == 'BRCA1'].copy()
    
    # 1. Load MAVE Metadata
    mave_meta = load_mave_metadata(MAVE_FILE)
    
    # 2. Load Cravat (Backbone)
    cravat_df = load_cravat(CRAVAT_FILE)
    print(f"Cravat variants loaded: {len(cravat_df)}")
    
    # 2b. Load Pillar
    pillar_df = load_pillar(PILLAR_FILE)
    
    # 3. Load Functional Datasets
    # Map dataset names from Mapping Sheet to File Paths/Loaders
    # Dataset names in sheet: "Findlay", "Adamovich_Hdr", "Adamovich_Cisplatin", etc.
    
    findlay_df = load_generic_dataset(FINDLAY_FILE, "FINDLAY_FILE")
    adamovich_hdr_df = load_generic_dataset(ADAMOVICH_HDR_FILE, "ADAMOVICH_HDR_FILE")
    adamovich_cisplatin_df = load_generic_dataset(ADAMOVICH_CISPLATIN_FILE, "ADAMOVICH_CISPLATIN_FILE")
    fernandes_df = load_generic_dataset(FERNANDES_FILE, "FERNANDES_FILE")
    bouwman_2013_df = load_generic_dataset(BOUWMAN_2013_FILE, "BOUWMAN_2013_FILE")
    bouwman_2020_df = load_generic_dataset(BOUWMAN_2020_FILE, "BOUWMAN_2020_FILE")
    caleca_df = load_generic_dataset(CALECA_FILE, "CALECA_FILE")
    gou_df = load_generic_dataset(GOU_FILE, "GOU_FILE")
    # fayer_df = load_generic_dataset(FAYER_FILE, "FAYER_FILE")
    bassi_df = load_generic_dataset(BASSI_FILE, "BASSI_FILE")
    langerud_df = load_generic_dataset(LANGERUD_FILE, "LANGERUD_FILE")

    lee_2010_df_1 = load_generic_dataset(LEE_FILE_1, "LEE_FILE_1")
    lee_2010_df_2 = load_generic_dataset(LEE_FILE_2, "LEE_FILE_2")
    lee_2010_df_3 = load_generic_dataset(LEE_FILE_3, "LEE_FILE_3")
    starita_df = load_generic_dataset(STARITA_FILE, "STARITA_FILE", file_type='tsv')
    hart_2018_df = load_generic_dataset(HART_FILE, "HART_FILE")

    dataset_dfs = {
        "CRAVAT_FILE": cravat_df,
        "PILLAR_FILE": pillar_df,
        "FINDLAY_FILE": findlay_df,
        "ADAMOVICH_HDR_FILE": adamovich_hdr_df,
        "ADAMOVICH_CISPLATIN_FILE": adamovich_cisplatin_df,
        "FERNANDES_FILE": fernandes_df,
        "BOUWMAN_2013_FILE": bouwman_2013_df,
        "BOUWMAN_2020_FILE": bouwman_2020_df,
        "CALECA_FILE": caleca_df,
        "GOU_FILE": gou_df,
        # "FAYER_FILE": fayer_df,
        "BASSI_FILE": bassi_df,
        "LANGERUD_FILE": langerud_df,
        "LEE_FILE_1": lee_2010_df_1,
        "LEE_FILE_2": lee_2010_df_2,
        "LEE_FILE_3": lee_2010_df_3,
        "STARITA_FILE": starita_df,
        "HART_FILE": hart_2018_df
    }

    # 4. Initialize Master DataFrame (start with Cravat)
    # We use Cravat as the base 'spine'
    master_df = cravat_df.copy()
    
    # Merge Pillar for ClinVar updates (if applicable)
    master_df = master_df.merge(pillar_df, on='join_key', how='left', suffixes=('', '_pillar'))
    
    # 5. Apply Functional Mappings
    master_df = apply_functional_mappings(master_df, gene_mapping, dataset_dfs)
    
    # 6. Fill Metadata Columns
    master_df['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:1100')
    
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
        
    # 7. Final Column Selection
    for col in TARGET_COLUMNS:
        if col not in master_df.columns:
            master_df[col] = np.nan
            
    # Filter extra columns not in target
    # Keep pillar clinvar if needed? The MSH2 script kept extra cols.
    # "We want ALL columns from Cravat + our Target Columns."
    other_cols = [c for c in master_df.columns if c not in TARGET_COLUMNS and not c.endswith('_pillar')]
    final_cols = TARGET_COLUMNS + other_cols
    
    final_df = master_df[final_cols]
    
    # 8. Write Output
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    # Parquet
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

    # 9. Write Unmapped Variants
    UNMAPPED_FILE = INPUT_DIR / "BRCA1_variants_not_mapped.csv"
    unmapped_mask = final_df['aa_pos'].isna() & final_df['aa_ref'].isna() & final_df['aa_alt'].isna()
    unmapped_df = final_df[unmapped_mask].copy()
    
    if unmapped_df.duplicated().any():
        unmapped_df = unmapped_df.drop_duplicates()
        
    print(f"Writing {len(unmapped_df)} unmapped variants to {UNMAPPED_FILE}...")
    unmapped_df.to_csv(UNMAPPED_FILE, index=False)

    # 10. Write Variants Not Merged
    NOT_MERGED_FILE = INPUT_DIR / "BRCA1_variants_not_merged.csv"
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
    
    # Iteratively merge all source datasets to capture data
    for name, df in dataset_dfs.items():
        if 'join_key' in df.columns:
            unmerged_full = unmerged_full.merge(df, on='join_key', how='left', suffixes=('', f'_{name}'))

    # Metadata
    unmerged_full['Gene'] = 'BRCA1'
    unmerged_full['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:1100')
    
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

    # Parse AA if missing
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

    # Ensure all target columns exist
    for col in final_df.columns:
        if col not in unmerged_full.columns:
            unmerged_full[col] = np.nan
            
    final_cols_ordered = ['datasets_present'] + list(final_df.columns)
    # Filter columns to match master schema + datasets_present
    # Note: This drops extra columns from individual datasets not in master schema
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

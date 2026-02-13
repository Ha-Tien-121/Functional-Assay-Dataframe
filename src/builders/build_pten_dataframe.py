import pandas as pd
import pathlib
import numpy as np
import re
import sys

# Setup paths for imports
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(ROOT_DIR / "src"))

from utils.variant_helpers import parse_hgvsp
from utils.dataset_loader import load_generic_dataset, apply_functional_mappings
from utils.mave_helpers import add_mave_intervals
from utils.splice_utils import compute_spliceai_max, get_spliceai_exclusion_mask, null_out_assay_for_splice

# Configuration
DATA_DIR = ROOT_DIR / "data" / "raw"
REF_DIR = ROOT_DIR / "data" / "reference"
OUTPUT_DIR = ROOT_DIR / "output"

OUTPUT_FILE = OUTPUT_DIR / "master_dataframes" / "PTEN_master_dataframe.csv"
MAPPING_FILE = REF_DIR / "Functional Assay Mapping - Sheet1.csv"

# File paths
CRAVAT_FILE = DATA_DIR / "PTEN/PTEN_annotated.csv.gz"
PILLAR_FILE = DATA_DIR / "PTEN/PTEN_pillar_data.csv.gz"
MAVE_FILE = REF_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "PTEN_Matreyek_2018_func_score",
    "PTEN_Matreyek_2018_reported_functional_class",
    "PTEN_Mighell_2018_func_score",
    "PTEN_Mighell_2018_reported_functional_class",
    "PP_auth_reported_rep_score",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
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

def normalize_pillar_date_column(df):
    """
    Helper function to normalize Pillar ClinVar date column names to clinvar_date_last_reviewed_2025.
    Checks for common alternate column names and renames if found.
    """
    alternate_names = [
        'clinvar_date_last_reviewed_2025',
        'clinvar_date_last_reviewed',
        'clinvar.date_last_reviewed',
        'clinvar_date_last_reviewed_2024',
        'clinvar_date_last_reviewed_2023'
    ]
    
    for alt_name in alternate_names:
        if alt_name in df.columns:
            if alt_name != 'clinvar_date_last_reviewed_2025':
                df = df.rename(columns={alt_name: 'clinvar_date_last_reviewed_2025'})
                print(f"  Renamed Pillar column '{alt_name}' to 'clinvar_date_last_reviewed_2025'")
            break
    
    return df

def load_pillar(filepath, usecols: list[str] | None = None):
    """
    Loads Pillar data for functional mappings and clinvar_date_last_reviewed_2025.
    Note: clinvar_sig_2025 and clinvar_star_2025 come from CRAVAT, but 
    clinvar_date_last_reviewed_2025 comes from Pillar.
    """
    print(f"Loading Pillar file from {filepath}...")
    try:
        # Load with usecols if provided, otherwise load all
        # pandas auto-detects .gz compression
        if usecols:
            df = pd.read_csv(filepath, low_memory=False, usecols=usecols)
        else:
            df = pd.read_csv(filepath, low_memory=False)
        
        # Normalize date column name (check all columns, not just usecols)
        df = normalize_pillar_date_column(df)

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
        
        print(f"  Pillar columns returned: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"Error loading Pillar: {e}")
        return pd.DataFrame()

ALLOWLIST_DATASETS = {
    "FINDLAY_FILE",
    "HUANG_FILE",
    "SAHU_2023_FILE",
    "SAHU_2025_FILE",
    "FUNK_FILE",
    "BUCKLEY_FILE",
}


def main():
    print("Starting PTEN pipeline...")
    
    if not MAPPING_FILE.exists():
        print(f"Error: Mapping file not found at {MAPPING_FILE}")
        return
        
    full_mapping = pd.read_csv(MAPPING_FILE)
    gene_mapping = full_mapping[full_mapping['Gene'] == 'PTEN'].copy()
    
    mave_meta = load_mave_metadata(MAVE_FILE)
    cravat_df = load_cravat(CRAVAT_FILE)
    cravat_df["spliceai_max"] = compute_spliceai_max(cravat_df)
    splice_mask = get_spliceai_exclusion_mask(cravat_df)
    exclude_keys = set(cravat_df.loc[splice_mask, "join_key"].dropna())
    print(f"Cravat variants loaded: {len(cravat_df)}")
    print(f"Variants with spliceai_max > 0.2: {len(exclude_keys)}")
    
    # Identify columns needed from Pillar based on mapping
    # Note: clinvar_sig_2025 and clinvar_star_2025 come from CRAVAT, but clinvar_date_last_reviewed_2025 comes from Pillar
    pillar_base_cols = {'hgvs_p', 'aa_pos', 'aa_ref', 'aa_alt', 'clinvar_date_last_reviewed_2025'}
    pillar_map_rows = gene_mapping[gene_mapping['Found in this dataset'] == 'PILLAR_FILE']
    
    pillar_extra_cols = set()
    for _, row in pillar_map_rows.iterrows():
        if pd.notna(row['Mapped to ']):
             pillar_extra_cols.add(str(row['Mapped to ']).strip())
             
    pillar_usecols = list(pillar_base_cols.union(pillar_extra_cols))
    pillar_df = load_pillar(PILLAR_FILE, usecols=pillar_usecols)
    

    dataset_dfs = {
        "CRAVAT_FILE": cravat_df,
        "PILLAR_FILE": pillar_df,
    }
    for name, df in dataset_dfs.items():
        if name in ALLOWLIST_DATASETS or name in {"CRAVAT_FILE"}:
            continue  # allowlisted assays bypass splice filtering
        nulled = null_out_assay_for_splice(df, exclude_keys)
        if nulled:
            print(f"  SpliceAI filter nulled {nulled} rows for {name}")

    master_df = cravat_df.copy()
    master_df = master_df.merge(pillar_df, on='join_key', how='left', suffixes=('', '_pillar'))

    rep_col = None
    if 'auth_reported_rep_score_pillar' in master_df.columns:
        rep_col = 'auth_reported_rep_score_pillar'
    elif 'auth_reported_rep_score' in master_df.columns:
        rep_col = 'auth_reported_rep_score'
    master_df['PP_auth_reported_rep_score'] = master_df[rep_col] if rep_col else np.nan
    rep_non_null = master_df['PP_auth_reported_rep_score'].notna().sum()
    print(f"  PP_auth_reported_rep_score non-null after Pillar merge: {rep_non_null}")
    
    # Load functional class columns from pillar data by dataset
    print("Loading functional classification columns from Pillar by dataset...")
    pillar_full = pd.read_csv(PILLAR_FILE, low_memory=False)
    
    # Build join_key for pillar_full if not present
    if 'join_key' not in pillar_full.columns:
        if {'aa_pos', 'aa_ref', 'aa_alt'}.issubset(pillar_full.columns):
            pillar_full['join_key'] = pillar_full.apply(
                lambda row: f"{row['aa_ref']}{int(row['aa_pos']) if pd.notna(row['aa_pos']) else ''}{row['aa_alt']}" 
                if pd.notna(row['aa_ref']) and pd.notna(row['aa_pos']) and pd.notna(row['aa_alt']) else np.nan, 
                axis=1
            )
        elif 'hgvs_p' in pillar_full.columns:
            parsed = pillar_full['hgvs_p'].apply(parse_hgvsp)
            pillar_full['join_key'] = parsed.apply(lambda x: x[3])
    
    # Matreyek functional class
    matreyek_df = pillar_full[pillar_full['Dataset'] == 'PTEN_Matreyek_2018'][['join_key', 'auth_reported_func_class']].copy()
    matreyek_df = matreyek_df.rename(columns={'auth_reported_func_class': 'PTEN_Matreyek_2018_reported_functional_class'})
    matreyek_df = matreyek_df.drop_duplicates(subset=['join_key'])
    master_df = master_df.merge(matreyek_df, on='join_key', how='left')
    matreyek_count = master_df['PTEN_Matreyek_2018_reported_functional_class'].notna().sum()
    print(f"  PTEN_Matreyek_2018_reported_functional_class non-null: {matreyek_count}")
    
    # Mighell functional class
    mighell_df = pillar_full[pillar_full['Dataset'] == 'PTEN_Mighell_2018'][['join_key', 'auth_reported_func_class']].copy()
    mighell_df = mighell_df.rename(columns={'auth_reported_func_class': 'PTEN_Mighell_2018_reported_functional_class'})
    mighell_df = mighell_df.drop_duplicates(subset=['join_key'])
    master_df = master_df.merge(mighell_df, on='join_key', how='left')
    mighell_count = master_df['PTEN_Mighell_2018_reported_functional_class'].notna().sum()
    print(f"  PTEN_Mighell_2018_reported_functional_class non-null: {mighell_count}")
    
    # Debug: Print counts before merge update
    cravat_count = master_df['clinvar_date_last_reviewed_2025'].notna().sum() if 'clinvar_date_last_reviewed_2025' in master_df.columns else 0
    pillar_count = master_df['clinvar_date_last_reviewed_2025_pillar'].notna().sum() if 'clinvar_date_last_reviewed_2025_pillar' in master_df.columns else 0
    
    # Update clinvar_date_last_reviewed_2025 from Pillar if available (only fill nulls, don't overwrite non-null)
    if 'clinvar_date_last_reviewed_2025_pillar' in master_df.columns:
        # Prefer CRAVAT if present, otherwise use Pillar
        if 'clinvar_date_last_reviewed_2025' not in master_df.columns:
            master_df['clinvar_date_last_reviewed_2025'] = master_df['clinvar_date_last_reviewed_2025_pillar']
        else:
            # Only fill where CRAVAT value is null
            master_df['clinvar_date_last_reviewed_2025'] = master_df['clinvar_date_last_reviewed_2025'].fillna(master_df['clinvar_date_last_reviewed_2025_pillar'])
        master_df = master_df.drop(columns=['clinvar_date_last_reviewed_2025_pillar'])
    
    # Debug: Print counts after merge update
    merged_count = master_df['clinvar_date_last_reviewed_2025'].notna().sum() if 'clinvar_date_last_reviewed_2025' in master_df.columns else 0
    print(f"  clinvar_date_last_reviewed_2025 counts - CRAVAT: {cravat_count}, Pillar: {pillar_count}, Final merged: {merged_count}")
    
    # Drop other ClinVar columns from Pillar (clinvar_sig and clinvar_star should come from CRAVAT only)
    clinvar_pillar_cols = [c for c in master_df.columns if c.endswith('_pillar') and any(x in c for x in ['clinvar_sig', 'clinvar_star'])]
    if clinvar_pillar_cols:
        master_df = master_df.drop(columns=clinvar_pillar_cols)
    
    master_df = apply_functional_mappings(master_df, gene_mapping, dataset_dfs)

    master_df['HGNC ID'] = mave_meta.get('HGNC ID', 'HGNC:9588')
    
    if 'Ensembl_transcript_ID' not in master_df.columns:
        master_df['Ensembl_transcript_ID'] = mave_meta.get('Ensembl_transcript_ID', np.nan)
    else:
        master_df['Ensembl_transcript_ID'] = master_df['Ensembl_transcript_ID'].fillna(mave_meta.get('Ensembl_transcript_ID', np.nan))
        
    if 'Ref_seq_transcript_ID' not in master_df.columns:
        master_df['Ref_seq_transcript_ID'] = mave_meta.get('Ref_seq_transcript_ID', np.nan)
    else:
        master_df['Ref_seq_transcript_ID'] = master_df['Ref_seq_transcript_ID'].fillna(mave_meta.get('Ref_seq_transcript_ID', np.nan))

    # Attach dataset-specific MAVE interval columns and drop legacy generic ones
    master_df = add_mave_intervals(master_df, MAVE_FILE, "PTEN")

    for col in TARGET_COLUMNS:
        if col not in master_df.columns:
            master_df[col] = np.nan
            
    other_cols = [c for c in master_df.columns if c not in TARGET_COLUMNS and not c.endswith('_pillar')]
    final_cols = TARGET_COLUMNS + other_cols
    
    final_df = master_df[final_cols]
    
    # Reorder columns: Interval columns immediately after clinvar_date_last_reviewed_2025
    anchor = "clinvar_date_last_reviewed_2025"
    interval_cols = [c for c in final_df.columns if c.startswith("Interval ")]
    
    def _interval_sort_key(col):
        parts = col.split()
        i = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 99
        field = parts[2] if len(parts) > 2 else ""
        dataset = " ".join(parts[3:]) if len(parts) > 3 else ""
        field_order = {"name": 0, "range": 1, "MaveDB": 2}
        return (i, field_order.get(field, 9), dataset)
    
    interval_cols = sorted(interval_cols, key=_interval_sort_key)
    cols = list(final_df.columns)
    
    if anchor not in cols:
        raise KeyError(f"{anchor} not found in dataframe columns")
    
    cols_no_intervals = [c for c in cols if c not in interval_cols]
    anchor_idx = cols_no_intervals.index(anchor) + 1
    new_cols = cols_no_intervals[:anchor_idx] + interval_cols + cols_no_intervals[anchor_idx:]
    final_df = final_df[new_cols]
    
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

    # 9. Write Unmapped Variants
    UNMAPPED_FILE = OUTPUT_DIR / "diagnostics" / "PTEN_variants_not_mapped.csv"
    unmapped_mask = final_df['aa_pos'].isna() & final_df['aa_ref'].isna() & final_df['aa_alt'].isna()
    unmapped_df = final_df[unmapped_mask].copy()
    
    if unmapped_df.duplicated().any():
        unmapped_df = unmapped_df.drop_duplicates()
        
    print(f"Writing {len(unmapped_df)} unmapped variants to {UNMAPPED_FILE}...")
    unmapped_df.to_csv(UNMAPPED_FILE, index=False)

    # 10. Write Variants Not Merged
    NOT_MERGED_FILE = OUTPUT_DIR / "diagnostics" / "PTEN_variants_not_merged.csv"
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

    def classify_not_merged_reason(datasets_present):
        ds = {d for d in str(datasets_present).split(";") if d}
        has_cravat = "CRAVAT_FILE" in ds or "Cravat" in ds
        if has_cravat:
            return "present_in_cravat_but_missing_from_master"
        if ds == {"PILLAR_FILE"}:
            return "not_in_cravat_backbone_pillar_only"
        if len(ds) == 1:
            return "not_in_cravat_backbone_single_assay"
        if len(ds) > 1:
            return "not_in_cravat_backbone_multi_assay"
        return "not_in_cravat_backbone_unknown_source"

    variants_not_merged["not_merged_reason"] = variants_not_merged["datasets_present"].apply(classify_not_merged_reason)

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

    # Attach dataset-specific MAVE interval columns and drop legacy generic ones
    unmerged_full = add_mave_intervals(unmerged_full, MAVE_FILE, "PTEN")

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
            
    final_cols_ordered = ['datasets_present', 'not_merged_reason'] + list(final_df.columns)
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

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
from utils.splice_utils import compute_spliceai_max, get_spliceai_exclusion_mask, null_out_assay_for_splice
from utils.dataset_loader import load_generic_dataset, apply_functional_mappings
from utils.mave_helpers import add_mave_intervals

# Configuration
DATA_DIR = ROOT_DIR / "data" / "raw"
REF_DIR = ROOT_DIR / "data" / "reference"
OUTPUT_DIR = ROOT_DIR / "output"

OUTPUT_FILE = OUTPUT_DIR / "master_dataframes" / "BRCA1_master_dataframe.csv"
MAPPING_FILE = REF_DIR / "Functional Assay Mapping - Sheet1.csv"

# File paths
CRAVAT_FILE = DATA_DIR / "BRCA1/BRCA1_annotated.csv.gz"
PILLAR_FILE = DATA_DIR / "BRCA1/BRCA1_pillar_data.csv"
FINDLAY_FILE = DATA_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_HDR_FILE = DATA_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_CISPLATIN_FILE = DATA_DIR / "BRCA1/BRCA1_pillar_data.csv"
FERNANDES_FILE = DATA_DIR / "BRCA1/supp_RA118.005274_139996_2_supp_282032_pm8gkd.xlsx - Supp. Table 3.csv"
BOUWMAN_2013_FILE = DATA_DIR / "BRCA1/supplementary_table_s3_mean_BRCA1_Bouwman_2013_PMID23867111.csv"
BOUWMAN_2020_FILE = DATA_DIR / "BRCA1/Supplementary Table S2_BRCA1_Bouwman_2020_PMID32546644.xlsx - Sheet1.csv"
CALECA_FILE = DATA_DIR / "BRCA1/table3_Caleca_BRCA1_BRCA2_2019_PMID30696104.csv"
GOU_FILE1 = DATA_DIR / "BRCA1/table1_Guo_BRCA1_2023_PMID37731132.csv"
GOU_FILE2 = DATA_DIR / "BRCA1/supplemental_table_S1_Guo_BRCA1_BRCA2_control_variants_2023_PMID37731132.csv"
# FAYER_FILE = DATA_DIR / "BRCA1/Fayer et al data.xlsx - Table_S7.csv"
BASSI_FILE = DATA_DIR / "BRCA1/table1_BRCA1_Bassi_2023_PMID37085799.csv"
LANGERUD_FILE = DATA_DIR / "BRCA1/table2_BRCA1_Langerud_2018_PMID30458859.csv"
LEE_FILE_1 = DATA_DIR / "BRCA1/Extracted_Supplement_BRCA1_Lee_2010_PMID20516115.csv"
LEE_FILE_2 = DATA_DIR / "BRCA1/Extracted_Table2_BRCA1_Lee_2010_PMID20516115.csv.csv"
LEE_FILE_3 = DATA_DIR / "BRCA1/Extracted_Fig3_BRCA1_Lee_2010_PMID20516115.csv"
STARITA_FILE1 = DATA_DIR / "BRCA1/TableS7 2_Starita_BRCA1_2018_PMID30219179.tsv"
STARITA_FILE2 = DATA_DIR / "BRCA1/urn_mavedb_00000081-a-1_scores.csv"
STARITA_FILE3 = DATA_DIR / "BRCA1/urn_mavedb_00000081-a-2_scores.csv"
HART_FILE = DATA_DIR / "BRCA1/Table S2.xls - Table S2.csv"
MAVE_FILE = REF_DIR / "MAVE Curation v3.csv"

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
    "BRCA1_Starita_HDR_predit", "BRCA1_Starita_2018HDR_depletion_score", "BRCA1_Starita_2018HDR_fluorescence_score",
    "BRCA1_Hart_2018_Functional_score", "BRCA1_Hart_2018_Functional_classification",
    "PP_auth_reported_rep_score",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
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
    Loads Pillar data for functional mappings.
    """
    print(f"Loading Pillar file from {filepath}...")
    try:
        # Load all columns (usecols not used, but kept for reference)
        # Note: ClinVar columns will be dropped after merge
        df = pd.read_csv(filepath, low_memory=False)

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

ALLOWLIST_DATASETS = {
    "FINDLAY_FILE",  # BRCA1_Findlay_2018 exempt from SpliceAI exclusion
    "HUANG_FILE",    # Huang 2025
    "SAHU_2023_FILE",
    "SAHU_2025_FILE",
    "FUNK_FILE",
    "BUCKLEY_FILE",
}


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
    cravat_df["spliceai_max"] = compute_spliceai_max(cravat_df)
    splice_mask = get_spliceai_exclusion_mask(cravat_df)
    exclude_keys = set(cravat_df.loc[splice_mask, "join_key"].dropna())
    # LANGERUD uses HGVSc. as join key, so build a separate set for it
    exclude_keys_hgvsc = set(cravat_df.loc[splice_mask, "HGVSc."].dropna())
    print(f"Cravat variants loaded: {len(cravat_df)}")
    print(f"Variants with spliceai_max > 0.2: {len(exclude_keys)}")
    
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
    gou_df = load_generic_dataset(GOU_FILE1, "GOU_FILE")
    gou_df2 = load_generic_dataset(GOU_FILE2, "GOU_FILE2")
    if gou_df2 is not None and not gou_df2.empty:
        if 'Gene' in gou_df2.columns:
            gou_df2 = gou_df2[gou_df2['Gene'].astype(str).str.upper() == 'BRCA1']
        rename_map = {}
        if 'Relative_HR_activity' in gou_df2.columns:
            rename_map['Relative_HR_activity'] = 'BRCA1_Gou_2023_Relative_HR_activity'
        if 'Clinical_significance' in gou_df2.columns:
            rename_map['Clinical_significance'] = 'BRCA1_Gou_2023_HR_function'
        if rename_map:
            gou_df2 = gou_df2.rename(columns=rename_map)
        keep_cols = ['join_key', 'BRCA1_Gou_2023_Relative_HR_activity', 'BRCA1_Gou_2023_HR_function']
        gou_df2 = gou_df2[[c for c in keep_cols if c in gou_df2.columns]]
        if 'join_key' in gou_df2.columns:
            gou_df2 = gou_df2.dropna(subset=['join_key']).drop_duplicates(subset=['join_key'])
    # fayer_df = load_generic_dataset(FAYER_FILE, "FAYER_FILE")
    bassi_df = load_generic_dataset(BASSI_FILE, "BASSI_FILE")
    langerud_df = load_generic_dataset(LANGERUD_FILE, "LANGERUD_FILE")

    lee_2010_df_1 = load_generic_dataset(LEE_FILE_1, "LEE_FILE_1")
    lee_2010_df_2 = load_generic_dataset(LEE_FILE_2, "LEE_FILE_2")
    lee_2010_df_3 = load_generic_dataset(LEE_FILE_3, "LEE_FILE_3")
    starita_df = load_generic_dataset(STARITA_FILE1, "STARITA_FILE", file_type='tsv')
    starita2_df = load_generic_dataset(STARITA_FILE2, "STARITA_FILE2")
    starita3_df = load_generic_dataset(STARITA_FILE3, "STARITA_FILE3")

    # Filter Starita to only Pass/Pass on E3 and Y2H before any renaming or keying
    if starita_df is not None and not starita_df.empty:
        # Prefer the pass/fail classification columns if present, else fall back to raw score columns
        e3_col = "Starita_E3_800_filter" if "Starita_E3_800_filter" in starita_df.columns else "Starita_E3_score"
        y2h_col = "Starita_Y2H_800_filter" if "Starita_Y2H_800_filter" in starita_df.columns else "Starita_Y2H_score"

        required_cols = [e3_col, y2h_col]
        missing = [c for c in required_cols if c not in starita_df.columns]
        if missing:
            raise KeyError(f"Starita filter requires columns {required_cols}, missing: {missing}")

        pre_count = len(starita_df)
        starita_df = starita_df[
            starita_df[e3_col].astype(str).str.strip().str.lower().eq("pass")
            & starita_df[y2h_col].astype(str).str.strip().str.lower().eq("pass")
        ]
        post_count = len(starita_df)
        print(f"Starita variants before filter: {pre_count}; after Pass/Pass filter: {post_count}")
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
        "GOU_FILE2": gou_df2,
        # "FAYER_FILE": fayer_df,
        "BASSI_FILE": bassi_df,
        "LANGERUD_FILE": langerud_df,
        "LEE_FILE_1": lee_2010_df_1,
        "LEE_FILE_2": lee_2010_df_2,
        "LEE_FILE_3": lee_2010_df_3,
        "STARITA_FILE": starita_df,
        "STARITA_FILE2": starita2_df,
        "STARITA_FILE3": starita3_df,
        "HART_FILE": hart_2018_df
    }

    # 4. Initialize Master DataFrame (start with Cravat)
    # We use Cravat as the base 'spine'
    master_df = cravat_df.copy()
    
    # Merge Pillar (for functional mappings, not ClinVar fields)
    master_df = master_df.merge(pillar_df, on='join_key', how='left', suffixes=('', '_pillar'))

    # Add replicate score from Pillar as PP_auth_reported_rep_score
    rep_col = None
    if 'auth_reported_rep_score_pillar' in master_df.columns:
        rep_col = 'auth_reported_rep_score_pillar'
    elif 'auth_reported_rep_score' in master_df.columns:
        rep_col = 'auth_reported_rep_score'

    master_df['PP_auth_reported_rep_score'] = master_df[rep_col] if rep_col else np.nan
    rep_non_null = master_df['PP_auth_reported_rep_score'].notna().sum()
    print(f"  PP_auth_reported_rep_score non-null after Pillar merge: {rep_non_null}")
    
    clinvar_pillar_cols = [c for c in master_df.columns if c.endswith('_pillar') and any(x in c for x in ['clinvar_sig', 'clinvar_star', 'clinvar_date'])]
    if clinvar_pillar_cols:
        master_df = master_df.drop(columns=clinvar_pillar_cols)
    
    # 5. Apply Functional Mappings
    # Apply SpliceAI exclusion per assay (allowlist bypassed)
    for name, df in dataset_dfs.items():
        if name in ALLOWLIST_DATASETS or name in {"CRAVAT_FILE"}:
            # Allowlisted assays bypass SpliceAI filtering
            continue
        # LANGERUD uses HGVSc. values in its join_key column, use exclude_keys_hgvsc
        keys_to_use = exclude_keys_hgvsc if name == "LANGERUD_FILE" else exclude_keys
        nulled = null_out_assay_for_splice(df, keys_to_use)
        if nulled:
            print(f"  SpliceAI filter nulled {nulled} rows for {name}")

    master_df = apply_functional_mappings(master_df, gene_mapping, dataset_dfs)

    # Merge Gou file 2 (control variants) without overwriting existing values
    if gou_df2 is not None and not gou_df2.empty:
        master_df = master_df.merge(gou_df2, on='join_key', how='left', suffixes=('', '_gou2_tmp'))
        for col in ['BRCA1_Gou_2023_Relative_HR_activity', 'BRCA1_Gou_2023_HR_function']:
            tmp = f"{col}_gou2_tmp"
            if tmp in master_df.columns:
                master_df[col] = master_df[col].fillna(master_df[tmp])
                master_df = master_df.drop(columns=[tmp])

    # 5b. Merge Starita MaveDB scores (new files)
    def _merge_starita(df, target_col):
        if df is None or df.empty:
            return
        if 'join_key' not in df.columns:
            return
        df = df.copy()
        df[target_col] = pd.to_numeric(df.get('score'), errors='coerce')
        df = df[['join_key', target_col]].dropna(subset=['join_key'])
        # Drop duplicates keeping first
        df = df.drop_duplicates(subset=['join_key'])
        non_null = df[target_col].notna().sum()
        print(f"  Merged {non_null} rows into {target_col} from Starita file")
        return df
        return df

    starita2_merge = _merge_starita(starita2_df, 'BRCA1_Starita_2018HDR_fluorescence_score')
    starita3_merge = _merge_starita(starita3_df, 'BRCA1_Starita_2018HDR_depletion_score')
    for merged_df in [starita2_merge, starita3_merge]:
        if merged_df is not None:
            master_df = master_df.merge(merged_df, on='join_key', how='left')
    
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

    # Attach dataset-specific MAVE interval columns and drop legacy generic ones
    master_df = add_mave_intervals(master_df, MAVE_FILE, "BRCA1")

    # Drop metadata columns that appear after join_key
    drop_after_join_key = [
        "ID", "Dataset", "HGNC_id", "STRAND", "hg19_pos",
        "auth_transcript_id", "transcript_pos", "transcript_ref", "transcript_alt",
        "hgvs_c", "hgvs_p", "consequence", "simplified_consequence",
        "auth_reported_score", "auth_reported_rep_score", "auth_reported_func_class",
        "splice_measure", "nucleotide_or_aa",
        "MaveDB Score Set URN", "Model_system", "Assay Type",
        "Phenotype Measured ontology term",
        "Molecular or Biological Process Investigated (GO term)",
        "IGVF_produced", "Flag",
        "REVEL", "AM_score", "AM_class",
        "spliceAI_DP_AG", "spliceAI_DP_AL", "spliceAI_DP_DG", "spliceAI_DP_DL",
        "REVEL_train", "MutPred2", "MP2_train",
        "clinvar_sig_2018", "clinvar_star_2018", "clinvar_date_last_reviewed_2018",
        "ClinVar Variation Id_ClinGen_repo", "Allele Registry Id_ClinGen_repo",
        "Disease_ClinGen_repo", "Mondo Id_ClinGen_repo",
        "Mode of Inheritance_ClinGen_repo", "Assertion_ClinGen_repo",
        "Applied Evidence Codes (Met)_ClinGen_repo",
        "Applied Evidence Codes (Not Met)_ClinGen_repo",
        "Summary of interpretation_ClinGen_repo", "PubMed Articles_ClinGen_repo",
        "Expert Panel_ClinGen_repo", "Guideline_ClinGen_repo",
        "Approval Date_ClinGen_repo", "Published Date_ClinGen_repo",
        "Retracted_ClinGen_repo", "Evidence Repo Link_ClinGen_repo",
        "Uuid_ClinGen_repo", "Updated_Classification_ClinGen_repo",
        "Updated_Evidence Codes_ClinGen_repo"
    ]
    
    if "join_key" in master_df.columns:
        join_idx = list(master_df.columns).index("join_key")
        cols_to_drop = [
            c for c in drop_after_join_key
            if c in master_df.columns and list(master_df.columns).index(c) > join_idx
        ]
        if cols_to_drop:
            master_df = master_df.drop(columns=cols_to_drop)
            print(f"  Dropped {len(cols_to_drop)} metadata columns appearing after 'join_key': {cols_to_drop[:5]}{'...' if len(cols_to_drop) > 5 else ''}")

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
    
    # Reorder columns: Interval columns immediately after clinvar_date_last_reviewed_2025
    anchor = "clinvar_date_last_reviewed_2025"
    interval_cols = [c for c in final_df.columns if c.startswith("Interval ")]
    
    # Optional: deterministic ordering
    def _interval_sort_key(col):
        parts = col.split()
        # "Interval {i} {field} {dataset...}"
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
    
    # 8. Write Output
    print(f"Writing {len(final_df)} rows to {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    # Parquet
    parquet_file = OUTPUT_FILE.with_suffix('.parquet')
    final_df.to_parquet(parquet_file, index=False)
    print(f"Written parquet to {parquet_file}")

    # 9. Write Unmapped Variants
    UNMAPPED_FILE = OUTPUT_DIR / "diagnostics" / "BRCA1_variants_not_mapped.csv"
    unmapped_mask = final_df['aa_pos'].isna() & final_df['aa_ref'].isna() & final_df['aa_alt'].isna()
    unmapped_df = final_df[unmapped_mask].copy()
    
    if unmapped_df.duplicated().any():
        unmapped_df = unmapped_df.drop_duplicates()
        
    print(f"Writing {len(unmapped_df)} unmapped variants to {UNMAPPED_FILE}...")
    unmapped_df.to_csv(UNMAPPED_FILE, index=False)

    # 10. Write Variants Not Merged
    NOT_MERGED_FILE = OUTPUT_DIR / "diagnostics" / "BRCA1_variants_not_merged.csv"
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

    # Attach dataset-specific MAVE interval columns and drop legacy generic ones
    unmerged_full = add_mave_intervals(unmerged_full, MAVE_FILE, "BRCA1")

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
            
    final_cols_ordered = ['datasets_present', 'not_merged_reason'] + list(final_df.columns)
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

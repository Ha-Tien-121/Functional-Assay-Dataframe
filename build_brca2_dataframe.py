import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "BRCA2_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "BRCA2/BRCA2_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "BRCA2/BRCA2_pillar_data.csv"
RICHARDSON_FILE = INPUT_DIR / "BRCA2/Supplemental_Table_1_Richardson_BRCA2_2021_PMID33609447.xlsx - Sheet1.csv"
HU_FILE = INPUT_DIR / "BRCA2/BRCA2_pillar_data.csv"
HUANG_FILE = INPUT_DIR / "BRCA2/SuppTables_Huang_BRCA2_2025_PMID39779857.xlsx - Table S1.csv"
IKEGAMI_FILE = INPUT_DIR / "BRCA2/SuppData5_Ikegami_BRCA2_2020_PMID32444794.xlsx - Results of Bayesian inference.csv"
HART_FILE = INPUT_DIR / "BRCA2/Table S2.xls - Table S2.csv"
BISWAS_FILE = INPUT_DIR / "BRCA2/Table 1 extracted - Sheet1.csv"
MESMAN_FILE = INPUT_DIR / "BRCA2/Table 1 and Table 2 extracted - Table 1.csv"
GUIDUGLI_FILE = INPUT_DIR / "BRCA2/table1_BRCA2_Guidugli_2018_PMID29394989.csv"
SAHU_2023_FILE = INPUT_DIR / "BRCA2/BRCA2_pillar_data.csv"
SAHU_2025_FILE = INPUT_DIR / "BRCA2/BRCA2_pillar_data.csv"
CALECA_FILE = INPUT_DIR / "BRCA2/table4_Caleca_BRCA1_BRCA2_2019_PMID30696104.csv"
GOU_FILE = INPUT_DIR / "BRCA2/table2_Guo_BRCA2_2023_PMID37731132.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "BRCA2_Richardson_2021_HDR_score", "BRCA2_Richardson_2021_HDR_Function",
    "BRCA2_Hu_2024_func_score", "BRCA2_Hu_2024_HDR_func_class",
    "BRCA2_Huang_2025_func_score", "BRCA2_Huang_2025_HDR_func_category",
    "BRCA2_Ikegami_2020_Olaparib_BF", "BRCA2_Ikegami_2020_Niraparib_BF", "BRCA2_Ikegami_2020_Rucaparib_BF", "BRCA2_Ikegami_2020_CBDCA_BF", "BRCA2_Ikegami_2020_Olaparib_fClass", "BRCA2_Ikegami_2020_Niraparib_fClass", "BRCA2_Ikegami_2020_Rucaparib_fClass", "BRCA2_Ikegami_2020_CBDCA_fClass",
    "BRCA2_Hart_2021_func_score", "BRCA2_Hart_2021_func_class",
    "BRCA2_Biswas_2020_PIF[HAT+DS]", "BRCA2_Biswas_2020_func_class",
    "BRCA2_Mesman_2021_Complementation", "BRCA2_Mesman_2021_HDR_capacity", "BRCA2_Mesman_2021_Cisplatin_sensitivity", 
    "BRCA2_Guidugli_2018_HDR_FC", "BRCA2_Guidugli_2018_HDR_annotation",
    "BRCA2_Sahu_2023_function_score", "BRCA2_Sahu_2023_functional_class",
    "BRCA2_Sahu_2025_function_score", "BRCA2_Sahu_2025_functional_class",
    "BRCA2_Caleca_2019_BRCA2_DSS1_Binding",
    "BRCA2_Gou_2023_Relative_HR_activity", "BRCA2_Gou_2023_HR_function",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

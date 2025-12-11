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
BISWAS_FILE = INPUT_DIR / "BRCA2/supplemental_table1_Biswas_2020_BRCA2_PMID33293522.csv"
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
    # "MSH2 Jia auth_func_score", "Jia_auth_reported_functional_class", 
    # "MSH2_Ollodart_auth_func_score", "Ollodart_auth_reported_functional_class", 
    # "MSH2_Bouvet_auth_func_score", "Bouvet_auth_reported_functional_class", 
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

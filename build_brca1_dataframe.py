import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "BRCA1_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "BRCA1/BRCA1_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
FINDLAY_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_HDR_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
ADAMOVICH_CISPLATIN_FILE = INPUT_DIR / "BRCA1/BRCA1_pillar_data.csv"
FERNANDES_FILE = INPUT_DIR / "BRCA1/supp_RA118.005274_139996_2_supp_282032_pm8gkd.xlsx - Supp. Table 1.csv"
BOUWMAN_2013_FILE = INPUT_DIR / "BRCA1/supplementary_table_s3_mean_BRCA1_Bouwman_2013_PMID23867111.csv"
BOUWMAN_2020_FILE = INPUT_DIR / "BRCA1/Supplementary Table S2_BRCA1_Bouwman_2020_PMID32546644.xlsx - Sheet1.csv"
CALECA_FILE = INPUT_DIR / "BRCA1/table3_Caleca_BRCA1_BRCA2_2019_PMID30696104.csv"
GOU_FILE = INPUT_DIR / "BRCA1/table1_Guo_BRCA1_2023_PMID37731132.csv"
FAYER_FILE = INPUT_DIR / "BRCA1/Fayer et al data.xlsx - Table_S1.csv"
BASSI_FILE = INPUT_DIR / "BRCA1/table1_BRCA1_Bassi_2023_PMID37085799.csv"
LANGERUD_FILE = INPUT_DIR / "BRCA1/table2_BRCA1_Langerud_2018_PMID30458859.csv"
LEE_FILE_1 = INPUT_DIR / "BRCA1/Extracted_Supplement_BRCA1_Lee_2010_PMID20516115.csv"
LEE_FILE_2 = INPUT_DIR / "BRCA1/Extracted_Table2_BRCA1_Lee_2010_PMID20516115.csv.csv"
LEE_FILE_3 = INPUT_DIR / "BRCA1/Extracted_Fig3_BRCA1_Lee_2010_PMID20516115.csv"
STARITA_FILE = INPUT_DIR / "BRCA1/TableS7 2_Starita_BRCA1_2018_PMID30219179.tsv"
HART_FILE = INPUT_DIR / "BRCA1/Table S2.xls - Table S2.csv"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"
# PETITALOT_FILE = INPUT_DIR / ""

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

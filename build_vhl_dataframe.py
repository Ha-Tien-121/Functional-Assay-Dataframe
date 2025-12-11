import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "VHL_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "VHL\VHL_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "VHL\VHL_pillar_data.csv"
BUCKLEY_FILE = INPUT_DIR / "VHL\VHL_pillar_data.csv"

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

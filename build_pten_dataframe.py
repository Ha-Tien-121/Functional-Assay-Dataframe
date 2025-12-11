import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "PTEN_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "PTEN/PTEN_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "PTEN/PTEN_pillar_data.csv.gz"
FAYER_FILE = INPUT_DIR / "PTEN/Fayer et al data.xlsx - Table_S1.csv"
MATREYEK_FILE = INPUT_DIR / "PTEN/PTEN_pillar_data.csv.gz"
MIGHELL_FILE = INPUT_DIR / "PTEN/PTEN_pillar_data.csv.gz"
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "PTEN_Fayer_2021_Activity_score", "PTEN_Fayer_2021_Activity_class", "PTEN_Fayer_2021_Abundance_score", "PTEN_Fayer_2021_Abudance_class",
    "PTEN_Matreyek_2018_func_score",
    "PTEN_Mighell_2018_func_score",
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

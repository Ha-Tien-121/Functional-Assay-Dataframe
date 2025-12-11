import pandas as pd
import pathlib
import numpy as np
import re
from variant_helpers import parse_hgvsp

# Configuration
INPUT_DIR = pathlib.Path(".")
OUTPUT_FILE = INPUT_DIR / "TP53_master_dataframe.csv"

# File paths
CRAVAT_FILE = INPUT_DIR / "TP53/TP53_annotated.csv.gz"
PILLAR_FILE = INPUT_DIR / "TP53/TP53_pillar_data.csv.gz"
FAYER_FILE = INPUT_DIR / "TP53/Fayer et al data.xlsx - Table_S1.csv"
FUNK_FILE = INPUT_DIR / "TP53/Funk et al. Supplementary tables.xlsx - Supp_Table_1.csv"
KOTLER_FILE = INPUT_DIR / "TP53/Kotler et al Supplemental table.xlsx - H1299_RFS_sequence_variants.csv"
KATO_FILE = INPUT_DIR / "TP53/TP53_pillar_data.csv.gz"
GIACOMELLI_FILE = INPUT_DIR / "TP53/TP53_pillar_data.csv.gz"
# KAWAGUCHI_FILE = INPUT_DIR / ""
MAVE_FILE = INPUT_DIR / "MAVE Curation v3.csv"

# Target columns
TARGET_COLUMNS = [
    "Gene", "HGNC ID", "Ensembl_transcript_ID", "Ref_seq_transcript_ID", "HGVSc.", "HGVSp.", "Chrom", 
    "hg38_start", "hg38_end", "ref_allele", "alt_allele", 
    "aa_pos", "aa_ref", "aa_alt", 
    "TP53_Fayer_2021_Classifier_Prob_func_normal", "TP53_Fayer_2021_Classifier_Prob_func_abnormal", "TP53_Fayer_2021_Classifier_prediction",
    "TP53_Funk_2025_rfs_median", "TP53_Funk_2025_func_class",
    "TP53_Kotler_RFS_H1299", "TP53_Kotler_func_class",
    "TP53_Kawaguchi_func_class",
    "TP53_Kato_func_score", "TP53_Kato_func_class",
    "TP53_Giacomelli_func_score", "TP53_Giacomelli_func_class", 
    "gnomad_MAF", 
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025", 
    "Interval 1 name", "Interval 1 range", "Interval 1 MaveDB class", 
    "Interval 2 name", "Interval 2 range", "Interval 2 MaveDB class", 
    "Interval 3 name", "Interval 3 range", "Interval 3 MaveDB class", 
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"
]

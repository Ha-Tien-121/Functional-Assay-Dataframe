# Functional Assay Dataframe

A pipeline for building comprehensive master dataframes of functional assay data for hereditary cancer gene variants, with likelihood ratio (LR) calibration for clinical interpretation.

## Supported Genes

- **BRCA1** - Breast Cancer 1
- **BRCA2** - Breast Cancer 2  
- **MSH2** - MutS Homolog 2
- **PTEN** - Phosphatase and Tensin Homolog
- **TP53** - Tumor Protein P53
- **VHL** - Von Hippel-Lindau

## Project Structure

```
Functional-Assay-Dataframe/
│
├── src/                          # Python source code
│   ├── builders/                 # Gene-specific dataframe builders
│   │   ├── build_brca1_dataframe.py
│   │   ├── build_brca2_dataframe.py
│   │   ├── build_msh2_dataframe.py
│   │   ├── build_pten_dataframe.py
│   │   ├── build_tp53_dataframe.py
│   │   └── build_vhl_dataframe.py
│   │
│   ├── calibration/              # Evidence strength calibration
│   │   ├── evidence_from_clinvar.py           # Original calibration
│   │   ├── evidence_from_clinvar_with_changes.py  # Hybrid calibration
│   │   └── run_all_genes.py      # Batch calibration runner
│   │
│   ├── utils/                    # Shared helper modules
│   │   ├── dataset_loader.py     # Data loading utilities
│   │   ├── mave_helpers.py       # MaveDB threshold parsing
│   │   ├── splice_utils.py       # SpliceAI filtering
│   │   └── variant_helpers.py    # Variant parsing utilities
│   │
│   └── scripts/                  # Utility scripts
│       ├── gzip_output_files.py  # Compress output files
│       └── ...
│
├── data/                         # Input data
│   ├── raw/                      # Gene-specific source data
│   │   ├── BRCA1/
│   │   ├── BRCA2/
│   │   ├── MSH2/
│   │   ├── PTEN/
│   │   ├── TP53/
│   │   └── VHL/
│   │
│   └── reference/                # Reference/mapping files
│       ├── MAVE Curation v3.csv
│       ├── Functional Assay Mapping - Sheet1.csv
│       └── final_pillar_data_with_clinvar_*.csv.gz
│
├── output/                       # Generated outputs
│   ├── master_dataframes/        # Main consolidated dataframes
│   ├── calibrations/             # LR calibration summaries
│   │   └── dataframes_with_evidence/  # Dataframes with LR columns
│   └── diagnostics/              # Unmapped/unmerged variants
│
├── notebooks/                    # Jupyter notebooks
├── docs/                         # Documentation & examples
└── README.md
```

## Quick Start

### 1. Build Master Dataframes

Run individual gene builders from the project root:

```bash
python src/builders/build_brca1_dataframe.py
python src/builders/build_brca2_dataframe.py
# ... etc for other genes
```

### 2. Calibrate Evidence Strength

Run LR calibration on all genes:

```bash
python src/calibration/run_all_genes.py
```

Or calibrate a single gene:

```bash
python src/calibration/evidence_from_clinvar_with_changes.py \
    --master-path output/master_dataframes/BRCA1_master_dataframe.csv.gz \
    --summary-output output/calibrations/brca1_calibration.csv
```

### 3. Compress Output Files

```bash
python src/scripts/gzip_output_files.py
```

## Calibration Methods

The pipeline supports two calibration approaches:

### Threshold-Based (MaveDB)
Uses curated score thresholds from MaveDB to classify variants into:
- **Abnormal** (loss of function)
- **Normal** (functional)
- **Indeterminate** (uncertain)

### Classification-Based
Uses author-reported functional classifications directly, calculating LR for each individual category.

### Evidence Strength Mapping

Likelihood ratios are mapped to ACMG/AMP evidence strength using Tavtigian thresholds:

| LR Range | Evidence Strength |
|----------|-------------------|
| > 350 | Pathogenic Very Strong |
| > 18.7 | Pathogenic Strong |
| > 4.3 | Pathogenic Moderate |
| > 2.1 | Pathogenic Supporting |
| 0.48 - 2.1 | Indeterminate |
| < 0.48 | Benign Supporting |
| < 0.23 | Benign Moderate |
| < 0.053 | Benign Strong |

## Output Files

### Master Dataframes
`output/master_dataframes/{GENE}_master_dataframe.csv.gz`

Contains:
- Variant identifiers (HGVSc, HGVSp, genomic coordinates)
- Functional assay scores and classifications
- Population frequencies (gnomAD)
- ClinVar annotations
- Computational predictions (REVEL, AlphaMissense, SpliceAI)

### Calibration Summary
`output/calibrations/assay_calibration_summary_with_changes.csv`

Contains per-category LR for each assay:
- Gene, assay name, category
- Calibration method used
- Pathogenic/benign counts
- Likelihood ratio and evidence strength

## Dependencies

- pandas
- numpy
- openpyxl (for Excel files)
- pyarrow (for Parquet files)

## License

[Add license information]

## Citation

[Add citation information]

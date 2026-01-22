"""Run calibration on all gene master dataframes and combine results."""
import pandas as pd
import subprocess
import os
import pathlib

# Setup paths
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent
MASTER_DIR = ROOT_DIR / "output" / "master_dataframes"
CALIBRATION_DIR = ROOT_DIR / "output" / "calibrations"

genes = ['BRCA1', 'BRCA2', 'MSH2', 'PTEN', 'TP53', 'VHL']
all_summaries = []

for gene in genes:
    master_file = MASTER_DIR / f'{gene}_master_dataframe.csv.gz'
    temp_summary = CALIBRATION_DIR / f'temp_calibration_{gene}.csv'
    
    if not master_file.exists():
        print(f'Skipping {gene}: {master_file} not found')
        continue
    
    print(f'Processing {gene}...')
    result = subprocess.run([
        'python', str(SCRIPT_DIR / 'evidence_from_clinvar_with_changes.py'),
        '--master-path', str(master_file),
        '--summary-output', str(temp_summary)
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f'Error for {gene}')
        print(result.stderr[-500:] if result.stderr else "No stderr")
        continue
    
    if temp_summary.exists():
        df = pd.read_csv(temp_summary)
        if len(df) > 0:
            all_summaries.append(df)
            assay_count = df['assay_name'].nunique()
            print(f'  {gene}: {len(df)} rows, {assay_count} assays')
        else:
            print(f'  {gene}: no assays calibrated')
        temp_summary.unlink()

# Combine all summaries
if all_summaries:
    combined = pd.concat(all_summaries, ignore_index=True)
    output_file = CALIBRATION_DIR / 'assay_calibration_summary_with_changes.csv'
    combined.to_csv(output_file, index=False)
    print(f'\n=== COMBINED SUMMARY ===')
    print(f'Output: {output_file}')
    print(f'Total rows: {len(combined)}')
    for gene in combined['Gene'].unique():
        subset = combined[combined['Gene'] == gene]
        assays = subset['assay_name'].nunique()
        print(f'  {gene}: {assays} assays, {len(subset)} categories')
else:
    print('No summaries generated')

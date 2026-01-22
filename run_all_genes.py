"""Run calibration on all gene master dataframes and combine results."""
import pandas as pd
import subprocess
import os

genes = ['BRCA1', 'BRCA2', 'MSH2', 'PTEN', 'TP53', 'VHL']
all_summaries = []

for gene in genes:
    master_file = f'{gene}_master_dataframe.csv.gz'
    temp_summary = f'temp_calibration_{gene}.csv'
    
    if not os.path.exists(master_file):
        print(f'Skipping {gene}: {master_file} not found')
        continue
    
    print(f'Processing {gene}...')
    result = subprocess.run([
        'python', 'evidence_from_clinvar_with_changes.py',
        '--master-path', master_file,
        '--summary-output', temp_summary
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f'Error for {gene}')
        print(result.stderr[-500:])
        continue
    
    if os.path.exists(temp_summary):
        df = pd.read_csv(temp_summary)
        if len(df) > 0:
            all_summaries.append(df)
            assay_count = df['assay_name'].nunique()
            print(f'  {gene}: {len(df)} rows, {assay_count} assays')
        else:
            print(f'  {gene}: no assays calibrated')
        os.remove(temp_summary)

# Combine all summaries
if all_summaries:
    combined = pd.concat(all_summaries, ignore_index=True)
    combined.to_csv('assay_calibration_summary_with_changes.csv', index=False)
    print(f'\n=== COMBINED SUMMARY ===')
    print(f'Total rows: {len(combined)}')
    for gene in combined['Gene'].unique():
        subset = combined[combined['Gene'] == gene]
        assays = subset['assay_name'].nunique()
        print(f'  {gene}: {assays} assays, {len(subset)} categories')
else:
    print('No summaries generated')

#!/usr/bin/env python3
"""
Script to combine VEP output, ClinVar data, and pillar curation sheet.
Fills empty columns in pillar data with information from VEP and ClinVar.
"""

import pandas as pd
import numpy as np
import re
import gzip
import os
from pathlib import Path

print("=== Loading Input Files ===")

# Load pillar data with curation
pillar_file = "./IGVF-cvfg-pillar-project/harmonized_data_v10/pillar_data_with_curation.csv"
pillar_data = pd.read_csv(pillar_file)
print(f"Loaded pillar data: {len(pillar_data)} rows")

# Load MAVE curation sheet for author_provided_Transcript ID
curation_file = "MAVE Curation v3.csv"
curation_data = pd.read_csv(curation_file)
print(f"Loaded curation sheet: {len(curation_data)} rows")

# Load VEP output files - only Ollodart for now
# Jia dataset will use existing final_pillar_data_with_clinvar CSV
vep_files = [
    # "VEP_output/Jia_VEP_output_11192025.txt",  # Commented out - using existing data instead
    "VEP_output/Ollodart_VEP_output_11192025.txt"
]

vep_data_list = []
for vep_file in vep_files:
    if Path(vep_file).exists():
        # VEP files have header starting with #, need to extract it manually
        with open(vep_file, 'r') as f:
            header_line = None
            for line in f:
                if line.startswith('#') and 'Location' in line:
                    # Remove the # and strip whitespace, then split by tab
                    header_line = line.lstrip('#').strip().split('\t')
                    break
        
        if header_line:
            # Read data skipping all comment lines, use extracted header
            vep_df = pd.read_csv(vep_file, sep='\t', comment='#', names=header_line, low_memory=False)
            print(f"Loaded VEP file: {vep_file} ({len(vep_df)} rows) with columns: {len(header_line)}")
        else:
            print(f"Warning: Could not find header in {vep_file}")
            vep_df = pd.DataFrame()
        
        if not vep_df.empty:
            vep_data_list.append(vep_df)

if vep_data_list:
    vep_data = pd.concat(vep_data_list, ignore_index=True)
    print(f"Total VEP rows: {len(vep_data)}")
else:
    print("Warning: No VEP files found!")
    vep_data = pd.DataFrame()

# Load ClinVar data - filter to MSH2 only to speed up processing
clinvar_file = "Clinvar_012025_variant_summary.txt.gz"
print(f"Loading ClinVar data from {clinvar_file}...")
print("Filtering to MSH2 variants only to optimize performance...")
with gzip.open(clinvar_file, 'rt') as f:
    # Read in chunks and filter for MSH2
    chunk_size = 100000
    clinvar_chunks = []
    for chunk in pd.read_csv(f, sep='\t', low_memory=False, chunksize=chunk_size):
        # Filter to MSH2 variants
        msh2_chunk = chunk[chunk['GeneSymbol'] == 'MSH2']
        if len(msh2_chunk) > 0:
            clinvar_chunks.append(msh2_chunk)
        if len(clinvar_chunks) % 10 == 0:
            print(f"  Processed {len(clinvar_chunks) * chunk_size:,} rows, found {sum(len(c) for c in clinvar_chunks):,} MSH2 variants...")
    
    clinvar_data = pd.concat(clinvar_chunks, ignore_index=True) if clinvar_chunks else pd.DataFrame()
print(f"Loaded ClinVar data: {len(clinvar_data)} MSH2 variants (filtered from full dataset)")

print("\n=== Processing VEP Data ===")

# Define SpliceAI DS columns (used for both VEP and existing data)
spliceai_ds_cols = ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 
                    'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']

def normalize_hgvs(hgvs_string):
    """Remove transcript ID prefix from HGVS notation"""
    if pd.isna(hgvs_string) or not isinstance(hgvs_string, str):
        return np.nan
    # Remove transcript ID prefix (e.g., "ENST00000233146.7:c.1_3delinsCAA" -> "c.1_3delinsCAA")
    # Pattern: transcript_id:notation or transcript_id.notation
    match = re.search(r':([cp]\..+)$', hgvs_string)
    if match:
        return match.group(1)
    return hgvs_string

def extract_location_info(location_string):
    """Extract chromosome, start, and end from Location column"""
    if pd.isna(location_string) or not isinstance(location_string, str):
        return np.nan, np.nan, np.nan
    
    # Pattern: "2:47403192-47403194"
    match = re.match(r'(\d+|X|Y|MT):(\d+)-(\d+)', location_string)
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return chrom, start, end
    return np.nan, np.nan, np.nan

def extract_alt_allele(uploaded_allele):
    """Extract alt allele from UPLOADED_ALLELE (part after /)"""
    if pd.isna(uploaded_allele) or not isinstance(uploaded_allele, str):
        return np.nan
    # Pattern: "ATG/CAA" -> "CAA"
    if '/' in uploaded_allele:
        return uploaded_allele.split('/')[1]
    return np.nan

if not vep_data.empty:
    # Process VEP data
    print("Processing VEP data...")
    
    # Extract location information
    location_info = vep_data['Location'].apply(extract_location_info)
    vep_data['Chrom'] = location_info.apply(lambda x: x[0] if isinstance(x, tuple) else np.nan)
    vep_data['hg38_start'] = location_info.apply(lambda x: x[1] if isinstance(x, tuple) else np.nan)
    vep_data['hg38_end'] = location_info.apply(lambda x: x[2] if isinstance(x, tuple) else np.nan)
    
    # Extract ref and alt alleles
    vep_data['ref_allele'] = vep_data['REF_ALLELE']
    vep_data['alt_allele'] = vep_data['UPLOADED_ALLELE'].apply(extract_alt_allele)
    
    # Normalize HGVS notation
    vep_data['hgvs_c_normalized'] = vep_data['HGVSc'].apply(normalize_hgvs)
    vep_data['hgvs_p_normalized'] = vep_data['HGVSp'].apply(normalize_hgvs)
    
    # Keep only relevant columns (spliceai_ds_cols already defined earlier)
    vep_columns = ['Chrom', 'hg38_start', 'hg38_end', 'ref_allele', 'alt_allele',
                   'hgvs_c_normalized', 'hgvs_p_normalized'] + spliceai_ds_cols
    vep_data_processed = vep_data[vep_columns].copy()
    
    print(f"Processed VEP data: {len(vep_data_processed)} rows")
    print(f"Unique HGVSp in VEP: {vep_data_processed['hgvs_p_normalized'].notna().sum()}")
else:
    vep_data_processed = pd.DataFrame()

print("\n=== Processing ClinVar Data ===")

def map_review_status_to_stars(review_status):
    """Map ClinVar ReviewStatus to star rating"""
    if pd.isna(review_status) or not isinstance(review_status, str):
        return 0
    
    review_status = review_status.strip()
    
    if review_status == 'practice guideline':
        return 4
    elif review_status == 'reviewed by expert panel':
        return 3
    elif review_status == 'criteria provided, multiple submitters, no conflicts':
        return 2
    elif review_status in ['criteria provided, single submitter', 
                          'criteria provided, conflicting classifications']:
        return 1
    else:
        return 0

# Process ClinVar data
clinvar_data['clinvar_sig'] = clinvar_data['ClinicalSignificance']
clinvar_data['clinvar_star'] = clinvar_data['ReviewStatus'].apply(map_review_status_to_stars)
clinvar_data['clinvar_date_last_reviewed'] = 'Jan 2025'

# Standardize chromosome format (convert to string, handle X, Y, MT)
clinvar_data['Chrom'] = clinvar_data['Chromosome'].astype(str)
clinvar_data['hg38_start'] = pd.to_numeric(clinvar_data['Start'], errors='coerce')
clinvar_data['hg38_end'] = pd.to_numeric(clinvar_data['Stop'], errors='coerce')
clinvar_data['ref_allele'] = clinvar_data['ReferenceAllele']
clinvar_data['alt_allele'] = clinvar_data['AlternateAllele']

# Keep only relevant columns for merging
clinvar_columns = ['Chrom', 'hg38_start', 'hg38_end', 'ref_allele', 'alt_allele',
                   'clinvar_sig', 'clinvar_star', 'clinvar_date_last_reviewed']
clinvar_data_processed = clinvar_data[clinvar_columns].copy()

print(f"Processed ClinVar data: {len(clinvar_data_processed)} rows")

print("\n=== Processing MAVE Curation Sheet ===")

# Extract author_provided_Transcript ID
curation_transcript = curation_data[['Dataset_tag', 'author_provided_Transcript ID']].copy()
curation_transcript.columns = ['Dataset', 'auth_transcript_id']
print(f"Extracted transcript IDs for {len(curation_transcript)} datasets")

print("\n=== Loading Existing Jia Data ===")

# Load existing final_pillar_data_with_clinvar for Jia dataset
existing_pillar_file = "final_pillar_data_with_clinvar_18_25_gnomad_wREVEL_wAM_wspliceAI_wMutpred2_wtrainvar_gold_standards_expanded_092925.csv"
if Path(existing_pillar_file).exists():
    existing_pillar_data = pd.read_csv(existing_pillar_file, low_memory=False)
    print(f"Loaded existing pillar data: {len(existing_pillar_data)} rows")
    
    # Filter for Jia dataset
    jia_existing = existing_pillar_data[existing_pillar_data['Dataset'] == 'MSH2_Jia_2021'].copy()
    print(f"Found {len(jia_existing)} rows for MSH2_Jia_2021 in existing data")
else:
    print(f"Warning: {existing_pillar_file} not found - skipping Jia data merge")
    jia_existing = pd.DataFrame()

print("\n=== Loading Existing Jia Data ===")

# Load existing final_pillar_data_with_clinvar for Jia dataset
existing_pillar_file = "final_pillar_data_with_clinvar_18_25_gnomad_wREVEL_wAM_wspliceAI_wMutpred2_wtrainvar_gold_standards_expanded_092925.csv"
if Path(existing_pillar_file).exists():
    existing_pillar_data = pd.read_csv(existing_pillar_file, low_memory=False)
    print(f"Loaded existing pillar data: {len(existing_pillar_data)} rows")
    
    # Filter for Jia dataset
    jia_existing = existing_pillar_data[existing_pillar_data['Dataset'] == 'MSH2_Jia_2021'].copy()
    print(f"Found {len(jia_existing)} rows for MSH2_Jia_2021 in existing data")
    
    # Normalize hgvs_p for matching
    jia_existing['hgvs_p_normalized'] = jia_existing['hgvs_p'].apply(normalize_hgvs)
else:
    print(f"Warning: {existing_pillar_file} not found - skipping Jia data merge")
    jia_existing = pd.DataFrame()

print("\n=== Merging Data ===")

# Start with pillar data
result_data = pillar_data.copy()
print(f"Starting with {len(result_data)} rows from pillar data")

# Merge author_provided_Transcript ID from curation sheet
result_data = pd.merge(result_data, curation_transcript, on='Dataset', how='left', suffixes=('', '_curation'))
# Fill auth_transcript_id only if it's blank
result_data['auth_transcript_id'] = result_data['auth_transcript_id'].fillna(result_data.get('auth_transcript_id_curation', pd.Series()))
if 'auth_transcript_id_curation' in result_data.columns:
    result_data = result_data.drop(columns=['auth_transcript_id_curation'])

# Normalize hgvs_p in pillar data for matching
result_data['hgvs_p_normalized'] = result_data['hgvs_p'].apply(normalize_hgvs)

# Merge Jia data from existing CSV (before VEP processing)
if not jia_existing.empty:
    print("\nMerging Jia data from existing CSV...")
    
    # Separate Jia and non-Jia rows
    jia_mask = result_data['Dataset'] == 'MSH2_Jia_2021'
    jia_rows = result_data[jia_mask].copy()
    non_jia_rows = result_data[~jia_mask].copy()
    
    # Merge Jia rows with existing data on normalized hgvs_p
    jia_merged = pd.merge(
        jia_rows,
        jia_existing,
        on='hgvs_p_normalized',
        how='left',
        suffixes=('', '_existing')
    )
    
    # Fill blank columns with existing data (only for Jia rows)
    # Map column names from existing file to pillar data columns
    column_mapping = {
        'Chrom': 'Chrom',
        'hg38_start': 'hg38_start',
        'hg38_end': 'hg38_end',
        'ref_allele': 'ref_allele',
        'alt_allele': 'alt_allele',
        'hgvs_c': 'hgvs_c',
        'auth_transcript_id': 'auth_transcript_id',
        # ClinVar columns (existing file has _2025 suffix)
        'clinvar_sig': 'clinvar_sig_2025',
        'clinvar_star': 'clinvar_star_2025',
        'clinvar_date_last_reviewed': 'clinvar_date_last_reviewed_2025',
        # SpliceAI columns (existing file uses spliceAI_DS_* not SpliceAI_pred_DS_*)
        'SpliceAI_pred_DS_AG': 'spliceAI_DS_AG',
        'SpliceAI_pred_DS_AL': 'spliceAI_DS_AL',
        'SpliceAI_pred_DS_DG': 'spliceAI_DS_DG',
        'SpliceAI_pred_DS_DL': 'spliceAI_DS_DL'
    }
    
    for pillar_col, existing_col in column_mapping.items():
        if pillar_col in jia_merged.columns and existing_col in jia_merged.columns:
            # Only fill if original is blank
            mask = jia_merged[pillar_col].isna() & jia_merged[existing_col].notna()
            jia_merged.loc[mask, pillar_col] = jia_merged.loc[mask, existing_col]
        
        # Also handle direct matches if column names match exactly
        if existing_col in jia_merged.columns and existing_col not in column_mapping.values():
            if existing_col in jia_rows.columns:
                existing_col_merged = f"{existing_col}_existing"
                if existing_col_merged in jia_merged.columns:
                    mask = jia_merged[existing_col].isna() & jia_merged[existing_col_merged].notna()
                    jia_merged.loc[mask, existing_col] = jia_merged.loc[mask, existing_col_merged]
    
    # Drop temporary columns
    existing_cols_to_drop = [c for c in jia_merged.columns if c.endswith('_existing')]
    jia_merged = jia_merged.drop(columns=[c for c in existing_cols_to_drop if c in jia_merged.columns])
    
    # Combine back
    result_data = pd.concat([jia_merged, non_jia_rows], ignore_index=True)
    print(f"After Jia merge: {len(result_data)} rows")
    print(f"Jia rows with data filled: {jia_merged[jia_merged['hgvs_c'].notna()].shape[0]} rows")

# Only process Ollodart with VEP (Jia already processed above)
if not vep_data_processed.empty:
    print("\nMerging VEP data for Ollodart dataset...")
    
    # Filter to only Ollodart rows for VEP processing
    ollodart_mask = result_data['Dataset'] == 'MSH2_Ollodart_2021human_amino_acid_numbering'
    ollodart_rows = result_data[ollodart_mask].copy()
    non_ollodart_rows = result_data[~ollodart_mask].copy()
    
    # Process only Ollodart rows with VEP
    result_data_to_process = ollodart_rows.copy()
    
    # Group VEP data by normalized HGVSp to handle multiple HGVSc
    vep_grouped = vep_data_processed.groupby('hgvs_p_normalized')
    
    # Create expanded rows for multiple HGVSc per HGVSp
    vep_expanded = []
    for hgvs_p_norm, group in vep_grouped:
        if pd.notna(hgvs_p_norm):
            # For each unique HGVSc, create a row
            for hgvs_c_norm in group['hgvs_c_normalized'].dropna().unique():
                # Take first row for genomic coordinates and other fields
                row = group.iloc[0].copy()
                row['hgvs_c_normalized'] = hgvs_c_norm
                vep_expanded.append(row)
    
    if vep_expanded:
        vep_expanded_df = pd.DataFrame(vep_expanded)
        print(f"Expanded VEP data to {len(vep_expanded_df)} rows (handling multiple HGVSc per HGVSp)")
        
        # Handle multiple HGVSc by creating separate rows for each match
        # Process each original pillar row and expand for multiple HGVSc
        result_rows = []
        matches_found = 0
        
        # Process each original pillar row (only Ollodart rows)
        for idx, pillar_row in result_data_to_process.iterrows():
            pillar_hgvs_p_norm = pillar_row['hgvs_p_normalized']
            
            # Find matching VEP rows (only where HGVSp matches hgvs_p)
            if pd.notna(pillar_hgvs_p_norm):
                matching_vep = vep_expanded_df[vep_expanded_df['hgvs_p_normalized'] == pillar_hgvs_p_norm]
            else:
                matching_vep = pd.DataFrame()
            
            if len(matching_vep) > 0:
                matches_found += 1
                # Create a row for each unique HGVSc that matches
                # Group by HGVSc to avoid duplicates
                for hgvs_c_norm in matching_vep['hgvs_c_normalized'].dropna().unique():
                    # Get first row with this HGVSc for genomic coordinates
                    vep_row = matching_vep[matching_vep['hgvs_c_normalized'] == hgvs_c_norm].iloc[0]
                    
                    new_row = pillar_row.copy()
                    
                    # Fill blank columns with VEP data
                    vep_fill_columns = ['Chrom', 'hg38_start', 'hg38_end', 'ref_allele', 'alt_allele']
                    for col in vep_fill_columns:
                        if pd.isna(new_row[col]) and pd.notna(vep_row[col]):
                            new_row[col] = vep_row[col]
                    
                    # Add hgvs_c (only if blank)
                    if pd.isna(new_row['hgvs_c']) and pd.notna(hgvs_c_norm):
                        new_row['hgvs_c'] = hgvs_c_norm
                    
                    # Add SpliceAI columns (take first non-null value if multiple)
                    for col in spliceai_ds_cols:
                        if col in vep_row.index and pd.notna(vep_row[col]):
                            if col not in new_row.index or pd.isna(new_row[col]):
                                new_row[col] = vep_row[col]
                    
                    result_rows.append(new_row)
            else:
                # No match - keep original row as is
                result_rows.append(pillar_row)
        
        ollodart_processed = pd.DataFrame(result_rows).reset_index(drop=True)
        
        # Count unique hgvs_p variants that were matched
        unique_hgvs_p_in_pillar = ollodart_processed['hgvs_p'].dropna().nunique()
        matched_hgvs_p = ollodart_processed[ollodart_processed['hgvs_c'].notna()]['hgvs_p'].dropna().nunique()
        
        print(f"Found {matches_found} Ollodart rows with matching VEP data")
        print(f"Unique hgvs_p variants in Ollodart data: {unique_hgvs_p_in_pillar}")
        print(f"Unique hgvs_p variants matched in VEP output: {matched_hgvs_p} ({matched_hgvs_p/unique_hgvs_p_in_pillar*100:.1f}%)")
        print(f"After VEP merge (with expansion for multiple HGVSc): {len(ollodart_processed)} Ollodart rows")
        
        # Clean up temporary columns
        if 'hgvs_p_normalized' in ollodart_processed.columns:
            ollodart_processed = ollodart_processed.drop(columns=['hgvs_p_normalized'])
        
        # Combine Ollodart (processed) with non-Ollodart rows (Jia already processed)
        result_data = pd.concat([ollodart_processed, non_ollodart_rows], ignore_index=True)
        print(f"Total rows after combining all datasets: {len(result_data)} rows")
        
        print(f"After VEP merge: {len(result_data)} rows")
    else:
        print("No VEP data to merge after expansion")
else:
    print("Skipping VEP merge - no VEP data available")
    # Still need to clean up normalized column if it exists
    if 'hgvs_p_normalized' in result_data.columns:
        result_data = result_data.drop(columns=['hgvs_p_normalized'])

# Merge ClinVar data
print("\nMerging ClinVar data...")
result_data = pd.merge(
    result_data,
    clinvar_data_processed,
    on=['Chrom', 'hg38_start', 'hg38_end', 'ref_allele', 'alt_allele'],
    how='left',
    suffixes=('', '_clinvar')
)

# Fill blank ClinVar columns
clinvar_fill_columns = ['clinvar_sig', 'clinvar_star', 'clinvar_date_last_reviewed']
for col in clinvar_fill_columns:
    clinvar_col = f"{col}_clinvar" if col in result_data.columns else col
    if clinvar_col in result_data.columns:
        # Only fill if original is blank
        mask = result_data[col].isna() & result_data[clinvar_col].notna()
        result_data.loc[mask, col] = result_data.loc[mask, clinvar_col]

# Clean up temporary columns
clinvar_cols_to_drop = [c for c in result_data.columns if c.endswith('_clinvar')]
result_data = result_data.drop(columns=[c for c in clinvar_cols_to_drop if c in result_data.columns])

print(f"After ClinVar merge: {len(result_data)} rows")

# Remove the temporary normalized column
if 'hgvs_p_normalized' in result_data.columns:
    result_data = result_data.drop(columns=['hgvs_p_normalized'])

print("\n=== Summary ===")
print(f"Final dataset: {len(result_data)} rows, {len(result_data.columns)} columns")
print(f"\nColumns filled from VEP:")
vep_filled = ['Chrom', 'hg38_start', 'hg38_end', 'ref_allele', 'alt_allele', 'hgvs_c']
for col in vep_filled:
    if col in result_data.columns:
        filled = result_data[col].notna().sum()
        print(f"  {col}: {filled} non-null values")

print(f"\nColumns filled from ClinVar:")
for col in clinvar_fill_columns:
    if col in result_data.columns:
        filled = result_data[col].notna().sum()
        print(f"  {col}: {filled} non-null values")

print(f"\nSpliceAI columns added:")
for col in spliceai_ds_cols:
    if col in result_data.columns:
        filled = result_data[col].notna().sum()
        print(f"  {col}: {filled} non-null values")

# Save output
output_folder = "../Final MSH2 assays"
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, "pillar_data_with_vep_clinvar.csv")
result_data.to_csv(output_file, index=False)
print(f"\n=== Saved output to: {output_file} ===")


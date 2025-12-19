import pandas as pd
import numpy as np


def add_mave_intervals(master_df: pd.DataFrame, mave_file, gene_symbol: str) -> pd.DataFrame:
    """
    For a given gene, read the MAVE curation file and attach dataset-specific
    interval annotations to the provided master_df.

    For each MAVE row for the gene:
      - For Interval 1–3, if *all three* fields (name, range, class) are null,
        that interval is ignored for that dataset.
      - Otherwise, create constant columns on master_df:
            Interval {i} name {Dataset Name}
            Interval {i} range {Dataset Name}
            Interval {i} MaveDB class {Dataset Name}
        where {Dataset Name} comes from the 'Dataset Name' column (or
        'Dataset_tag' as fallback).

    All legacy generic interval columns like:
        Interval 1 name, Interval 1 range, Interval 1 MaveDB class, ...
    are dropped from master_df before returning.
    """
    if master_df is None or master_df.empty:
        return master_df

    try:
        mave_path = getattr(mave_file, "resolve", None)
        if callable(mave_path):
            mave_path = mave_file
        else:
            mave_path = mave_file

        mave_df = pd.read_csv(mave_file, encoding="latin1")
    except Exception as e:
        print(f"Warning: add_mave_intervals could not read MAVE file '{mave_file}': {e}")
        # Still drop legacy interval columns if present
        for i in range(1, 7):
            for field in ["name", "range", "MaveDB class"]:
                col = f"Interval {i} {field}"
                if col in master_df.columns:
                    master_df = master_df.drop(columns=[col])
        return master_df

    if "Gene (HGNC symbol)" not in mave_df.columns:
        print("Warning: 'Gene (HGNC symbol)' column not found in MAVE file; skipping interval attachment.")
        return master_df

    gene_rows = mave_df[mave_df["Gene (HGNC symbol)"] == gene_symbol]
    interval_cols_added = []
    cols_seen = set()  # Track columns we've already added to prevent duplicates
    
    if gene_rows.empty:
        print(f"Warning: No MAVE rows found for gene {gene_symbol}; skipping interval attachment.")
    else:
        for _, row in gene_rows.iterrows():
            # Explicitly get Dataset Name from MAVE file row
            dataset_name = None
            if "Dataset Name" in row.index:
                dataset_name_val = row["Dataset Name"]
                if pd.notna(dataset_name_val):
                    dataset_name = str(dataset_name_val).strip()
            
            # Fallback to Dataset_tag if Dataset Name is missing
            if not dataset_name and "Dataset_tag" in row.index:
                dataset_tag_val = row["Dataset_tag"]
                if pd.notna(dataset_tag_val):
                    dataset_name = str(dataset_tag_val).strip()
            
            # Final fallback if both are missing
            if not dataset_name:
                dataset_name = "MAVE_dataset"
                print(f"  Warning: No Dataset Name or Dataset_tag found for MAVE row; using '{dataset_name}'")

            for i in range(1, 4):
                name_val = row.get(f"Interval {i} name")
                range_val = row.get(f"Interval {i} range")
                class_val = row.get(f"Interval {i} MaveDB class")

                # Ignore this interval if all three fields are null/NaN
                if all(pd.isna(v) for v in [name_val, range_val, class_val]):
                    continue

                base = f"Interval {i}"
                name_col = f"{base} name {dataset_name}"
                range_col = f"{base} range {dataset_name}"
                class_col = f"{base} MaveDB class {dataset_name}"
                
                # Only add columns if they haven't been added before
                if name_col not in cols_seen:
                    master_df[name_col] = name_val
                    interval_cols_added.append(name_col)
                    cols_seen.add(name_col)
                if range_col not in cols_seen:
                    master_df[range_col] = range_val
                    interval_cols_added.append(range_col)
                    cols_seen.add(range_col)
                if class_col not in cols_seen:
                    master_df[class_col] = class_val
                    interval_cols_added.append(class_col)
                    cols_seen.add(class_col)

    # Drop legacy generic interval columns (1–6) if present
    for i in range(1, 7):
        for field in ["name", "range", "MaveDB class"]:
            col = f"Interval {i} {field}"
            if col in master_df.columns:
                master_df = master_df.drop(columns=[col])
    
    # Reorder columns to place interval columns after clinvar_date_last_reviewed_2025
    if interval_cols_added:
        cols = list(master_df.columns)
        
        # Remove interval columns from their current positions
        for col in interval_cols_added:
            if col in cols:
                cols.remove(col)
        
        # Find the position of clinvar_date_last_reviewed_2025
        if "clinvar_date_last_reviewed_2025" in cols:
            insert_idx = cols.index("clinvar_date_last_reviewed_2025") + 1
            # Insert interval columns after clinvar_date_last_reviewed_2025
            cols[insert_idx:insert_idx] = interval_cols_added
        else:
            # If clinvar_date_last_reviewed_2025 doesn't exist, append at end
            cols.extend(interval_cols_added)
        
        master_df = master_df[cols]

    return master_df



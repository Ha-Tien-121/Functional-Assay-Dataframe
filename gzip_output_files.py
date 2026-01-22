import gzip
import shutil
import os

files_to_gzip = [
    # "VHL_master_dataframe.csv",
    # "VHL_master_dataframe.parquet", 
    # "VHL_variants_not_merged.csv", 
    # "VHL_variants_not_mapped.csv",
    "BRCA1_master_dataframe.csv",
    "BRCA1_master_dataframe.parquet", 
    "BRCA1_variants_not_merged.csv", 
    "BRCA1_variants_not_mapped.csv",
    # "BRCA2_master_dataframe.csv",
    # "BRCA2_master_dataframe.parquet", 
    # "BRCA2_variants_not_merged.csv", 
    # "BRCA2_variants_not_mapped.csv",
    # "PTEN_master_dataframe.csv",
    # "PTEN_master_dataframe.parquet", 
    # "PTEN_variants_not_merged.csv", 
    # "PTEN_variants_not_mapped.csv",
    # "TP53_master_dataframe.csv",
    # "TP53_master_dataframe.parquet", 
    # "TP53_variants_not_merged.csv", 
    # "TP53_variants_not_mapped.csv",
    # "MSH2_master_dataframe.csv",
    # "MSH2_master_dataframe.parquet", 
    # "MSH2_variants_not_merged.csv", 
    # "MSH2_variants_not_mapped.csv",
]

for filename in files_to_gzip:
    if os.path.exists(filename):
        try:
            print(f"Compressing {filename}...")
            with open(filename, 'rb') as f_in:
                output_filename = f"{filename}.gz"
                print(f"Compressing to {output_filename}...")
                with gzip.open(output_filename, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Created {output_filename}")
            # Removing original file to mimic gzip behavior
            os.remove(filename)
            print(f"Removed original {filename}")
            
            # Try to rename new gz to original gz name if possible, otherwise keep .gz
            target_gz = f"{filename}.gz"
            
            # If output_filename is already target_gz, we are done
            if os.path.abspath(output_filename) == os.path.abspath(target_gz):
                print(f"File is already named {target_gz}")
                continue

            try:
                if os.path.exists(target_gz):
                    os.remove(target_gz)
                os.rename(output_filename, target_gz)
                print(f"Renamed to {target_gz}")
            except Exception as e:
                print(f"Could not rename to {target_gz} (locked?): {e}")
                print(f"Kept as {output_filename}")
        except PermissionError:
            print(f"Permission denied: Could not access {filename}. Skipping...")
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    else:
        print(f"File not found: {filename}")


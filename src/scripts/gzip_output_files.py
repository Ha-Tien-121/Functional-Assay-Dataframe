"""Compress output files with gzip."""
import gzip
import shutil
import os
import pathlib

# Setup paths
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent
OUTPUT_DIR = ROOT_DIR / "output"

# Files organized by output subdirectory
FILES_TO_GZIP = {
    "master_dataframes": [
        "BRCA1_master_dataframe.csv",
        "BRCA1_master_dataframe.parquet",
        "BRCA2_master_dataframe.csv",
        "BRCA2_master_dataframe.parquet",
        "MSH2_master_dataframe.csv",
        "MSH2_master_dataframe.parquet",
        "PTEN_master_dataframe.csv",
        "PTEN_master_dataframe.parquet",
        "TP53_master_dataframe.csv",
        "TP53_master_dataframe.parquet",
        "VHL_master_dataframe.csv",
        "VHL_master_dataframe.parquet",
    ],
    "diagnostics": [
        "BRCA1_variants_not_merged.csv",
        "BRCA1_variants_not_mapped.csv",
        "BRCA2_variants_not_merged.csv",
        "BRCA2_variants_not_mapped.csv",
        "MSH2_variants_not_merged.csv",
        "MSH2_variants_not_mapped.csv",
        "PTEN_variants_not_merged.csv",
        "PTEN_variants_not_mapped.csv",
        "TP53_variants_not_merged.csv",
        "TP53_variants_not_mapped.csv",
        "VHL_variants_not_merged.csv",
        "VHL_variants_not_mapped.csv",
    ],
}

def compress_file(filepath: pathlib.Path):
    """Compress a single file with gzip."""
    if not filepath.exists():
        print(f"File not found: {filepath}")
        return False
    
    try:
        print(f"Compressing {filepath.name}...")
        output_path = filepath.with_suffix(filepath.suffix + '.gz')
        
        with open(filepath, 'rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        print(f"  Created {output_path.name}")
        
        # Remove original file
        filepath.unlink()
        print(f"  Removed original {filepath.name}")
        return True
        
    except PermissionError:
        print(f"  Permission denied: {filepath}")
        return False
    except Exception as e:
        print(f"  Error: {e}")
        return False

def main():
    """Compress all output files."""
    total_compressed = 0
    
    for subdir, filenames in FILES_TO_GZIP.items():
        dir_path = OUTPUT_DIR / subdir
        print(f"\n=== {subdir} ===")
        
        for filename in filenames:
            filepath = dir_path / filename
            if compress_file(filepath):
                total_compressed += 1
    
    print(f"\n=== DONE ===")
    print(f"Compressed {total_compressed} files")

if __name__ == "__main__":
    main()

import pandas as pd
import pathlib

def inspect_file(filepath):
    print(f"\n--- Inspecting {filepath} ---")
    path = pathlib.Path(filepath)
    try:
        if path.suffix == '.xlsx':
            df = pd.read_excel(path, nrows=5)
        elif path.suffix == '.tsv':
            df = pd.read_csv(path, sep='\t', nrows=5)
        elif path.suffix == '.csv':
            try:
                df = pd.read_csv(path, nrows=5)
            except Exception:
                # Try with different encoding if default fails
                df = pd.read_csv(path, nrows=5, encoding='latin1')
        else:
            print("Unknown extension")
            return

        print("Columns:", list(df.columns))
        # Print first row as dictionary, but handle encoding issues
        try:
            row_dict = df.iloc[0].to_dict() if not df.empty else "Empty DataFrame"
            print(f"First row: {str(row_dict)}")
        except Exception as e:
            print(f"Could not print first row due to encoding/display: {e}")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")

files = [
    "Example dataframe for MSH2.xlsx",
    "cravat_msh2_cleaned.tsv",
    "Jia MSH2 2021.xlsx",
    "Ollodart2020.xlsx",
    "MSH2_Bouvet.csv",
    "MAVE Curation v3.csv"
]

for f in files:
    inspect_file(f)


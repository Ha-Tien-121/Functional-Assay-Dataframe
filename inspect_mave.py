import pandas as pd

try:
    df = pd.read_csv("MAVE Curation v3.csv", encoding='latin1') # 'charmap' failed, trying latin1
    msh2 = df[df['Gene (HGNC symbol)'] == 'MSH2']
    print(f"Found {len(msh2)} rows for MSH2")
    cols_of_interest = [c for c in df.columns if 'Interval' in c]
    print("Interval columns:", cols_of_interest)
    
    for i, row in msh2.iterrows():
        print(f"\nRow {i}: {row['Dataset_tag']}")
        for c in cols_of_interest:
            val = row[c]
            if pd.notna(val):
                print(f"  {c}: {val}")
except Exception as e:
    print(e)


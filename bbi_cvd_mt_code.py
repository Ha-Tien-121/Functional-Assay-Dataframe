
import pandas as pd
import numpy as np
import re

def data_harmonization(dc):

    aa_dict = {'A': 'Ala','R': 'Arg','N': 'Asn','D': 'Asp','C': 'Cys','E': 'Glu','Q': 'Gln','G': 'Gly','H': 'His',
           'I': 'Ile','L': 'Leu','K': 'Lys','M': 'Met','F': 'Phe','P': 'Pro','S': 'Ser','T': 'Thr','W': 'Trp',
           'Y': 'Tyr','V': 'Val', '=': '=','*': '*', "~": 'del', '-':'del','STOP': '*','Del':'del', "X": 'del'}

    reversed_aa_dict = {v: k for k, v in aa_dict.items()}

    dataset = dc['Dataset']

    x = dc['x']

    y = dc['y']
    
    # Handle both Excel and CSV files
    if x.endswith('.csv'):
        # For CSV files, skiprows may be specified in y
        if isinstance(y, int) and y > 0:
            df1 = pd.read_csv(x, skiprows=y, header=0)
        else:
            df1 = pd.read_csv(x, header=y)
    else:
        df1 = pd.read_excel(x, header = y)

    #create an empty dataframe with column names as follows

    pillar_data = pd.DataFrame(columns=["Dataset", "Gene", "HGNC_id","Chrom","hg19_pos","hg38_start","hg38_end","ref_allele",
                                    "alt_allele","auth_transcript_id","transcript_pos","transcript_ref","transcript_alt",
                                   "aa_pos","aa_ref","aa_alt","hgvs_c","hgvs_p","consequence","auth_reported_score",
                                   "auth_reported_rep_score","auth_reported_func_class","splice_measure","gnomad_MAF",
                                    "clinvar_sig","clinvar_star","clinvar_date_last_reviewed","nucleotide_or_aa",
                                    "rfs_1","rfs_2","rfs_3"])

    # Use df1 only to determine the number of rows to create for pillar_data
    num_rows = len(df1)

    pillar_data = pd.DataFrame(index=range(num_rows), columns=pillar_data.columns)

    for key, value in dc.items():
        if key in pillar_data.columns:
            if pd.api.types.is_scalar(value):
                pillar_data[key] = value
            else:
                # Handle different data structures properly
                if isinstance(value, np.ndarray):
                    pillar_data[key] = value
                elif isinstance(value, pd.Series):
                    pillar_data[key] = value
                elif isinstance(value, pd.DataFrame):
                    if value.shape[1] == 1:
                        pillar_data[key] = value.iloc[:, 0]
                    else:
                        raise ValueError(f"Expected 1D data for {key}, but got shape {value.shape}")
                else:
                    pillar_data[key] = np.asarray(value)


    #convert syntax like M1L or p.M1L to hgvs format with three letter amino acid representations

    def convert_hgvs_p(hgvs_p):
        if pd.isna(hgvs_p):
            return np.nan
        try:
            # Strip 'p.' if present
            hgvs_p_clean = hgvs_p[2:] if hgvs_p.startswith('p.') else hgvs_p

            # Regex pattern: ref_aa (1 letter), position (digits), alt_aa (1 letter or '=')
            match = re.match(r'^([A-Z]|\*)(\d+)([A-Z*=\-~X]|Ter|STOP)$', hgvs_p_clean)
            if not match:
                # If does not match pattern, return original
                return hgvs_p

            ref_a, pos, alt_a = match.groups()

            # For synonymous, alt = ref (or alt = '=')
            if alt_a == '=':
                alt_a = ref_a

            return f"p.{aa_dict[ref_a]}{pos}{aa_dict[alt_a]}"
        except KeyError as e:
            print(f"KeyError converting {dataset}-{hgvs_p}: {e}")
            return hgvs_p

    # Apply
    if dc['hgvs_p_conversion'] == 'Yes':
        pillar_data['hgvs_p'] = pillar_data['hgvs_p'].apply(convert_hgvs_p)


    def construct_hgvs_p_from_aa(row):
        try:
            # Handle missing ref/alt cleanly
            if pd.isna(row['aa_ref']) and pd.isna(row['aa_alt']):
                return np.nan

            ref_aa = row['aa_ref']
            alt_aa = row['aa_alt']
            pos = row['aa_pos']

            # Handle missing aa_pos
            if pd.isna(pos):
                return np.nan

            # Construct using your aa_dict mapping
            return f"p.{aa_dict[ref_aa]}{pos}{aa_dict[alt_aa]}"

        except KeyError as e:
            print(f"KeyError for row {dataset}-{row.name} with aa_ref={row.get('aa_ref')}, aa_alt={row.get('aa_alt')}, aa_pos={row.get('aa_pos')}: {e}")
            # Fallback: return current hgvs_p if it exists
            return row.get('hgvs_p', np.nan)

    # Apply
    if dc['hgvs_from_aa'] == 'Yes':
        pillar_data['hgvs_p'] = pillar_data.apply(construct_hgvs_p_from_aa, axis=1)



    def extract_aa_from_hgvs(row):
        hgvs_p_pattern = re.compile(r'^p\.([A-Z][a-z]{2}|\*|del|Del)(\d+)([A-Z][a-z]{2}|=|\*|Ter|del|Del)$')

        hgvs_p = row['hgvs_p']
        if pd.isna(hgvs_p):
            return pd.Series({'aa_ref': np.nan, 'aa_pos': np.nan, 'aa_alt': np.nan})
        try:
            match = hgvs_p_pattern.match(hgvs_p)
            if not match:
                # If pattern doesn't match, return NaNs for safety
                print(f"No match for row {dataset}-{row.name} with value {hgvs_p}")
                return pd.Series({'aa_ref': np.nan, 'aa_pos': np.nan, 'aa_alt': np.nan})

            ref_aa_3, pos, alt_aa_3 = match.groups()

            # Map ref_aa from 3-letter to 1-letter
            ref_aa = reversed_aa_dict.get(ref_aa_3, np.nan)

            # Handle special cases
            if alt_aa_3 in ['=', 'Ter', '*']:
                if alt_aa_3 == '=':
                    alt_aa = ref_aa
                else:
                    alt_aa = '*'
            else:
                alt_aa = reversed_aa_dict.get(alt_aa_3, np.nan)

            return pd.Series({'aa_ref': ref_aa, 'aa_pos': pos, 'aa_alt': alt_aa})

        except Exception as e:
            print(f"Error for {dataset}-row {row.name} with value {hgvs_p}: {e}")
            return pd.Series({'aa_ref': np.nan, 'aa_pos': np.nan, 'aa_alt': np.nan})

    # Apply
    if dc['aa_from_hgvs'] == 'Yes':
        pillar_data[['aa_ref', 'aa_pos', 'aa_alt']] = pillar_data.apply(extract_aa_from_hgvs, axis=1)

    #if hgvs_c is provided, it uses that information to populate transcript information

    def parse_hgvs_c(row):
        hgvs_c = row['hgvs_c']

        output = {
            'transcript_pos': np.nan,
            'transcript_ref': np.nan,
            'transcript_alt': np.nan
        }

        if pd.isna(hgvs_c):
            return pd.Series(output)

        try:
            # Substitution: c.123A>G, c.*37T>G, c.123+1G>A
            sub_match = re.match(r'c\.([\*\-\d_+\d]+)([ACGT])>([ACGT])$', hgvs_c)
            if sub_match:
                output['transcript_pos'] = sub_match.group(1)
                output['transcript_ref'] = sub_match.group(2)
                output['transcript_alt'] = sub_match.group(3)
                return pd.Series(output)

            # Delins: c.123_125delinsG
            delins_match = re.match(r'c\.([\*\d_\+\-]+)delins([ACGT]+)$', hgvs_c)
            if delins_match:
                output['transcript_pos'] = delins_match.group(1)
                output['transcript_ref'] = 'del'
                output['transcript_alt'] = delins_match.group(2)
                return pd.Series(output)

            # Deletion: c.123_125del or c.123del or c.123_125delAG
            del_match = re.match(r'c\.([\*\d_\+\-]+)del([ACGT]*)$', hgvs_c)
            if del_match:
                output['transcript_pos'] = del_match.group(1)
                output['transcript_ref'] = del_match.group(2) if del_match.group(2) else 'del'
                output['transcript_alt'] = ''
                return pd.Series(output)

            # Insertion: c.123_124insG
            ins_match = re.match(r'c\.([\*\d_\+\-]+)ins([ACGT]+)$', hgvs_c)
            if ins_match:
                output['transcript_pos'] = ins_match.group(1)
                output['transcript_ref'] = ''
                output['transcript_alt'] = ins_match.group(2)
                return pd.Series(output)

            # Duplication: c.942dupG or c.995_998dupAGGC
            dup_match = re.match(r'c\.([\*\d_\+\-]+)dup([ACGT]*)$', hgvs_c)
            if dup_match:
                output['transcript_pos'] = dup_match.group(1)
                output['transcript_ref'] = ''
                output['transcript_alt'] = 'dup' + dup_match.group(2)
                return pd.Series(output)

            # If no pattern matched
            print(f"[WARN] Unparsed HGVS_c: {hgvs_c}")
            return pd.Series(output)

        except Exception as e:
            print(f"[ERROR] Parsing {hgvs_c} failed: {e}")
            return pd.Series(output)

    # Apply
    if dc.get('transcript_from_hgvs_c') == 'Yes':
        pillar_data[['transcript_pos', 'transcript_ref', 'transcript_alt']] = (
            pillar_data.apply(parse_hgvs_c, axis=1)
        )

    return pillar_data

def data_harmonization_loop(dc_list):

    combined_dataframes = []


    for dc in dc_list:
            harmonized_df = data_harmonization(dc)


            if harmonized_df is not None and not harmonized_df.empty:
                combined_dataframes.append(harmonized_df)


    combined_df = pd.concat(combined_dataframes, ignore_index=True) if combined_dataframes else pd.DataFrame()

    return combined_df

# ============================================================================
# MSH2 DATASETS (Ollodart and Jia)
# ============================================================================

# Load Ollodart2020 data
df1 = pd.read_excel("Ollodart2020.xlsx", header=0)

# Configure Ollodart2020 dataset
dc_1 = {"x": "Ollodart2020.xlsx",
               "y": 0,
               "Dataset" :'MSH2_Ollodart_2021human_amino_acid_numbering',
               "Gene" : "MSH2",
               "HGNC_id" : 7325,
               "Chrom" : np.nan,
               "hg19_pos" : np.nan,
               "hg38_start" : np.nan,
               "ref_allele": np.nan,
               "alt_allele": np.nan,
               "auth_transcript_id" : np.nan,
               "transcript_pos" : np.nan,
               "transcript_ref" : np.nan,
               "transcript_alt" : np.nan,
               "aa_pos": np.nan,
               "aa_ref": np.nan,
               "aa_alt": np.nan,
               "hgvs_c": np.nan,
               "hgvs_p" : df1['Human Msh2 Genotypeb'],
               "consequence": np.nan,
               "auth_reported_score": df1['Scoreg'],
               "auth_reported_rep_score": np.nan,
               "auth_reported_func_class": np.nan,
               "splice_measure": 'No',
               "gnomad_MAF": np.nan,
               "clinvar_sig": np.nan,
               "clinvar_star": np.nan,
               "clinvar_date_last_reviewed": np.nan,
               "nucleotide_or_aa": 'aa',
               "hgvs_p_conversion": 'Yes',
               "hgvs_from_aa": 'No',
               "aa_from_hgvs": 'Yes',
               "transcript_from_hgvs_c": 'No'}

# Load Jia MSH2 2021 data
df2 = pd.read_excel("Jia MSH2 2021.xlsx", header=0)

# Configure Jia MSH2 2021 dataset
dc_2 = {"x": "Jia MSH2 2021.xlsx",
               "y": 0,
               "Dataset" :'MSH2_Jia_2021',
               "Gene" : "MSH2",
               "HGNC_id" : 7325,
               "Chrom" : np.nan,
               "hg19_pos" : np.nan,
               "hg38_start" : np.nan,
               "ref_allele": np.nan,
               "alt_allele": np.nan,
               "auth_transcript_id" : np.nan,
               "transcript_pos" : np.nan,
               "transcript_ref" : np.nan,
               "transcript_alt" : np.nan,
               "aa_pos": np.nan,
               "aa_ref": np.nan,
               "aa_alt": np.nan,
               "hgvs_c": np.nan,
               "hgvs_p" : df2['Variant'],
               "consequence": np.nan,
               "auth_reported_score": df2['LOF score'],
               "auth_reported_rep_score": np.nan,
               "auth_reported_func_class": np.nan,
               "splice_measure": 'No',
               "gnomad_MAF": np.nan,
               "clinvar_sig": np.nan,
               "clinvar_star": np.nan,
               "clinvar_date_last_reviewed": np.nan,
               "nucleotide_or_aa": 'aa',
               "hgvs_p_conversion": 'Yes',
               "hgvs_from_aa": 'No',
               "aa_from_hgvs": 'Yes',
               "transcript_from_hgvs_c": 'No'}

# Load Bouvet MSH2 2019 data
df_bouvet = pd.read_csv("MSH2_Bouvet.csv", header=0)

# Extract percentage and functional class from MT assay result
def extract_mt_assay_score(mt_result):
    """Extract percentage from MT assay result (e.g., '77.12% Potentially damaging' -> 77.12)"""
    if pd.isna(mt_result) or not isinstance(mt_result, str):
        return np.nan
    # Extract percentage (number before %)
    match = re.search(r'([\d.]+)%', mt_result)
    if match:
        try:
            return float(match.group(1))
        except:
            return np.nan
    return np.nan

def extract_mt_assay_class(mt_result):
    """Extract functional class from MT assay result (e.g., '77.12% Potentially damaging' -> 'Potentially damaging')"""
    if pd.isna(mt_result) or not isinstance(mt_result, str):
        return np.nan
    # Extract text after percentage
    match = re.search(r'[\d.]+%\s*(.+)', mt_result)
    if match:
        return match.group(1).strip()
    return np.nan

bouvet_score = df_bouvet['MT assay result'].apply(extract_mt_assay_score)
bouvet_class = df_bouvet['MT assay result'].apply(extract_mt_assay_class)

# Configure Bouvet MSH2 2019 dataset
dc_bouvet = {"x": "MSH2_Bouvet.csv",
               "y": 0,
               "Dataset" :'MSH2_Bouvet_2019',
               "Gene" : "MSH2",
               "HGNC_id" : 7325,
               "Chrom" : np.nan,
               "hg19_pos" : np.nan,
               "hg38_start" : np.nan,
               "ref_allele": np.nan,
               "alt_allele": np.nan,
               "auth_transcript_id" : np.nan,
               "transcript_pos" : np.nan,
               "transcript_ref" : np.nan,
               "transcript_alt" : np.nan,
               "aa_pos": np.nan,
               "aa_ref": np.nan,
               "aa_alt": np.nan,
               "hgvs_c": np.nan,
               "hgvs_p" : df_bouvet['Protein'],
               "consequence": np.nan,
               "auth_reported_score": bouvet_score,
               "auth_reported_rep_score": np.nan,
               "auth_reported_func_class": bouvet_class,
               "splice_measure": 'No',
               "gnomad_MAF": np.nan,
               "clinvar_sig": np.nan,
               "clinvar_star": np.nan,
               "clinvar_date_last_reviewed": np.nan,
               "nucleotide_or_aa": 'aa',
               "hgvs_p_conversion": 'Yes',  # Convert single-letter AA codes to 3-letter
               "hgvs_from_aa": 'No',
               "aa_from_hgvs": 'Yes',  # Extract AA info from hgvs_p
               "transcript_from_hgvs_c": 'No'}

# Process MSH2 datasets
print("\n=== Processing MSH2 Datasets ===")
msh2_dc_list = [dc_1, dc_2, dc_bouvet]
msh2_combined_df = data_harmonization_loop(msh2_dc_list)

# ============================================================================
# TP53 FUNK DATASET
# ============================================================================

# Load TP53 Funk et al data
df3 = pd.read_csv("TP53 Funk et al.csv", skiprows=2, header=0)

# Filter to only include missense variants (type_p == 'mis')
df3 = df3[df3['type_p'] == 'mis'].copy()
print(f"\n=== Processing TP53 Funk Dataset ===")
print(f"Filtered TP53 Funk data to {len(df3)} missense variants (type_p == 'mis')")

# Extract HGVS notation from columns that include transcript IDs
# hg38_cDNA format: "NM_000546.6:c.376-12_376-10del" -> extract "c.376-12_376-10del"
# hg38_protein format: "NP_000537.3:p.Tyr126Val" -> extract "p.Tyr126Val" (or "NP_000537.3:p.?" -> np.nan)
def extract_hgvs_c(hgvs_string):
    if pd.isna(hgvs_string) or not isinstance(hgvs_string, str):
        return np.nan
    # Pattern: transcript_id:notation or transcript_id.notation
    match = re.search(r':([cp]\..+)$', hgvs_string)
    if match:
        return match.group(1)
    return np.nan

def extract_hgvs_p(hgvs_string):
    if pd.isna(hgvs_string) or not isinstance(hgvs_string, str) or hgvs_string.endswith('p.?'):
        return np.nan
    # Pattern: transcript_id:notation or transcript_id.notation
    match = re.search(r':([cp]\..+)$', hgvs_string)
    if match:
        return match.group(1)
    return np.nan

hgvs_c_extracted = df3['hg38_cDNA'].apply(extract_hgvs_c)
hgvs_p_extracted = df3['hg38_protein'].apply(extract_hgvs_p)

# Configure TP53 Funk et al 2025 dataset
dc_3 = {"x": "TP53 Funk et al.csv",
               "y": 2,  # skiprows=2 for CSV
               "Dataset" :'TP53_Funk_2025',
               "Gene" : "TP53",
               "HGNC_id" : 11998,
               "Chrom" : np.nan,
               "hg19_pos" : np.nan,
               "hg38_start" : np.nan,
               "ref_allele": np.nan,
               "alt_allele": np.nan,
               "auth_transcript_id" : np.nan,
               "transcript_pos" : np.nan,
               "transcript_ref" : np.nan,
               "transcript_alt" : np.nan,
               "aa_pos": np.nan,
               "aa_ref": np.nan,
               "aa_alt": np.nan,
               "hgvs_c": hgvs_c_extracted,
               "hgvs_p" : hgvs_p_extracted,
               "consequence": np.nan,
               "auth_reported_score": df3['rfs_median'],
               "auth_reported_rep_score": np.nan,
               "auth_reported_func_class": np.nan,
               "splice_measure": 'No',
               "gnomad_MAF": np.nan,
               "clinvar_sig": np.nan,
               "clinvar_star": np.nan,
               "clinvar_date_last_reviewed": np.nan,
               "nucleotide_or_aa": 'aa',
               "hgvs_p_conversion": 'No',  # Already in 3-letter format
               "hgvs_from_aa": 'No',
               "aa_from_hgvs": 'Yes',  # Extract AA info from hgvs_p
               "transcript_from_hgvs_c": 'Yes',  # Extract transcript info from hgvs_c
               "rfs_1": df3['rfs_1'],
               "rfs_2": df3['rfs_2'],
               "rfs_3": df3['rfs_3']}

# Process TP53 Funk dataset
tp53_dc_list = [dc_3]
tp53_combined_df = data_harmonization_loop(tp53_dc_list)

pd.set_option('display.max_columns', None)

# ============================================================================
# PROCESS AND SAVE MSH2 DATASETS
# ============================================================================

def process_and_save_dataset(df, dataset_name, output_folder):
    """Process a dataset and save to separate files"""
    import os
    # Drop lines without a functional score
    df = df[
        ((df['Dataset'] != 'F9_Popp_2025_model') &
         (~df['auth_reported_score'].isna())) |
        ((df['Dataset'] == 'F9_Popp_2025_model') &
         (~df['auth_reported_func_class'].isna()))
    ]
    
    # for OTC_Lo drop variants in the SMG loop, referred to authors as not useable
    df = df[df['auth_reported_func_class'] != 'SMG Loop']
    
    # PRKN VEP not working, dropping for now
    df = df[df['Dataset'] != 'PRKN_Clausen_2024']
    
    # Give each variant a unique identifier
    df['ID'] = df.apply(lambda row: f"{row['Dataset']}_var{row.name + 1}", axis=1)
    
    # Save the file
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f'{dataset_name}_pillar_data_combined_df.csv')
    df.to_csv(output_file, index=False)
    print(f"Saved {dataset_name} to {output_file}")
    
    # Load curation sheet and merge
    curation_data = pd.read_csv("MAVE Curation v3.csv", header = 0)
    
    # Merge and add additional information from the big curation sheet
    df_merge = pd.merge(df, curation_data[['Dataset_tag', 'MaveDB Score Set URN', 'Ensembl_transcript_ID',
                                               "Ref_seq_transcript_ID","Model_system","Assay Type",
                                              "Phenotype Measured ontology term","Molecular or Biological Process Investigated (GO term)","IGVF_produced",'Interval 1 name', 'Interval 1 range', 'Interval 1 MaveDB class',
                                             'Interval 2 name', 'Interval 2 range', 'Interval 2 MaveDB class', 'Interval 3 name',
                                             'Interval 3 range', 'Interval 3 MaveDB class', 'Interval 4 name', 'Interval 4 range',
                                             'Interval 4 MaveDB class', 'Interval 5 name', 'Interval 5 range', 'Interval 5 MaveDB class',
                                             'Interval 6 name', 'Interval 6 range', 'Interval 6 MaveDB class']],
                                               left_on='Dataset', right_on = "Dataset_tag", how='left')
    
    # Manually fill transcript IDs for Bouvet (not in curation sheet, use same as Jia)
    bouvet_mask = df_merge['Dataset'] == 'MSH2_Bouvet_2019'
    if bouvet_mask.any():
        df_merge.loc[bouvet_mask, 'Ensembl_transcript_ID'] = 'ENST00000233146.7'
        df_merge.loc[bouvet_mask, 'Ref_seq_transcript_ID'] = 'NM_000251.2'
    
    df_merge = df_merge.drop_duplicates()
    
    curation_output_file = os.path.join(output_folder, f'{dataset_name}_pillar_data_with_curation.csv')
    df_merge.to_csv(curation_output_file, index = False)
    print(f"Saved {dataset_name} with curation to {curation_output_file}")
    
    return df_merge

# Process MSH2 datasets
print("\n=== Processing MSH2 Datasets ===")
msh2_processed = process_and_save_dataset(msh2_combined_df, 'MSH2', './IGVF-cvfg-pillar-project/harmonized_data_v10')

# Skip TP53 Funk processing - not needed for Bouvet addition
# print("\n=== Processing TP53 Funk Dataset ===")
# tp53_processed = process_and_save_dataset(tp53_combined_df, 'TP53_Funk', './IGVF-cvfg-pillar-project/harmonized_data_v10')

#some datasets need to have their ref/alt and transcript ref/alt set to the same allele given they are on
#the positive strand

# condition = (p_data['Dataset'] == 'BRCA2_Hu_2024') | (p_data['Dataset'] == 'BRCA2_Sahu_2023_exon13_SGE') | (p_data['Dataset'] == 'BRCA2_Sahu_2023_exon13_Olaparib') | (p_data['Dataset'] == 'BRCA2_Sahu_2023_exon13_Cisplatin')

# #setting ref/alt and transcript ref/alt to the same allele
# p_data.loc[condition,'ref_allele'] = p_data.loc[condition,'transcript_ref']

# p_data.loc[condition,'alt_allele'] = p_data.loc[condition,'transcript_alt']

# ============================================================================
# GENERATE VEP INPUT FILES
# ============================================================================

import os
import re
import pandas as pd

# Specify SGE datasets
chromosome_datasets = [
    "BARD1_unpublished", "CTCF_unpublished", 'BRCA2_unpublished',
    'PALB2_unpublished','RAD51D_unpublished','SFPQ_unpublished','XRCC2_unpublished',
    'BAP1_Waters_2024','DDX3X_Radford_2023_cLFC_day15',"BRCA2_Sahu_2025_HDR",
    'BRCA2_Huang_2025_HDR',"RAD51C_Olvera-León_2024_z_score_D4_D14",'NPC1_Erwood_2022_RPE1',
    "NPC1_Erwood_2022_HEK293T",'LARGE1_Ma_2024','FKRP_Ma_2024', 'JAG1_Gilbert_2024'
]

#some datasets have a c.provided already

cdot_datasets = ['BRCA1_Findlay_2018',"VHL_Buckley_2024","BRCA2_Sahu_2023_exon13_SGE", 'BRCA2_Hu_2024',
                 "BRCA2_Sahu_2023_exon13_Cisplatin","BRCA2_Sahu_2023_exon13_Olaparib","RHO_Wan_2019",
                "CBS_Sun_2020_high_B6","CBS_Sun_2020_low_B6",'CARD11_Meitlis_2020_DMSO_no_introns',
                'CARD11_Meitlis_2020_Ibrutinib_no_introns']

amino_acid_to_codon = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'],
    'Ter': ['TAA', 'TAG', 'TGA'], '*': ['TAA', 'TAG', 'TGA'],
    'C': ['TGT', 'TGC'], 'W': ['TGG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '-': ['del']

}

output_folder = "./VEP_input/"
os.makedirs(output_folder, exist_ok=True)

# Process VEP input files - only generate for Bouvet (skip Jia and Ollodart)
# Only process Bouvet from MSH2 datasets, skip TP53 Funk
bouvet_data = msh2_processed[msh2_processed['Dataset'] == 'MSH2_Bouvet_2019']

if not bouvet_data.empty:
    # Process datasets
    for dataset, group in bouvet_data.groupby('Dataset'):
        output_file = os.path.join(output_folder, f"{dataset}_vep_input.txt")

        with open(output_file, "w") as file:
            # Check if the dataset requires chromosome-based notation
            if dataset in chromosome_datasets:
                for i, (chrom, pos, ref, alt) in enumerate(zip(group['Chrom'], group['hg38_start'], group['ref_allele'], group['alt_allele'])):
                    try:
                        pos = int(pos)
                        chrom_str = str(chrom).upper()
                        if chrom_str == "X":
                            chrom_fmt = "23"
                        elif chrom_str == "Y":
                            chrom_fmt = "24"
                        elif chrom_str in ["MT", "M"]:
                            chrom_fmt = "MT"
                        else:
                            chrom_int = int(chrom_str)
                            chrom_fmt = f"{chrom_int:02}"  # leading zeros for 1–9

                        chro_s = f"NC_0000{chrom_fmt}:g.{pos}{ref}>{alt}"
                        file.write(chro_s + "\n")
                    except KeyError as e:
                        print(f"KeyError for row {i}{dataset} with value {chrom}: {e}")
                        continue
                    except ValueError as e:
                        print(f"ValueError for row {i} with position {pos}: {e}")
                        continue

            elif dataset in cdot_datasets:
                for i, (hgvs_c,ID) in enumerate(zip(group['hgvs_c'],group['Ensembl_transcript_ID'])):
                    if pd.notna(hgvs_c):
                        transcript = re.sub(r"(ENST\d+)\.\d+", r"\1", ID)
                        file.write(f"{transcript}:{hgvs_c}\n")
            else:
                # Use transcript-based notation for other datasets
                for i, (value, pos, ID, aa_ref) in enumerate(zip(group['aa_alt'], group['aa_pos'], group['Ensembl_transcript_ID'], group['aa_ref'])):
                    try:
                        amino_acid = aa_ref if value in ['0', '='] else value

                        if pd.notna(ID) and isinstance(ID, str) and amino_acid in amino_acid_to_codon:
                            transcript = re.sub(r"(ENST\d+)\.\d+", r"\1", ID)
                            codons = amino_acid_to_codon[amino_acid]
                            nucleotide_pos = (int(float(pos)) * 3) - 2
                            for codon in codons:
                                if amino_acid == '-':
                                    file.write(f"{transcript}:c.{int(nucleotide_pos)}_{int(nucleotide_pos + 2)}{codon}\n")
                                else:
                                    file.write(f"{transcript}:c.{int(nucleotide_pos)}_{int(nucleotide_pos+2)}delins{codon}\n")
                    except KeyError as e:
                        print(f"KeyError for row {i} with value {amino_acid}: {e}")
                        continue
                    except ValueError as e:
                        print(f"ValueError for row {i} with position {pos}: {e}")
                        continue
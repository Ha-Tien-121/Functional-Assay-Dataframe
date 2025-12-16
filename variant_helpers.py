import re
import pandas as pd
import pathlib

# Three letter to one letter map
AA_3_TO_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'Sec': 'U', # Selenocysteine
    'Stop': '*'
}

AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items()}

def parse_hgvsp(hgvsp):
    """
    Parses HGVSp string (e.g., 'p.Met779Arg', 'M779R', 'p.M779R', 'p.Pro283=', 'p.Cys218Ter') 
    into (pos, ref, alt, short_code).
    
    Returns (None, None, None, None) if parsing fails.
    """
    if not isinstance(hgvsp, str) or pd.isna(hgvsp):
        return None, None, None, None

    clean = hgvsp.strip().replace('p.', '')
    
    # Handle synonymous cases explicitly if not caught by regex
    # e.g. p.Pro283=
    
    # Pattern 1: 3-letter (e.g., Met779Arg, Cys218Ter, Pro283=)
    # Matches Ref(3)Pos(N)Alt(3/Special)
    match_3 = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop|\*|=|del|Del)$', clean)
    if match_3:
        ref_3, pos, alt_token = match_3.groups()
        ref = AA_3_TO_1.get(ref_3, '?')
        
        if alt_token in ['Ter', 'Stop', '*']:
            alt = '*'
        elif alt_token in ['del', 'Del']:
            alt = 'del'
        elif alt_token == '=':
            alt = ref # Synonymous
        else:
            alt = AA_3_TO_1.get(alt_token, '?')
        
        return int(pos), ref, alt, f"{ref}{pos}{alt}"

    # Pattern 2: 1-letter (e.g., M779R, P283=, C218*)
    match_1 = re.match(r'^([A-Z])(\d+)([A-Z]|\*|=|del)$', clean)
    if match_1:
        ref, pos, alt_token = match_1.groups()
        if alt_token == '*':
            alt = '*'
        elif alt_token == '=':
            alt = ref          # synonymous
        elif alt_token == 'del':
            alt = 'del'
        else:
            alt = alt_token

        return int(pos), ref, alt, f"{ref}{pos}{alt}"
    
    # Pattern 3: Complex cases (fs, del range, ext, etc.) - retained from original
    # e.g. Gly669LysfsTer?, Ala2_Met26del, Met1?, Ter935ArgextTer7
    match_complex = re.match(r'^([A-Z][a-z]{2})(\d+)(.*)$', clean)
    if match_complex:
        ref_3, pos, remainder = match_complex.groups()
        ref = AA_3_TO_1.get(ref_3, '?')
        
        alt = '?'
        if 'fs' in remainder:
            alt = 'fs'
        elif 'del' in remainder:
            alt = 'del'
        elif 'ins' in remainder:
            alt = 'ins'
        elif 'dup' in remainder:
            alt = 'dup'
        elif 'ext' in remainder:
            alt = 'ext'
        elif remainder == '?':
            alt = '?'
        else:
             # Try to find an alt AA at start of remainder e.g. Lys...
             match_alt = re.match(r'^([A-Z][a-z]{2})', remainder)
             if match_alt:
                 alt = AA_3_TO_1.get(match_alt.group(1), '?')
             elif remainder in ['Ter', '*', 'Stop']:
                 alt = '*'

        return int(pos), ref, alt, f"{ref}{pos}{alt}"

    return None, None, None, None

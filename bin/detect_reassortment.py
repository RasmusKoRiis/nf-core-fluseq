#!/usr/bin/env python3
"""
Summarise BLAST hits for an influenza sample and flag possible reassortment.

• Subject header is expected as:  >STRAIN|LINEAGE|SEGMENT|ACCESSION
  e.g.  A/Victoria/2570/2019|pdm09|PB1|EPI_ISL_528951

Reassortment flag:
  No       – every high‑quality segment points to the *same* STRAIN
  Yes      – different STRAINs among high‑quality segments
  Less Likely - at least one segment < IDENTITY_THRESHOLD
  Unknown  – at least one segment missing 
"""

import argparse
import re
import pandas as pd

EXPECTED_SEGMENTS   = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
IDENTITY_THRESHOLD  = 80.0  # %

SEGMENT_ALIASES = {
    'M': 'MP',
    'M1': 'MP',
    'M2': 'MP',
    'HA1': 'HA',
    'HA2': 'HA',
    'NS1': 'NS',
    'NS2': 'NS',
}

FALLBACK_SEGMENT_PATTERNS = [
    'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS',
    'M2', 'M1', 'M', 'HA2', 'HA1', 'NS2', 'NS1',
]


def canonicalize_segment(token) -> str:
    """Normalize observed segment labels to report schema labels."""
    if token is None or pd.isna(token):
        return pd.NA
    t = str(token).strip().upper()
    if not t:
        return pd.NA
    if t in EXPECTED_SEGMENTS:
        return t
    if t in SEGMENT_ALIASES:
        return SEGMENT_ALIASES[t]
    return pd.NA


def infer_segment_from_qseqid(qseqid) -> str:
    """
    Infer segment from heterogeneous query IDs.
    Supports legacy IDs like sample_HA and more relaxed forms like
    sample_HA_flumut or sample|01-HA-H5N1.
    """
    if qseqid is None or pd.isna(qseqid):
        return pd.NA

    raw = str(qseqid).strip()
    if not raw:
        return pd.NA
    upper = raw.upper()

    # 1) Legacy behavior: last token after underscore.
    legacy = re.search(r'_([^_]+)$', upper)
    if legacy:
        seg = canonicalize_segment(legacy.group(1))
        if not pd.isna(seg):
            return seg

    # 2) Token scan across common delimiters.
    tokens = [tok for tok in re.split(r'[^A-Z0-9]+', upper) if tok]
    for tok in tokens:
        seg = canonicalize_segment(tok)
        if not pd.isna(seg):
            return seg

    # 3) Boundary-based fallback scan in full ID string.
    for pat in FALLBACK_SEGMENT_PATTERNS:
        if re.search(rf'(?<![A-Z0-9]){pat}(?![A-Z0-9])', upper):
            seg = canonicalize_segment(pat)
            if not pd.isna(seg):
                return seg

    return pd.NA

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('--blast',   required=True, help='BLAST outfmt 6 table')
    p.add_argument('--output',  required=True, help='Single‑line CSV result')
    p.add_argument('--sample',  required=True, help='Sample ID')
    return p.parse_args()

def main() -> None:
    args = parse_args()

    # Load BLAST tabular file
    cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
            'qstart','qend','sstart','send','evalue','bitscore']
    df = pd.read_csv(args.blast, sep='\t', names=cols)

    # Infer segment from query id using tolerant parsing.
    df['segment'] = df['qseqid'].apply(infer_segment_from_qseqid)

    # Strain name is EVERYTHING before the first “|” in the *subject* id
    df['strain']  = df['sseqid'].str.split('|').str[0]

    # Best hit per segment = highest % identity
    best = (df.dropna(subset=['segment'])
              .sort_values('pident', ascending=False)
              .drop_duplicates('segment'))

    # Build one‑line output
    row = {'Sample': args.sample}
    different_strains = set()
    missing_segment   = False         
    low_identity_hit  = False        

    for seg in EXPECTED_SEGMENTS:
        hit = best[best['segment'] == seg]
        if hit.empty:
            row[seg] = 'Missing'
            missing_segment = True    
            continue

        hit = hit.iloc[0]
        pid = round(hit['pident'], 1)
        if pid < IDENTITY_THRESHOLD:
            row[seg] = f"TooLow({pid})"
            low_identity_hit = True   
            continue

        strain = hit['strain']
        row[seg] = f"{strain}({pid})"
        different_strains.add(strain)

    # Decide flag
    if missing_segment:                           
        row['Reassortment'] = 'Unknown'           
    elif low_identity_hit:                        
        row['Reassortment'] = 'Unknown (Less Likely)'  
    else:
        row['Reassortment'] = 'No' if len(different_strains) == 1 else 'Yes'

    # Save
    pd.DataFrame([row]).to_csv(args.output, index=False)

if __name__ == '__main__':
    main()

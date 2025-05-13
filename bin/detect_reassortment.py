#!/usr/bin/env python3
"""
Summarise BLAST hits for an influenza sample and flag possible reassortment.

• Subject header is expected as:  >STRAIN|LINEAGE|SEGMENT|ACCESSION
  e.g.  A/Victoria/2570/2019|pdm09|PB1|EPI_ISL_528951

Reassortment flag:
  No       – every high‑quality segment points to the *same* STRAIN
  Yes      – different STRAINs among high‑quality segments
  Unknown  – at least one segment missing or < IDENTITY_THRESHOLD
"""

import argparse
import pandas as pd

EXPECTED_SEGMENTS   = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
IDENTITY_THRESHOLD  = 80.0  # %

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

    # Segment name comes after the last “_” in the *query* id
    df['segment'] = df['qseqid'].str.extract(r'_([^_]+)$')

    # Strain name is EVERYTHING before the first “|” in the *subject* id
    df['strain']  = df['sseqid'].str.split('|').str[0]

    # Best hit per segment = highest % identity
    best = (df.sort_values('pident', ascending=False)
              .drop_duplicates('segment'))

    # Build one‑line output
    row = {'Sample': args.sample}
    different_strains = set()
    quality_ok = True

    for seg in EXPECTED_SEGMENTS:
        hit = best[best['segment'] == seg]
        if hit.empty:
            row[seg] = 'Missing'
            quality_ok = False
            continue

        hit = hit.iloc[0]
        pid = round(hit['pident'], 1)
        if pid < IDENTITY_THRESHOLD:
            row[seg] = f"TooLow({pid})"
            quality_ok = False
            continue

        strain = hit['strain']
        row[seg] = f"{strain}({pid})"
        different_strains.add(strain)

    # Decide flag
    if not quality_ok:
        row['Reassortment'] = 'Unknown'
    else:
        row['Reassortment'] = 'No' if len(different_strains) == 1 else 'Yes'

    # Save
    pd.DataFrame([row]).to_csv(args.output, index=False)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', required=True, help="BLAST tabular result")
    parser.add_argument('--output', required=True, help="Output CSV (one-line)")
    parser.add_argument('--sample', required=True, help="Sample ID")
    return parser.parse_args()

def main():
    args = parse_args()

    EXPECTED_SEGMENTS = ['PB1', 'PB2', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    IDENTITY_THRESHOLD = 80.0  # percent

    # Load BLAST result
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_csv(args.blast, sep='\t', names=cols)

    # Extract segment from query (e.g., 2495513-INFB_PB2 → PB2)
    df['segment'] = df['qseqid'].str.extract(r'_([^_]+)$')

    # Extract lineage from subject (e.g., Ref_H3N2_PB1 → H3N2)
    df['lineage'] = df['sseqid'].str.extract(r'_([^_]+)_')

    # Sort to get best hit per segment
    best_hits = df.sort_values('pident', ascending=False).drop_duplicates('segment')

    # Dictionary for all segment calls
    output = {'Sample': args.sample}
    lineage_list = []
    segments_flagged = False

    for seg in EXPECTED_SEGMENTS:
        segment_row = best_hits[best_hits['segment'] == seg]
        if segment_row.empty:
            output[seg] = 'Missing'
            segments_flagged = True
        else:
            row = segment_row.iloc[0]
            lineage = row['lineage']
            identity = round(row['pident'], 1)

            if identity < IDENTITY_THRESHOLD:
                output[seg] = f"TooLow({identity})"
                segments_flagged = True
            else:
                output[seg] = f"{lineage}({identity})"
                lineage_list.append(lineage)

    # Determine reassortment flag
    if segments_flagged:
        output['Reassortment'] = 'Unknown'
    else:
        output['Reassortment'] = 'No' if len(set(lineage_list)) == 1 else 'Yes'

    pd.DataFrame([output]).to_csv(args.output, index=False)

if __name__ == "__main__":
    main()

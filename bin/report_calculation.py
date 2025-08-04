
#!/usr/bin/env python3
"""
Process influenza NGS QC CSV files.

Adds two columns:

1. **NGS_QC_Sum** – segment‑wise QC issues in the format:
      HA:MS|PB1:LC,FS|NP:LC
2. **GISAID_Comment** – "Review" when any QC issues are present (i.e. NGS_QC_Sum not empty).

Segments considered
-------------------
HA, NA, MP, NP, NS, PA, PB1, PB2

Issue codes
-----------
FS – frameshift detected
LC – low coverage (<80 %)
MS – mixed sites (>3)

Usage
-----
$ python process_qc.py input.csv              # writes input_processed.csv
$ python process_qc.py input.csv -o output.csv
"""

import argparse
from pathlib import Path
import pandas as pd

# -----------------------------------------
# Column mappings
# -----------------------------------------
COVERAGE_COL = {
    'HA': 'Coverage-HA',
    'NA': 'Coverage-NA',
    'MP': 'Coverage-M',
    'NP': 'Coverage-NP',
    'NS': 'Coverage-NS',
    'PA': 'Coverage-PA',
    'PB1': 'Coverage-PB1',
    'PB2': 'Coverage-PB2',
}

FRAMESHIFT_COLS = {
    'HA': ['frameShifts HA1', 'frameShifts HA2'],
    'NA': ['frameShifts NA'],
    'MP': ['frameShifts M1', 'frameShifts M2'],
    'NP': ['frameShifts NP'],
    'NS': ['frameShifts NS'],
    'PA': ['frameShifts PA'],
    'PB1': ['frameShifts PB1'],
    'PB2': ['frameShifts PB2'],
}

MIXED_COLS = {
    'HA': ['Nextclade Mixed Sites HA1', 'Nextclade Mixed Sites HA2'],
    'NA': ['Nextclade Mixed Sites NA'],
    'MP': ['Nextclade Mixed Sites M1', 'Nextclade Mixed Sites M2'],
    'NP': ['Nextclade Mixed Sites NP'],
    'NS': ['Nextclade Mixed Sites NS'],
    'PA': ['Nextclade Mixed Sites PA'],
    'PB1': ['Nextclade Mixed Sites PB1'],
    'PB2': ['Nextclade Mixed Sites PB2'],
}

SEGMENT_ORDER = ['HA', 'NA', 'MP', 'NP', 'NS', 'PA', 'PB1', 'PB2']


# -----------------------------------------
# Core summarisation logic
# -----------------------------------------
def qc_summary(row: pd.Series) -> str:
    """Return QC summary string for one row."""
    segments_out = []
    for seg in SEGMENT_ORDER:
        issues = []

        # Frameshift: any column not equal 'No frameShifts'
        for col in FRAMESHIFT_COLS[seg]:
            val = row.get(col)
            if pd.notna(val) and str(val).strip().lower() != 'no frameshifts':
                issues.append('FS')
                break

        # Low coverage
        cov_val = row.get(COVERAGE_COL[seg])
        if pd.notna(cov_val) and cov_val < 80:
            issues.append('LC')

        # Mixed sites
        ms_sum = sum(
            float(row.get(col, 0) or 0)
            for col in MIXED_COLS[seg]
            if pd.notna(row.get(col))
        )
        if ms_sum > 3:
            issues.append('MS')

        if issues:
            segments_out.append(f"{seg}:{','.join(sorted(issues))}")

    return '|'.join(segments_out)


def process_file(in_csv: Path, out_csv: Path) -> None:
    df = pd.read_csv(in_csv)
    df['NGS_QC_Sum'] = df.apply(qc_summary, axis=1)
    # New column: flag for GISAID comment
    df['GISAID_Comment'] = df['NGS_QC_Sum'].apply(lambda x: 'Review' if str(x).strip() else '')
    df.to_csv(out_csv, index=False)
    print(f"Wrote processed file to {out_csv}")  # noqa: T201


# -----------------------------------------
# CLI
# -----------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add NGS_QC_Sum and GISAID_Comment columns to influenza QC CSV files"
    )
    parser.add_argument("input", type=Path, help="Input CSV file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output CSV file (default: <input>_processed.csv)",
    )
    args = parser.parse_args()

    out_path = args.output or args.input.with_name(
        f"{args.input.stem}_processed.csv"
    )
    process_file(args.input, out_path)


if __name__ == "__main__":
    main()

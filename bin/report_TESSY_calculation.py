#!/usr/bin/env python3
"""
Process influenza NGS QC CSV files.

Adds two columns:

1. NGS_QC_Sum – segment-wise QC issues, e.g. HA:MS|PB1:LC,FS|NP:LC
2. GISAID_Comment – "Review" when any QC issues are present, else "NA".
"""
import argparse
from pathlib import Path
import pandas as pd

# -----------------------------------------
# Column mappings
# -----------------------------------------
COVERAGE_COL = {
    "HA": "Coverage-HA",
    "NA": "Coverage-NA",
    "MP": "Coverage-M",
    "NP": "Coverage-NP",
    "NS": "Coverage-NS",
    "PA": "Coverage-PA",
    "PB1": "Coverage-PB1",
    "PB2": "Coverage-PB2",
}

COVERAGE_ALIASES = {
    "MP": ["Coverage-MP", "Coverage-M"],
}

FRAMESHIFT_COLS = {
    "HA": ["frameShifts HA1", "frameShifts HA2"],
    "NA": ["frameShifts NA"],
    "MP": ["frameShifts M1", "frameShifts M2", "frameShifts MP1", "frameShifts MP2"],
    "NP": ["frameShifts NP"],
    "NS": ["frameShifts NS"],
    "PA": ["frameShifts PA"],
    "PB1": ["frameShifts PB1"],
    "PB2": ["frameShifts PB2"],
}

MIXED_COLS = {
    "HA": ["Nextclade Mixed Sites HA1", "Nextclade Mixed Sites HA2"],
    "NA": ["Nextclade Mixed Sites NA"],
    "MP": [
        "Nextclade Mixed Sites M1",
        "Nextclade Mixed Sites M2",
        "Nextclade Mixed Sites MP1",
        "Nextclade Mixed Sites MP2",
    ],
    "NP": ["Nextclade Mixed Sites NP"],
    "NS": ["Nextclade Mixed Sites NS"],
    "PA": ["Nextclade Mixed Sites PA"],
    "PB1": ["Nextclade Mixed Sites PB1"],
    "PB2": ["Nextclade Mixed Sites PB2"],
}

SEGMENT_ORDER = ["HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2"]


# -----------------------------------------
# Core summarisation logic
# -----------------------------------------
MISSING_TOKENS = {"", "NA", "NAN", "NONE", "NULL", "N/A", "-"}

SUBTYPE_RESULT_MAP = {
    "H3N2": "A/H3N2",
    "H1N1": "A/H1N1",
    "VICVIC": "B/Victoria",
    "VIC": "B/Victoria",
    "YAMYAM": "B/Yamagata",
    "YAM": "B/Yamagata",
}


def is_missing_value(value) -> bool:
    if value is None or pd.isna(value):
        return True
    return str(value).strip().upper() in MISSING_TOKENS


def get_coverage_value(row: pd.Series, seg: str):
    for col in COVERAGE_ALIASES.get(seg, [COVERAGE_COL[seg]]):
        value = row.get(col)
        if not is_missing_value(value):
            return value
    return row.get(COVERAGE_COL[seg])


def has_low_or_missing_coverage(row: pd.Series, seg: str) -> bool:
    cov_raw = get_coverage_value(row, seg)
    if is_missing_value(cov_raw):
        return True

    cov_val = pd.to_numeric(cov_raw, errors="coerce")
    if pd.isna(cov_val):
        return True

    return float(cov_val) < 80


def qc_summary(row: pd.Series) -> str:
    """Return QC summary string for one row."""
    segments_out = []
    for seg in SEGMENT_ORDER:
        issues = []

        # Frameshift: any present value not equal to "No frameShifts".
        for col in FRAMESHIFT_COLS[seg]:
            val = row.get(col)
            if is_missing_value(val):
                continue
            if str(val).strip().lower() != "no frameshifts":
                issues.append("FS")
                break

        # Low coverage includes absent, non-numeric, zero, and <80% coverage.
        if has_low_or_missing_coverage(row, seg):
            issues.append("LC")

        # Mixed sites: coerce each to numeric and sum.
        ms_sum = 0.0
        for col in MIXED_COLS[seg]:
            v = pd.to_numeric(row.get(col), errors="coerce")
            if pd.notna(v):
                ms_sum += float(v)
        if ms_sum > 3:
            issues.append("MS")

        if issues:
            segments_out.append(f"{seg}:{','.join(sorted(issues))}")

    return "|".join(segments_out)


def has_minimum_result_coverage(row: pd.Series) -> bool:
    for seg in ("HA", "NA"):
        cov_raw = get_coverage_value(row, seg)
        if is_missing_value(cov_raw):
            return False
        cov_val = pd.to_numeric(cov_raw, errors="coerce")
        if pd.isna(cov_val) or float(cov_val) < 30:
            return False
    return True


def strict_sekvens_resultat(row: pd.Series) -> str:
    subtype = str(row.get("Subtype", "")).strip()
    if not subtype or subtype.upper() in MISSING_TOKENS:
        return "NA"

    if not has_minimum_result_coverage(row):
        return "NA"

    return SUBTYPE_RESULT_MAP.get(subtype, subtype)


def process_file(in_csv: Path, out_csv: Path) -> None:
    # Read and normalize blanks to <NA>
    df = pd.read_csv(in_csv, low_memory=False)
    obj_cols = df.select_dtypes(include="object").columns
    if len(obj_cols):
        # strip spaces and normalize common empty-like tokens to NA
        df[obj_cols] = df[obj_cols].apply(lambda s: s.str.replace(r"[\u00A0\u200B\uFEFF]", " ", regex=True).str.strip())
        df[obj_cols] = df[obj_cols].replace(to_replace=r"(?i)^(?:na|nan|none|null|n/?a|-)?$", value=pd.NA, regex=True)
        df = df.replace(r"^\s*$", pd.NA, regex=True)

    # Build QC summary
    df["NGS_QC_Sum"] = df.apply(qc_summary, axis=1)

    # If no issues, make it an empty string instead of "NA"
    df["NGS_QC_Sum"] = df["NGS_QC_Sum"].replace(r"^\s*$", "", regex=True)

    # GISAID comment: "Review" if any issues, else empty string
    df["GISAID_Comment"] = df["NGS_QC_Sum"].apply(lambda x: "Review" if str(x).strip() else "")

    # Only call a sequence result when HA and NA coverage are both at least 30%.
    if "Subtype" in df.columns:
        df["Sekvens_Resultat"] = df.apply(strict_sekvens_resultat, axis=1)

    # Write with NA shown explicitly
    df.to_csv(out_csv, index=False, na_rep="NA")
    print(f"Wrote processed file to {out_csv}")  # noqa: T201


# -----------------------------------------
# CLI
# -----------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Add NGS_QC_Sum and GISAID_Comment columns to influenza QC CSV files")
    parser.add_argument("input", type=Path, help="Input CSV file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output CSV file (default: <input>_processed.csv)",
    )
    args = parser.parse_args()

    out_path = args.output or args.input.with_name(f"{args.input.stem}_processed.csv")
    process_file(args.input, out_path)


if __name__ == "__main__":
    main()

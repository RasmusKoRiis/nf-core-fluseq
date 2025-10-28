#!/usr/bin/env python3
import os, glob
import pandas as pd

# Exclude pipeline outputs to avoid feedback loops
EXCLUDE = {"merged_report.csv", "fluseq_merged_report.csv", "id_map.csv"}
csv_files = [f for f in glob.glob("*.csv") if os.path.basename(f) not in EXCLUDE]

KEY_CANDIDATES = ["Sample", "SampleID", "SequenceID", "sample_id", "id", "ID", "Name"]

def infer_sample_from_filename(path: str) -> str:
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    # UID is prefix before first '|' or '_'
    for delim in ("|", "_"):
        if delim in stem:
            return stem.split(delim, 1)[0]
    return stem

def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # drop unnamed index columns and trim headers
    df = df.loc[:, ~df.columns.str.match(r"^Unnamed(:\s*\d+)?$")]
    df.columns = [c.strip() for c in df.columns]
    # all strings, strip whitespace
    return df.apply(lambda s: s.astype(str).str.strip())

frames = []
for f in csv_files:
    try:
        df = pd.read_csv(f, dtype=str, keep_default_na=False)
    except Exception:
        continue
    if df is None or df.empty:
        continue

    df = normalize_df(df)

    # Ensure we have a 'Sample' column
    sample_col = next((c for c in KEY_CANDIDATES if c in df.columns), None)
    if sample_col is None:
        df["Sample"] = infer_sample_from_filename(f)
    elif sample_col != "Sample":
        df = df.rename(columns={sample_col: "Sample"})

    df["Sample"] = df["Sample"].fillna("").astype(str).str.strip()
    df = df[df["Sample"] != ""]
    if not df.empty:
        frames.append(df)

# If nothing to merge, output a safe skeleton
if not frames:
    pd.DataFrame(columns=["Sample"]).to_csv("merged_report.csv", index=False)
    raise SystemExit(0)

merged = pd.concat(frames, ignore_index=True, sort=False)

# For each column per sample, keep the first non-empty value
def first_non_empty(series: pd.Series):
    for x in series:
        if pd.notna(x) and str(x) != "":
            return x
    return ""

agg = {col: first_non_empty for col in merged.columns if col != "Sample"}
merged = merged.groupby("Sample", as_index=False).agg(agg)

# Nice column ordering: Sample, then coverages (if any), then the rest sorted
front = ["Sample"]
priority = [c for c in merged.columns if c.startswith("Coverage-")]
rest = sorted([c for c in merged.columns if c not in front + priority])
merged = merged[front + priority + rest]

merged.to_csv("merged_report.csv", index=False)

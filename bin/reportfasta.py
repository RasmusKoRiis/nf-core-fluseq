#!/usr/bin/env python3
import os, glob
import pandas as pd
import numpy as np

# -------- inputs / exclusions --------
EXCLUDE = {"merged_report.csv", "fluseq_merged_report.csv", "id_map.csv"}
csv_files = [f for f in glob.glob("*.csv") if os.path.basename(f) not in EXCLUDE]

KEY_CANDIDATES = ["Sample", "SampleID", "SequenceID", "sample_id", "id", "ID", "Name"]

def infer_sample_from_filename(path: str) -> str:
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
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
        df = pd.read_csv(f, dtype=str, keep_default_na=False, low_memory=False)
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

# NIPH tweak: replace "!" with "-" in Sample
merged["Sample"] = merged["Sample"].astype(str).str.replace("!", "-", regex=False)

# -------- de-dup: prefer real values over NA/blank --------
def first_best(series: pd.Series):
    for x in series:
        s = str(x).strip()
        if s and s.upper() not in {"NA", "NAN", "NONE"}:
            return s
    for x in series:
        s = str(x).strip()
        if s:
            return s
    return ""

agg = {col: first_best for col in merged.columns if col != "Sample"}
merged = merged.groupby("Sample", as_index=False).agg(agg)

# -------- ensure required cols & compute DR_* + Sekvens_Resultat --------
def ensure_column(df: pd.DataFrame, col: str, fill="NA"):
    if col not in df.columns:
        df[col] = fill

# Columns referenced by DR and QC logic
required = [
    "Subtype",
    "M2 inhibtion mutations",
    "NA inhibtion mutations",
    "PA inhibtion mutations",
    # QC/coverage columns may not always exist but we’ll normalize them later
]
for c in required:
    ensure_column(merged, c, fill="NA")

# Vectorized resistance classification
def classify_res(df: pd.DataFrame, src_col: str, ok_code: str) -> pd.Series:
    s = df[src_col].fillna("NA").astype(str)
    return np.where(
        s.eq("NA"),
        "NA",
        np.where(s.str.contains("No matching mutations", na=False), ok_code, "Review")
    )

# Normalize mutation columns
for src in ["M2 inhibtion mutations", "NA inhibtion mutations", "PA inhibtion mutations"]:
    merged[src] = merged[src].fillna("NA").astype(str)

# DR_Res_* outputs
merged["DR_Res_Adamantine"] = classify_res(merged, "M2 inhibtion mutations", "AANI")
merged["DR_Res_Oseltamivir"] = classify_res(merged, "NA inhibtion mutations", "AANI")
merged["DR_Res_Zanamivir"]   = classify_res(merged, "NA inhibtion mutations", "AANI")
merged["DR_Res_Peramivir"]   = classify_res(merged, "NA inhibtion mutations", "AANI")
merged["DR_Res_Laninamivir"] = classify_res(merged, "NA inhibtion mutations", "AANI")
merged["DR_Res_Baloxavir"]   = classify_res(merged, "PA inhibtion mutations", "AANS")

# DR mutation detail columns
merged["DR_M2_Mut"] = np.where(
    merged["M2 inhibtion mutations"].eq("NA"),
    "NA",
    np.where(merged["DR_Res_Adamantine"].eq("Review"),
             merged["M2 inhibtion mutations"],
             "No Mutations")
)
any_review_na = (
    merged["DR_Res_Oseltamivir"].eq("Review") |
    merged["DR_Res_Zanamivir"].eq("Review")   |
    merged["DR_Res_Peramivir"].eq("Review")   |
    merged["DR_Res_Laninamivir"].eq("Review")
)
merged["DR_NA_Mut"] = np.where(
    merged["NA inhibtion mutations"].eq("NA"),
    "NA",
    np.where(any_review_na, merged["NA inhibtion mutations"], "No Mutations")
)
merged["DR_PA_Mut"] = np.where(
    merged["PA inhibtion mutations"].eq("NA"),
    "NA",
    np.where(merged["DR_Res_Baloxavir"].eq("Review"),
             merged["PA inhibtion mutations"],
             "No Mutations")
)

# Sekvens_Resultat from Subtype
_sub_map = {
    "H3N2": "A/H3N2",
    "H1N1": "A/H1N1",
    "VICVIC": "B/Victoria",
    "VIC": "B/Victoria",
    "YAMYAM": "B/Yamagata",
    "YAM": "B/Yamagata",
}
merged["Subtype"] = merged["Subtype"].astype(str)
merged["Sekvens_Resultat"] = merged["Subtype"].map(_sub_map).fillna(merged["Subtype"])

# Drop any “mammalian” columns if present
mammal_cols = list(merged.filter(regex="mammalian", axis=1).columns)
if mammal_cols:
    merged = merged.drop(columns=mammal_cols)

# -------- numeric cleanup (coverage etc.) --------
coverage_cols = [c for c in merged.columns if c.startswith("Coverage")]
for c in coverage_cols:
    merged[c] = pd.to_numeric(merged[c], errors="coerce").round(2)
if coverage_cols:
    merged[coverage_cols] = merged[coverage_cols].where(~merged[coverage_cols].isna(), other="NA")

# IRMA_noise if present → numeric round(5)
if "IRMA_noise" in merged.columns:
    merged["IRMA_noise"] = pd.to_numeric(merged["IRMA_noise"], errors="coerce").round(5)

# Any other numeric cols → round(5)
num_cols = merged.select_dtypes(include="number").columns.tolist()
other_numeric = [c for c in num_cols if c not in coverage_cols]
if other_numeric:
    merged[other_numeric] = merged[other_numeric].round(5)

# Normalize blanks to 'NA'
merged = merged.replace(r"^\s*$", pd.NA, regex=True).fillna("NA")

# -------- final column ordering --------
front = ["Sample", "Sekvens_Resultat"]
priority = [c for c in merged.columns if c.startswith("Coverage-")]
# Put DR_* columns together next
dr_cols = [c for c in merged.columns if c.startswith("DR_")]
others = [c for c in merged.columns if c not in (front + priority + dr_cols)]
merged = merged[front + priority + dr_cols + sorted(others)]

merged.to_csv("merged_report.csv", index=False)

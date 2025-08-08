#!/usr/bin/env python3
import sys
import glob
import pandas as pd
import numpy as np

# ────────────────────────── helpers ──────────────────────────
def ensure_column(df: pd.DataFrame, column_name: str, fill='NA'):
    """Ensure a column exists; if not, create it filled with `fill`."""
    if column_name not in df.columns:
        df[column_name] = fill

def classify_res(df: pd.DataFrame, src_col: str, ok_code: str) -> pd.Series:
    """
    Vectorized resistance classification:
    - 'NA' -> 'NA'
    - contains 'No matching mutations' -> ok_code
    - else -> 'Review'
    """
    s = df[src_col].fillna('NA').astype(str)
    return np.where(
        s.eq('NA'),
        'NA',
        np.where(s.str.contains('No matching mutations', na=False), ok_code, 'Review')
    )

# ────────────────────────── I/O ─────────────────────────────
# CLI
samplesheet = sys.argv[1]

# Load sample sheet (TSV)
samplesheet_df = pd.read_csv(samplesheet, sep='\t', dtype=str)
samplesheet_df.rename(columns={'SequenceID': 'Sample'}, inplace=True)
# Drop Barcode if present
if 'Barcode' in samplesheet_df.columns:
    samplesheet_df = samplesheet_df.drop(columns=['Barcode'])

# Collect all CSVs in CWD (data reports)
csv_files = glob.glob('*.csv')

# Read all CSVs as strings to avoid mixed dtypes; convert later where needed
dataframes = []
for f in csv_files:
    try:
        df = pd.read_csv(f, dtype=str, low_memory=False)
        dataframes.append(df)
    except Exception as e:
        # Skip unreadable CSVs but keep going
        sys.stderr.write(f"[WARN] Skipping {f}: {e}\n")

if dataframes:
    merged_data = pd.concat(dataframes, ignore_index=True)
else:
    # No CSVs found → start empty
    merged_data = pd.DataFrame()

# ────────────────────────── harmonize ───────────────────────
# If Sample column is missing entirely in data, create it to proceed
ensure_column(merged_data, 'Sample')

# Keep first occurrence per Sample (same as groupby(...).first, but faster here)
# If Sample has NA, drop them before de-dup
merged_data = merged_data.dropna(subset=['Sample']).copy()
merged_data = merged_data.groupby('Sample', as_index=False, sort=False).first()

# NIPH specific: replace "!" with "-" in Sample
merged_data['Sample'] = merged_data['Sample'].astype(str).str.replace('!', '-', regex=False)

# ────────────────────────── ensure required columns ─────────
required_columns = [
    'Sample', 'Sekvens_Resultat', 'Coverage-HA', 'Coverage-M', 'Coverage-NA', 'Coverage-NP',
    'Coverage-NS', 'Coverage-PA', 'Coverage-PB1', 'Coverage-PB2', 'DEPTH_HA', 'DEPTH_MP',
    'DEPTH_NA', 'DEPTH_NP', 'DEPTH_NS', 'DEPTH_PA', 'DEPTH_PB1', 'DEPTH_PB2', 'DR_M2_Mut',
    'DR_NA_Mut', 'DR_PA_Mut', 'DR_Res_Adamantine', 'DR_Res_Baloxavir', 'DR_Res_Oseltamivir',
    'DR_Res_Peramivir', 'DR_Res_Zanamivir', 'HA1 Differences human', 'HA1 Differences human_vaccine',
    'HA2 Differences human', 'HA2 Differences human_vaccine', 'IRMA_altmatch', 'IRMA_chimeric',
    'IRMA_failQC', 'IRMA_initial', 'IRMA_match', 'IRMA_nomatch', 'IRMA_passQC', 'M1 Differences human',
    'M2 Differences human', 'M2 Differences inhibition_human', 'M2 inhibtion mutations',
    'NA Differences human', 'NA Differences human_vaccine', 'NA Differences inhibition_human',
    'NA inhibtion mutations', 'NP Differences human', 'NS1 Differences human', 'Nextclade QC HA1',
    'Nextclade QC M1', 'Nextclade QC NA', 'Nextclade QC NP', 'Nextclade QC NS', 'Nextclade QC PA',
    'Nextclade QC PB1', 'Nextclade QC PB2', 'PA Differences human', 'PA Differences inhibition_human',
    'PA inhibtion mutations', 'PB1 Differences human', 'PB2 Differences human', 'SigPep Differences human',
    'Subtype', 'aaDeletions HA1', 'aaDeletions M1', 'aaDeletions M2', 'aaDeletions NA', 'aaDeletions NP',
    'aaDeletions NS', 'aaDeletions PA', 'aaDeletions PB1', 'aaDeletions PB2', 'aaInsertions HA1',
    'aaInsertions M1', 'aaInsertions M2', 'aaInsertions NA', 'aaInsertions NP', 'aaInsertions NS',
    'aaInsertions PA', 'aaInsertions PB1', 'aaInsertions PB2', 'clade', 'clade NA', 'frameShifts HA1',
    'frameShifts M1', 'frameShifts M2', 'frameShifts NA', 'frameShifts NP', 'frameShifts NS',
    'frameShifts PA', 'frameShifts PB1', 'frameShifts PB2', 'glycosylation', 'subclade'
]
for col in required_columns:
    ensure_column(merged_data, col)

# Also referenced later but not in required_columns in original script
ensure_column(merged_data, 'IRMA_noise')

# ────────────────────────── RESISTANCE (vectorized) ────────
# Normalize the three source columns to strings with 'NA' for missing
for src in ['M2 inhibtion mutations', 'NA inhibtion mutations', 'PA inhibtion mutations']:
    merged_data[src] = merged_data[src].fillna('NA').astype(str)

# DR_Res_* columns
merged_data['DR_Res_Adamantine'] = classify_res(merged_data, 'M2 inhibtion mutations', 'AANI')
merged_data['DR_Res_Oseltamivir'] = classify_res(merged_data, 'NA inhibtion mutations', 'AANI')
merged_data['DR_Res_Zanamivir']   = classify_res(merged_data, 'NA inhibtion mutations', 'AANI')
merged_data['DR_Res_Peramivir']   = classify_res(merged_data, 'NA inhibtion mutations', 'AANI')
merged_data['DR_Res_Laninamivir'] = classify_res(merged_data, 'NA inhibtion mutations', 'AANI')
merged_data['DR_Res_Baloxavir']   = classify_res(merged_data, 'PA inhibtion mutations', 'AANS')

# DR_M2_Mut
merged_data['DR_M2_Mut'] = np.where(
    merged_data['M2 inhibtion mutations'].eq('NA'),
    'NA',
    np.where(merged_data['DR_Res_Adamantine'].eq('Review'),
             merged_data['M2 inhibtion mutations'],
             'No Mutations')
)

# DR_PA_Mut
merged_data['DR_PA_Mut'] = np.where(
    merged_data['PA inhibtion mutations'].eq('NA'),
    'NA',
    np.where(merged_data['DR_Res_Baloxavir'].eq('Review'),
             merged_data['PA inhibtion mutations'],
             'No Mutations')
)

# DR_NA_Mut (any NA drug flagged Review → include mutation list)
any_review_na = (
    merged_data['DR_Res_Oseltamivir'].eq('Review') |
    merged_data['DR_Res_Zanamivir'].eq('Review')   |
    merged_data['DR_Res_Peramivir'].eq('Review')   |
    merged_data['DR_Res_Laninamivir'].eq('Review')
)
merged_data['DR_NA_Mut'] = np.where(
    merged_data['NA inhibtion mutations'].eq('NA'),
    'NA',
    np.where(any_review_na, merged_data['NA inhibtion mutations'], 'No Mutations')
)

# ────────────────────────── SUBTYPE column ─────────────────
# Map as per original logic; otherwise keep original value
_sub_map = {
    'H3N2': 'A/H3N2',
    'H1N1': 'A/H1N1',
    'VICVIC': 'B/Victoria',
    'VIC': 'B/Victoria',
    'YAMYAM': 'B/Yamagata',
    'YAM': 'B/Yamagata',
}
merged_data['Subtype'] = merged_data['Subtype'].astype(str)
mapped = merged_data['Subtype'].map(_sub_map)
merged_data['Sekvens_Resultat'] = mapped.fillna(merged_data['Subtype'])

# ────────────────────────── Remove 'mammalian' cols ────────
mammalian_cols = list(merged_data.filter(regex='mammalian', axis=1).columns)
if mammalian_cols:
    merged_data = merged_data.drop(columns=mammalian_cols)

# ────────────────────────── Column order & new samples ─────
# Put Sample + Sekvens_Resultat first, rest after (alphabetical stable order)
front = ['Sample', 'Sekvens_Resultat']
others = [c for c in merged_data.columns if c not in front]
merged_data = merged_data[front + sorted(others)]

# Add samples present in samplesheet but not yet in merged_data
new_samples = samplesheet_df[~samplesheet_df['Sample'].isin(merged_data['Sample'])]
merged_data = pd.concat([merged_data, new_samples], ignore_index=True)

# Ensure all columns from both frames are present and ordered: keep current order + any extras from new_samples
extra_cols = [c for c in new_samples.columns if c not in merged_data.columns]
for c in extra_cols:
    merged_data[c] = 'NA'
# Reindex to front + sorted others again to keep consistent layout
others = [c for c in merged_data.columns if c not in front]
merged_data = merged_data[front + sorted(others)]

# ────────────────────────── Numeric handling & rounding ────
# Convert IRMA_noise to numeric safely and round 5
merged_data["IRMA_noise"] = pd.to_numeric(merged_data.get("IRMA_noise", pd.Series(index=merged_data.index)), errors='coerce').round(5)

# Coverage columns: numeric → round(2) → fill 'NA' for missing
coverage_columns = [col for col in merged_data.columns if "Coverage" in col]
for col in coverage_columns:
    merged_data[col] = pd.to_numeric(merged_data[col], errors='coerce').round(2)
# After rounding, replace NaN with 'NA' strings as in original
merged_data[coverage_columns] = merged_data[coverage_columns].where(~merged_data[coverage_columns].isna(), other='NA')

# Round any other numeric cols to 5 decimals (keeps coverage at 2 as already done)
numeric_cols = merged_data.select_dtypes(include='number').columns.tolist()
other_numeric = [c for c in numeric_cols if c not in coverage_columns]
if other_numeric:
    merged_data[other_numeric] = merged_data[other_numeric].round(5)

# ── FINAL SWEEP: normalize empties across the whole DF ──
# 1) Trim whitespace in all string columns
obj_cols = merged_data.select_dtypes(include='object').columns
merged_data[obj_cols] = merged_data[obj_cols].apply(lambda s: s.str.strip())

# 2) Turn common empty-like tokens into real <NA>
merged_data[obj_cols] = merged_data[obj_cols].replace(
    to_replace=r'^(?i)(?:na|nan|none|null|n/?a|-)?$',  # case-insensitive
    value=pd.NA,
    regex=True
)

# 3) Also convert pure blanks to <NA>, then fill everything with "NA"
merged_data = merged_data.replace(r'^\s*$', pd.NA, regex=True).fillna("NA")


# ────────────────────────── Write ──────────────────────────
merged_data.to_csv('merged_report.csv', index=False)

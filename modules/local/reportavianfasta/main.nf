process REPORTAVIANFASTA {
    label 'process_single'

    // Python helper scripts live in /project-bin on the container (mounted from repo bin/)
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    /*
      Use `path` so Nextflow STAGES the files into CWD.
      These can be lists of paths (from .toList()) — Nextflow will copy them.
      For context, the source data come from N:\Virologi\BioNumerics\RARI\Avian Influensa SQLite.
    */
    input:
    path subtype                   // list or single; subtype summary tables
    path coverage                  // list; coverage QC per sample
    path mutation_mamalian         // list; mammalian mutation annotations
    path lookup_1                  // list; lookup tables (part 1)
    path lookup_2                  // list; lookup tables (part 2)
    path nextclade_summary_ha      // list; Nextclade summary for HA
    path nextclade_sample          // list; per-sample Nextclade CSVs
    path id_map                    // single TSV (SampleID \t OriginalName)
    val  runid                     // used to name outputs / metadata
    val  release_version
    val  filtered_fasta            // ignored here (passed along for parity with other reports)
    val  seq_instrument
    val  samplesheet               // ignored here (useful context when debugging)
    path genin2                    // staged so report scripts can pick up GENIN2 outputs
    path flumut                    // staged FLUMUT outputs

    output:
    path("${runid}.csv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
"""
set -euo pipefail

echo "[REPORTHUMANFASTA] staged CSVs:"
ls -1 *.csv || true

# STEP 1: merge all staged CSVs (reportfasta.py handles de-duplication logic)
python /project-bin/reportfasta.py

# STEP 2: enrich merged report with ID mapping + metadata before QC
export RUNID='${runid}'
export SEQINST='${seq_instrument}'
export RELVER='${release_version}'
export IDMAP='${id_map}'

python - <<'PY'
import os
import pandas as pd
from datetime import date

runid   = os.environ.get("RUNID","unknown")
seqinst = os.environ.get("SEQINST","unknown")
relver  = os.environ.get("RELVER","unknown")
idmap_p = os.environ["IDMAP"]

# --- helpers (regex-free split on '|' or '_') ---
def cut_prefix(x: str) -> str:
    x = str(x)
    for d in ('|','_'):
        i = x.find(d)
        if i != -1:
            return x[:i]
    return x

def norm_uid_series(s: pd.Series) -> pd.Series:
    return s.astype(str).map(cut_prefix)

def dedup_on_key(df: pd.DataFrame, key: str) -> pd.DataFrame:
    if df.empty or key not in df.columns: 
        return df
    df2  = df.copy()
    cols = [c for c in df2.columns if c != key]
    if not cols:
        return df2.drop_duplicates(subset=[key], keep="first")
    non_empty = (df2[cols].notna()) & (df2[cols].astype(str) != "")
    df2["__score"] = non_empty.sum(axis=1)
    df2 = df2.sort_values("__score", ascending=False).drop(columns="__score")
    return df2.drop_duplicates(subset=[key], keep="first")

# --- id_map clean (no backticks / no  examples) ---
idmap = pd.read_csv(idmap_p, sep="\t", dtype=str, keep_default_na=False)
idmap.columns = ["SampleID","OriginalName"]
idmap["SampleID"]     = idmap["SampleID"].astype(str).str.strip()
idmap["OriginalName"] = idmap["OriginalName"].astype(str).str.strip()
mask = idmap["SampleID"].str.len().gt(0) & idmap["SampleID"].str.isalnum()
idmap = idmap[mask].drop_duplicates(subset=["SampleID"], keep="first")

# --- load merged ---
try:
    merged = pd.read_csv("merged_report.csv", dtype=str, keep_default_na=False)
except Exception:
    merged = pd.DataFrame(columns=["Sample"])

uids  = set(idmap["SampleID"])
names = set(idmap["OriginalName"])

# choose best join key: UID exact -> UID normalized -> OriginalName
best_key, best_hits, mode = None, -1, "uid_exact"
for c in merged.columns:
    hits = merged[c].astype(str).isin(uids).sum()
    if hits > best_hits:
        best_hits, best_key, mode = hits, c, "uid_exact"

if best_hits == 0 and not merged.empty:
    for c in merged.columns:
        hits = norm_uid_series(merged[c]).isin(uids).sum()
        if hits > best_hits:
            merged["__UID"] = norm_uid_series(merged[c])
            best_hits, best_key, mode = hits, "__UID", "uid_norm"

if best_hits == 0 and not merged.empty:
    for c in merged.columns:
        hits = merged[c].astype(str).isin(names).sum()
        if hits > best_hits:
            best_hits, best_key, mode = hits, c, "name_exact"

# de-dup merged on chosen key (prefer "richer" rows)
if best_key:
    merged = dedup_on_key(merged, best_key)

# left-join so every SampleID appears once
if best_key:
    if mode in ("uid_exact","uid_norm"):
        out = idmap.merge(merged, how="left", left_on="SampleID", right_on=best_key)
    else:
        out = idmap.merge(merged, how="left", left_on="OriginalName", right_on=best_key)
    if "__UID" in out.columns:
        out = out.drop(columns=["__UID"])
else:
    out = idmap.copy()

# final guards: unique rows, drop mirror 'Sample'
out = out.drop_duplicates()
out = out.drop_duplicates(subset=["SampleID"], keep="first")
if "Sample" in out.columns:
    same_raw  = out["Sample"].astype(str).eq(out["SampleID"]).all()
    same_norm = out["Sample"].astype(str).map(cut_prefix).eq(out["SampleID"]).all()
    if same_raw or same_norm:
        out = out.drop(columns=["Sample"])

# meta columns
out["RunID"]           = runid
out["Instrument ID"]   = seqinst
out["Date"]            = date.today().isoformat()
out["Release Version"] = relver

# order
front = ["SampleID","OriginalName","RunID","Instrument ID","Date","Release Version"]
rest  = [c for c in out.columns if c not in front]
out = out[front + rest]

# hand off to QC calc
out.to_csv(f"{runid}_qc_input.csv", index=False)
PY

# STEP 3: run QC calculation, producing the final CSV
python /project-bin/report_QC_calculation.py ${runid}_qc_input.csv -o ${runid}.csv

# STEP 4: post-process HA differences column for readability on downstream systems
python - <<'PY'
import os
import pandas as pd

runid = os.environ.get("RUNID", "unknown")
fname = f"{runid}.csv"

col  = "HA Differences mamailian"
col1 = f"{col}_1"
col2 = f"{col}_2"
MAX_LEN = 148

def split_safely(s: str, max_len: int = MAX_LEN):
    s = "" if s is None else str(s)
    if s == "":
        return "", ""

    # Find the last ';' such that the left part length <= max_len
    cut = s.rfind(";", 0, max_len)  # searches in s[0:max_len)

    if cut != -1:
        left = s[:cut + 1]                  # include ';' so _1 ends with ';'
        right = s[cut + 1:].lstrip()        # rest to _2
        return left, right

    # Edge case: no ';' before max_len (can't split without cutting a mutation)
    # If it fits, we can still put it in _1; optionally add ';' if there's room.
    if len(s) < max_len and not s.endswith(";"):
        return s + ";", ""
    if len(s) <= max_len and s.endswith(";"):
        return s, ""

    # Otherwise, keep everything in _2
    return "", s

df = pd.read_csv(fname, dtype=str, keep_default_na=False)

# Ensure columns exist even if source column is missing
if col not in df.columns:
    df[col] = ""

parts = df[col].apply(lambda x: pd.Series(split_safely(x)))
df[col1] = parts[0]
df[col2] = parts[1]

# Place new columns right after the original column
cols = [c for c in df.columns if c not in (col1, col2)]
insert_at = cols.index(col) + 1 if col in cols else len(cols)
cols = cols[:insert_at] + [col1, col2] + cols[insert_at:]
df = df[cols]

tmp = fname + ".tmp"
df.to_csv(tmp, index=False)
os.replace(tmp, fname)
PY




"""
}

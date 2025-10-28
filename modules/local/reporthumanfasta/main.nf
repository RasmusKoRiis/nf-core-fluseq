process REPORTHUMANFASTA {
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    /*
      Use `path` so Nextflow STAGES the files into CWD.
      These can be lists of paths (from .toList()) — Nextflow will copy them.
    */
    input:
    path subtype                   // list or single
    path coverage                  // list
    path mutation_human            // list
    path mutation_inhibtion        // list
    path lookup                    // list
    path nextclade_summary_ha      // list
    path nextclade_sample          // list
    path mutation_vaccine          // list
    path id_map                    // single TSV (SampleID \t OriginalName)
    val  runid
    val  release_version
    val  filtered_fasta            // ignored here
    val  seq_instrument
    val  samplesheet               // ignored here

    output:
    path("${runid}.csv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
"""
set -euo pipefail

echo "[REPORTHUMANFASTA] staged CSVs:"
ls -1 *.csv || true

# 1) Merge staged CSVs into merged_report.csv (dedup by Sample, keep first non-empty per column)
python /project-bin/reportfasta.py

# 2) Join OriginalName + add meta, robust de-dup (no .str.split)
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

# 3) QC calculation -> FINAL report (adds DR_* + Sekvens_Resultat)
python /project-bin/report_QC_calculation.py ${runid}_qc_input.csv -o ${runid}.csv
"""
}

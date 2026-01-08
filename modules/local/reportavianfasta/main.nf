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
    path mutation_vaccine          // list; vaccine mutation annotations
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

# STEP 4: post-process Difference columns for readability on downstream systems
python - <<'PY'
import os
import pandas as pd

runid = os.environ.get("RUNID", "unknown")
fname = f"{runid}.csv"

df = pd.read_csv(fname, dtype=str, keep_default_na=False)

CHUNK_COUNT = 3
MAX_LEN = 140
TARGET_SUFFIX = "Differences mamailian"

def normalize_string(val: str) -> str:
    if val is None:
        return ""
    text = str(val).strip()
    if text and not text.endswith(";"):
        text += ";"
    return text

def chunk_mutations(val: str):
    text = normalize_string(val)
    if not text:
        return [""] * CHUNK_COUNT

    tokens = [t.strip() for t in text.split(";") if t.strip()]
    token_chunks = [f"{t};" for t in tokens]

    chunks = []
    idx = 0
    while len(chunks) < CHUNK_COUNT and idx < len(token_chunks):
        chunk = ""
        while idx < len(token_chunks):
            piece = token_chunks[idx]
            if not chunk:
                chunk = piece
                idx += 1
                if len(chunk) > MAX_LEN:
                    break
                continue
            if len(chunk) + len(piece) <= MAX_LEN:
                chunk += piece
                idx += 1
            else:
                break
        chunks.append(chunk)

    if idx < len(token_chunks):
        tail = "".join(token_chunks[idx:])
        if chunks:
            chunks[-1] = (chunks[-1] + tail)
        else:
            chunks.append(tail)

    while len(chunks) < CHUNK_COUNT:
        chunks.append("")

    return chunks[:CHUNK_COUNT]

target_columns = [c for c in list(df.columns) if c.endswith(TARGET_SUFFIX)]

for col in target_columns:
    parts = pd.DataFrame(df[col].apply(chunk_mutations).tolist(),
                         columns=[f"{col}_{i+1}" for i in range(CHUNK_COUNT)])
    insert_at = df.columns.get_loc(col) + 1 if col in df.columns else len(df.columns)
    for offset, new_col in enumerate(parts.columns):
        if new_col in df.columns:
            df.drop(columns=[new_col], inplace=True)
        df.insert(insert_at + offset, new_col, parts[new_col])

tmp = fname + ".tmp"
df.to_csv(tmp, index=False)
os.replace(tmp, fname)
PY




"""
}

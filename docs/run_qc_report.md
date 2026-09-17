# Human FASTQ HTML QC report

The `REPORT_QC_HTML` module runs automatically after `REPORTHUMAN` in the
`human-fastq` workflow. It reads the final CSV directly from the report output
channel and publishes these files beside it:

- `reporthuman/<runid>_qc.html`: a self-contained, searchable run assessment;
- `reporthuman/<runid>_qc.json`: the same assessment and sample review data in a
  versioned JSON structure.

Open the HTML file in a browser; no server or internet connection is needed.
Filter by QC assessment, subtype or free text, switch segment cells between
coverage and IRMA read counts, and expand sample details. Print/save PDF uses
the currently displayed samples and segment metric; reset filters to print
all samples. Run totals always refer to the full CSV. The report remains
readable with JavaScript disabled.

## Assessment

The top cards reproduce the source `NGS_QC_Sum` counts: an explicitly empty
cell has no reported flags, a populated summary prompts review, and a missing
column or NA-like cell is unavailable. The sample assessment also checks:

- coverage below the existing report threshold of 80%;
- absent, non-numeric, non-finite or out-of-range coverage;
- missing NGS summary;
- non-good or unavailable Nextclade QC for the ten expected protein columns;
- a populated source `GISAID_Comment`.

The eight segments are HA, NA, MP, NP, NS, PA, PB1 and PB2. `Coverage-MP` takes
precedence over `Coverage-M`; a missing MP value permits the M alias. Missing
values remain distinct from zero. Medians exclude unavailable measurements,
and tables show availability denominators. An entirely absent NGS summary
and entirely absent coverage give an unavailable sample assessment.

Any sample requiring review or unavailable assessment, or incomplete or mixed
run metadata, prompts run review. A header-only CSV is not assessable.
“No QC flags detected” is not formal run acceptance: control outcomes, the
expected sample count, and laboratory acceptance rules are not available from
this CSV. Samples omitted from the CSV cannot be detected by this module.

The 80% boundary follows `report_QC_calculation.py`, which currently uses a
fixed threshold independently of `seq_quality_threshold`. The renderer does
not recalculate or modify source mutation flags. It preserves the original CSV.

## Findings and metric definitions

The overview counts reported subtypes and characterisation statuses. VIC and
VICVIC are both displayed as B/Victoria; raw labels remain in JSON. Raw subtype
counts are distinct from coverage-gated `Sekvens_Resultat`, available in sample
details. Subclade distributions distinguish provisional calls, exact rule
matches and calls without assessed match evidence. Incomplete match fractions,
missing defining-state evidence or source review statuses retain a provisional
label. Genetic classification labels are not antigenic phenotype predictions.

`DEPTH_*` means IRMA segment read count, **not mean per-base depth**. Reported
coverage normally measures the non-N percentage of observed consensus sequence;
it does not establish recovery of the full expected segment. Missing mixed-site
and frame-shift inputs can leave the upstream QC summary empty. See the
[data dictionary](report_data_dictionary.md) for source metric limitations.

Sample details preserve reported reassortment conclusions and drug lookup
codes without adding biological interpretations. Similarity-screen flags do
not establish reassortment, and lookup codes are not clinical susceptibility
results.

## Render an existing CSV

From the repository root, with Python 3.9 or newer:

```bash
python bin/report_qc_html.py /private/path/run.csv --outdir /private/path/qc
```

The standalone renderer uses only the Python standard library. The Nextflow
module uses the same digest-pinned Python container as `REPORTHUMAN`. Both the
script and template are staged as inputs, so edits invalidate the task cache.
The new module does not run in FASTA or avian workflows.

JSON includes the generator version, schema version, source filename and SHA-256,
run metadata, assessment rules, warnings, counts, segment summaries and sample
review data. Identical input and renderer produce identical output; the HTML
does not include a render timestamp. Retain the original CSV and pipeline and
reference provenance. Both output files contain the sample identifiers from
the source and should be stored with the run's other reports.

## Verification

```bash
python -m pytest tests/test_report_qc_html.py tests/test_module_metadata.py
```

Tests use invented records and cover QC boundaries, missing values, aliases,
provisional labels, CSV validation, HTML escaping, deterministic output and
actual Nextflow module execution/publication when Nextflow is installed.

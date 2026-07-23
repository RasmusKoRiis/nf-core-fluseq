# Report data dictionary

This document describes the CSV produced by the human report workflow
(`modules/local/reporthuman/main.nf`).  The file is first assembled by
[`bin/report.py`](../bin/report.py), then receives `RunID`, `Instrument ID`,
`Date`, and `Release Version`, and is finally processed by
[`bin/report_QC_calculation.py`](../bin/report_QC_calculation.py).  The avian
and FASTA workflows use closely related report builders; fields not listed
below are pass-through columns from their staged input CSVs and should be
documented by the producing module.

## Conventions

* One row represents one `Sample`.  Duplicate sample rows are collapsed; the
  first non-null value selected by the merge is not necessarily the most recent
  or highest-quality result.
* `NA` means unavailable, not calculated, or not applicable.  It is not a
  negative result.  Empty QC summaries mean that no issue was detected.
* Coverage is a percentage, rounded to two decimals.  Other numeric values are
  rounded to five decimals; `IRMA_noise` is rounded to five decimals.
* “Source” identifies the upstream report/table or the calculation in the
  named script. Exact software/database versions should be read from the run's
  `software_versions.yml` and parameters from `params.json`.

## Columns

| Column | Meaning and units | How it is made | Source | Important weakness / interpretation |
|---|---|---|---|---|
| `Sample` | Sample identifier | Sample-sheet `SequenceID` (or upstream sample key); `!` changed to `-` | Input samplesheet and staged CSVs | Identifier normalisation can merge distinct names; duplicates are collapsed |
| `Sekvens_Resultat` | Human-readable influenza type/subtype | Maps `H3N2`→`A/H3N2`, `H1N1`→`A/H1N1`, `VIC(VIC)`→`B/Victoria`, `YAM(YAM)`→`B/Yamagata`; otherwise retains `Subtype` (QC script only reports it when HA and NA coverage ≥30%) | `Subtype` from Nextclade/subtype analysis | A subtype label is not proof of complete genome quality; low coverage can be reported as `NA` |
| `RunID` | Pipeline run identifier | Constant supplied to the Nextflow process | Workflow parameter | Describes the run, not sample-level provenance |
| `Instrument ID` | Sequencing instrument identifier | Constant supplied to the process | Workflow parameter | May be missing, manually entered, or shared by many samples |
| `Date` | Report generation date (`YYYY-MM-DD`) | Shell `date` at report-process execution | Execution environment clock | Not the collection, sequencing, or analysis date; clock/time-zone errors are possible |
| `Release Version` | Pipeline/reference release label | Constant supplied to the process | Workflow parameter | Does not by itself identify every database/tool version |
| `Coverage-HA`, `Coverage-M`, `Coverage-NA`, `Coverage-NP`, `Coverage-NS`, `Coverage-PA`, `Coverage-PB1`, `Coverage-PB2` | Percent of each influenza segment meeting the coverage criterion | Numeric conversion and rounding in `report.py`; `M` denotes the MP segment | Coverage module output | Threshold and reference depend on upstream configuration; missing/non-numeric values become `NA` |
| `DEPTH_HA`, `DEPTH_MP`, `DEPTH_NA`, `DEPTH_NP`, `DEPTH_NS`, `DEPTH_PA`, `DEPTH_PB1`, `DEPTH_PB2` | Per-segment read depth (depth units) | Passed through from depth analysis | `depth_analysis*` outputs | Mean/median and position handling depend on the producing script; depth alone does not establish correctness |
| `Nextclade QC HA1`, `Nextclade QC M1`, `Nextclade QC NA`, `Nextclade QC NP`, `Nextclade QC NS`, `Nextclade QC PA`, `Nextclade QC PB1`, `Nextclade QC PB2` | Nextclade quality-control result for a segment/protein | Passed through from Nextclade summaries | Nextclade summary output | Tool/reference version and QC thresholds can change; a QC flag is not a diagnosis |
| `Subtype` | Subtype/type reported by upstream analysis | Passed through; used to derive `Sekvens_Resultat` | Nextclade/subtype analysis | Ambiguous or contaminated samples may receive an unstable label |
| `clade`, `clade NA`, `subclade` | Nextclade clade/subclade assignments | Passed through from Nextclade and nomenclature steps | Nextclade and `subclade_nomenclature.py` | Nomenclature is reference-release dependent and may be `NA` for novel/divergent viruses |
| `aaDeletions HA1`, `aaDeletions M1`, `aaDeletions M2`, `aaDeletions NA`, `aaDeletions NP`, `aaDeletions NS`, `aaDeletions PA`, `aaDeletions PB1`, `aaDeletions PB2` | Amino-acid deletion positions | Passed through from Nextclade amino-acid mutation output | Nextclade | Calls are limited to covered/reference-aligned positions; mixed/low-quality bases can hide events |
| `aaInsertions HA1`, `aaInsertions M1`, `aaInsertions M2`, `aaInsertions NA`, `aaInsertions NP`, `aaInsertions NS`, `aaInsertions PA`, `aaInsertions PB1`, `aaInsertions PB2` | Amino-acid insertion positions/sequences | Passed through from Nextclade | Nextclade | Same alignment and coverage limitations as deletion calls |
| `frameShifts HA1`, `frameShifts M1`, `frameShifts M2`, `frameShifts NA`, `frameShifts NP`, `frameShifts NS`, `frameShifts PA`, `frameShifts PB1`, `frameShifts PB2` | Detected coding-frame shifts | Passed through; interpreted by QC as an `FS` issue unless `No frameShifts` | Nextclade | Sequencing errors and partial segments can create false calls; absence is not evidence of intact biology |
| `glycosylation` | Predicted glycosylation-site changes | Passed through from mutation analysis | Mutation/Nextclade output | Prediction depends on complete local amino-acid context and reference; it is not an experimentally measured glycan |
| `HA1 Differences human`, `HA1 Differences human_vaccine`, `HA2 Differences human`, `HA2 Differences human_vaccine`, `M1 Differences human`, `M2 Differences human`, `NA Differences human`, `NA Differences human_vaccine`, `NA Differences inhibition_human`, `NP Differences human`, `NS1 Differences human`, `PA Differences human`, `PA Differences inhibition_human`, `PB1 Differences human`, `PB2 Differences human`, `SigPep Differences human` | Amino-acid differences from the named human reference, vaccine, or inhibition reference; positions/changes as encoded by upstream output | Passed through from mutation comparison tables | `mutation_human.py`, vaccine and inhibition comparison outputs | “Difference” is relative to a reference, not necessarily a functional mutation; reference choice and alignment affect results |
| `M2 inhibtion mutations` | M2 mutations relevant to adamantane resistance | Passed through from inhibition mutation analysis | `mutation_finder.py` / lookup table | Column name contains the historical `inhibtion` typo; interpretation depends on lookup database |
| `NA inhibtion mutations` | NA mutations relevant to neuraminidase-inhibitor resistance | Passed through from inhibition mutation analysis | Mutation finder / resistance lookup | Same reference/database and coverage limitations; typo is retained for compatibility |
| `PA inhibtion mutations` | PA mutations relevant to baloxavir resistance | Passed through from inhibition mutation analysis | Mutation finder / resistance lookup | Same limitations; a detected mutation is not a phenotypic susceptibility measurement |
| `DR_Res_Adamantine` | Adamantane resistance classification | `NA` if source missing; `AANI` if source says “No matching mutations”; otherwise `Review` | Derived in `report.py` | Heuristic string classification; `AANI` is a workflow code and should not be read as susceptibility |
| `DR_Res_Oseltamivir`, `DR_Res_Zanamivir`, `DR_Res_Peramivir`, `DR_Res_Laninamivir` | Neuraminidase-inhibitor resistance classification | Same rule using `NA inhibtion mutations` | Derived in `report.py` | Does not model genotype–phenotype uncertainty, mixtures, or drug-specific evidence |
| `DR_Res_Baloxavir` | Baloxavir resistance classification | `NA`/`AANS`/`Review` rule using `PA inhibtion mutations` | Derived in `report.py` | Heuristic and reference dependent; `AANS` is a workflow code |
| `DR_M2_Mut` | M2 resistance mutation detail | Source missing→`NA`; flagged review→mutation text; otherwise `No Mutations` | Derived in `report.py` | “No Mutations” means no match in the lookup, not proof of susceptibility |
| `DR_NA_Mut` | NA resistance mutation detail | Source missing→`NA`; any NA drug marked `Review`→mutation text; otherwise `No Mutations` | Derived in `report.py` | Collapses four drugs into one detail field and can obscure drug-specific evidence |
| `DR_PA_Mut` | PA resistance mutation detail | Source missing→`NA`; baloxavir review→mutation text; otherwise `No Mutations` | Derived in `report.py` | Lookup-based and not a phenotypic assay |
| `IRMA_altmatch` | IRMA alternative-match indicator/count | Passed through from IRMA report | IRMA | Meaning and coding are IRMA-version dependent; not comparable without that version |
| `IRMA_chimeric` | IRMA chimeric-read/segment flag | Passed through from IRMA | IRMA | Algorithmic flag; may reflect technical artefact or true reassortment |
| `IRMA_failQC`, `IRMA_passQC` | IRMA overall QC fail/pass indicators | Passed through from IRMA | IRMA | The two flags may both be absent/ambiguous; thresholds are tool/configuration dependent |
| `IRMA_initial`, `IRMA_match`, `IRMA_nomatch` | IRMA initial/matched/non-matched segment or reference counts/statuses | Passed through from IRMA | IRMA | Exact semantics depend on IRMA report schema and version |
| `IRMA_noise` | IRMA noise metric | Numeric conversion and rounding to five decimals | IRMA | Scale and recommended cutoff are version/configuration dependent |
| `NGS_QC_Sum` | Compact segment QC summary, e.g. `HA:MS|PB1:LC,FS` | `report_QC_calculation.py`: `LC` for missing/non-numeric/<80% coverage, `FS` for frame shift, `MS` for >3 mixed sites; segments joined with `|` | Coverage, Nextclade frame-shift and mixed-site fields | Thresholds are hard-coded; an empty value means no listed issue, not complete validation |
| `GISAID_Comment` | Submission review suggestion | `Review` when `NGS_QC_Sum` is non-empty; otherwise empty | Derived in `report_QC_calculation.py` | Administrative recommendation only; it does not replace curator review or GISAID validation |

## Provenance and reproducibility

The report is a merge of every readable CSV staged in the report process working
directory.  Therefore, the complete provenance of a value includes the input
file, pipeline parameters, reference datasets, and software versions.  Retain
the run's `samplesheet.valid.csv`, `params.json`, `software_versions.yml`,
Nextflow execution report, and the staged upstream CSVs with the report.

The merge is intentionally permissive: unreadable CSVs are skipped, missing
columns are created as `NA`, and one row per sample is retained.  These choices
make a report robust to incomplete runs but can conceal missing inputs.  A
report should consequently be reviewed together with process logs and QC
outputs, especially when a high proportion of fields are `NA`.


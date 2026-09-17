# Human report data dictionary

This document describes the sample-level CSV produced by the human FASTQ
workflow (`modules/local/reporthuman/main.nf`). It documents columns as
**families** because segment-specific columns use the same calculation and only
differ by segment or protein name.

## Pattern notation and dimensions

The dictionary uses two reusable tokens. They are deliberately defined once so
new segment columns can be added by extending a domain list rather than copying
a row:

| Token         | Allowed values                                                 | Examples                                                          |
| ------------- | -------------------------------------------------------------- | ----------------------------------------------------------------- |
| `<SEGMENT>`   | `HA`, `NA`, `M`/`MP`, `NP`, `NS`, `PA`, `PB1`, `PB2`           | `Coverage-<SEGMENT>`, `DEPTH_<SEGMENT>`, reassortment `<SEGMENT>` |
| `<PROTEIN>`   | `HA1`, `HA2`, `M1`, `M2`, `NA`, `NP`, `NS`, `PA`, `PB1`, `PB2` | `aaDeletions <PROTEIN>`, `Nextclade QC <PROTEIN>`                 |
| `<REFERENCE>` | `human`, `human_vaccine`, `inhibition_human`                   | `<PROTEIN> Differences <REFERENCE>`                               |

`M` and `MP` are aliases for the matrix (MP) genome segment; the spelling in a
specific output is retained for compatibility. `NS1` and `SigPep` are special
protein/reference names used only by the human difference outputs. A literal
`GISAID_<SEGMENT>_*` family, if introduced by a report variant, follows the same
rule: the token is the segment and the suffix defines the GISAID attribute; it
must not be interpreted as a new algorithm without a corresponding source row.

## Compact column inventory

| Pattern or field                                                                                         | Meaning                                                                                                    | Produced from                                                         | Main limitation                                                                                                                   |
| -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `Sample`; `RunID`; `Instrument ID`; `Date`; `Release Version`                                            | Sample and run provenance                                                                                  | Samplesheet, workflow parameters, and report execution date           | Run metadata does not fully identify reference/database revisions                                                                 |
| `Subtype`; `Sekvens_Resultat`                                                                            | Raw subtype call and display-ready result                                                                  | BLASTN classification of HA and NA; final coverage gate               | A best database hit is not a phylogenetic or phenotypic classification                                                            |
| `Coverage-<SEGMENT>`                                                                                     | Percent of selected denominator containing a non-`N` consensus character                                   | IRMA amended consensus and `coverage_finder.py`                       | Segment detection is applied to sequence text rather than its header; denominator normally falls back to observed sequence length |
| `IRMA_<STAT>`; `DEPTH_<SEGMENT>`                                                                         | IRMA read-flow counts, alternate-match noise ratio, and stage-4 segment read counts                        | IRMA `READ_COUNTS.txt` and `sequence_quality.py`                      | `DEPTH_<SEGMENT>` is not mean per-base depth; absent IRMA records become zero                                                     |
| `clade`; `clade NA`; `legacy-clade`; `subclade`; `glycosylation`; `HA_glycosylation_{1,2,3}`             | Nextclade clade and motif annotations                                                                      | Nextclade dataset-specific analysis                                   | Meaning changes with dataset and revision                                                                                         |
| `Nextclade QC <PROTEIN>`; `Nextclade Mixed Sites <PROTEIN>`                                              | Segment/protein QC status and mixed-site count                                                             | Nextclade `qc.overallStatus` and `qc.mixedSites.totalMixedSites`      | HA and M values are copied into both protein columns, so downstream mixed-site sums double them                                   |
| `aaDeletions <PROTEIN>`; `aaInsertions <PROTEIN>`; `frameShifts <PROTEIN>`                               | Reference-relative amino-acid indels and coding-frame disruptions                                          | Nextclade alignment, translation, and mutation calling                | Calls depend on consensus quality, reference coordinates, and alignment                                                           |
| `<PROTEIN> Differences <REFERENCE>`; `NS1 Differences human`; `SigPep Differences human`                 | Amino-acid substitutions relative to a named reference family                                              | Nextclade translations followed by `mutation_finder.py`               | These are reference differences, not necessarily important or causal mutations                                                    |
| `M2 inhibtion mutations`; `NA inhibtion mutations`; `PA inhibtion mutations`; `DR_Res_*`; `DR_*_Mut`     | Lookup matches, drug review codes, and mutation detail                                                     | Local Excel lookup table, `table_lookup.py`, and `report.py`          | Genotypic lookup is not a phenotypic susceptibility test; historical `inhibtion` spelling is retained                             |
| `<SEGMENT>` reassortment fields; `Reassortment`; `Conclusion`; `Origins`; `Subtypes`; `ReferenceStrains` | Best annotated reference/identity per segment and an overall screening interpretation                      | BLASTN, accession metadata, and `detect_reassortment.py`              | Similarity-based screen, not a phylogenetic reassortment analysis                                                                 |
| `Subclade_Nomenclature_<ATTRIBUTE>`                                                                      | Seasonal HA clade/subclade rule match, evidence, closest match, and source                                 | Local rule caller and influenza-clade-nomenclature YAML definitions   | Rules are downloaded from unpinned `main` branches                                                                                |
| `Characterisation_<ATTRIBUTE>`                                                                           | Guideline-based reference-virus category, closest guideline reference, and additional amino-acid mutations | Seasonal subclade call and local NH 2025/2026 characterisation tables | Genetic classification only; incomplete subclade calls are provisional                                                            |
| `NGS_QC_Sum`; `GISAID_Comment`                                                                           | Compact review flags and submission suggestion                                                             | `report_QC_calculation.py`                                            | Hard-coded review thresholds do not replace manual QC or submission validation                                                    |

`<STAT>` expands to `initial`, `failQC`, `passQC`, `chimeric`, `nomatch`,
`match`, `altmatch`, or `noise`. `<ATTRIBUTE>` expands to the suffixes listed
in the subclade section below. The report merger can also preserve arbitrary
columns from readable staged CSVs or the samplesheet; those are outside this
stable contract and must be documented by their producing module.

## How the values are calculated

### Report assembly and missing values

`report.py` reads every CSV staged in the report process, concatenates them, and
groups by `Sample`. For each column it retains the first available value in
staging order. It creates expected-but-missing columns as `NA`, changes `!` to
`-` in sample identifiers, converts coverage values to numbers rounded to two
decimals, and rounds `IRMA_noise` to five decimals. Samples present only in the
samplesheet are appended with unavailable analysis fields.

This design allows incomplete runs to produce a report, but it can hide a failed
upstream process. CSV files that cannot be read are skipped with a warning, and
multiple conflicting results are not ranked by date or quality.

### Subtype and sequence result

`Subtype` is called with NCBI BLAST+ 2.15.0 from the IRMA amended HA and NA
consensus sequences:

1. HA and NA are searched separately against the configured local reference
   FASTA files using `blastn -task blastn -outfmt 6 -max_target_seqs 2`.
2. Hits are ordered by bit score and then percent identity; the first hit is
   retained.
3. A call is accepted when alignment length is at least 20% of the expected
   coding length: 340.2 nt for HA (`0.2 × 1701`) and 282 nt for NA
   (`0.2 × 1410`). There is no explicit minimum percent-identity threshold.
4. The second underscore-separated token of each reference identifier supplies
   the HA and NA labels; the two labels are concatenated into `Subtype`.
5. `Sekvens_Resultat` maps `H3N2` to `A/H3N2`, `H1N1` to `A/H1N1`, `VIC` or
   `VICVIC` to `B/Victoria`, and `YAM` or `YAMYAM` to `B/Yamagata`. The result is
   forced to `NA` unless both `Coverage-HA` and `Coverage-NA` are at least 30%.

The database content and identifier convention are therefore part of the
method. A short local alignment can pass the call threshold, and concatenating
independent HA and NA calls can obscure mixed infections.

### Consensus coverage

Coverage is calculated from the IRMA consensus, not from aligned read
depth:

```text
Coverage-segment = 100 × (consensus length − number of N characters)
                         / selected denominator
```

The implementation receives the segment name but uses it only for the output
column. To choose a denominator, it searches the **sequence characters** for a
literal marker such as `HA-`, `NA-`, or `PB2-`; it does not inspect the FASTA
header or the segment argument. If a marker is found, the hard-coded denominator
is HA 1800, NA 1450, PB2/PB1 2400, PA 2300, NP 1600, NS 920, or M 1100 nt. If
no marker is found—as expected for an ordinary nucleotide sequence—the
denominator is the observed consensus length. In that usual case, the value is
simply the percentage of consensus characters that are not `N`, rather than the
percentage of the expected full segment recovered.

Only `N`/`n` is treated as uncovered; every other IUPAC ambiguity code counts as
covered. Terminally missing sequence is invisible when observed length is used.
The metric can exceed 100% in a hard-coded-denominator branch because it is not
capped. Values are rounded to two decimals. These behaviors are important
limitations of the current implementation and should be corrected before the
field is interpreted as segment completeness.

### IRMA assembly statistics

IRMA v1.3.5 (`IRMA FLU-minion`) assembles the filtered FASTQ reads into amended
consensus sequences. `sequence_quality.py` extracts the `Reads` value from
specific records in IRMA `READ_COUNTS.txt`:

| Pattern                                                        | Exact meaning in this report             |
| -------------------------------------------------------------- | ---------------------------------------- |
| `IRMA_initial`                                                 | Reads at record `1-initial`              |
| `IRMA_failQC`; `IRMA_passQC`                                   | Reads at `2-failQC` and `2-passQC`       |
| `IRMA_chimeric`; `IRMA_nomatch`; `IRMA_match`; `IRMA_altmatch` | Reads at the corresponding `3-*` records |
| `DEPTH_{segment}`                                              | `Reads` from each `4-*_<segment>` record |
| `IRMA_noise`                                                   | `IRMA_altmatch / IRMA_match`             |

If an expected record is absent, most IRMA counters are set to zero. A zero can
therefore mean either a true zero or a missing record. Division by zero can
produce an unavailable or infinite noise value. Interpretation depends on the
IRMA module and configuration used for the run.

### Nextclade annotations and structural mutations

Only consensus segments that pass the pipeline's configurable sequence-quality
threshold are sent to Nextclade. The module selects a subtype- and
segment-specific Nextclade dataset, aligns each sequence to that dataset's
reference, translates annotated coding sequences, calls reference-relative
mutations, assigns clades where supported, and applies the dataset's QC rules.

The report transformer reads these native Nextclade fields:

| Report family      | Native Nextclade field                           |
| ------------------ | ------------------------------------------------ |
| clade fields       | `clade`, `legacy-clade`, and `subclade`          |
| QC                 | `qc.overallStatus`                               |
| mixed sites        | `qc.mixedSites.totalMixedSites`                  |
| structural changes | `aaDeletions`, `aaInsertions`, and `frameShifts` |
| glycosylation      | dataset-provided `glycosylation` annotation      |

HA and M annotations are split into protein-specific columns (`HA1`/`HA2` and
`M1`/`M2`). The converter copies the same segment-level QC status, mixed-site
count, frameshift, insertion, and deletion value into both protein columns.
Consequently, final QC sums the HA and M mixed-site count twice: two mixed sites
in HA or M become a sum of four and pass the `MS > 3` trigger.
`HA_glycosylation_1` through `_3` are storage chunks of at most 22
comma-separated glycosylation entries each; they are not three biological
classes. Missing calls are rendered as strings such as `No frameShifts` or
`No aaDeletions`.

Nextclade positions are 1-based reference coordinates. Its QC output is a
screening aid: a warning warrants review but does not prove that a sequence is
incorrect. Dataset identity and revision must be retained because references,
clade definitions, motifs, and QC thresholds are dataset-specific.

### Reference-difference mutation families

The `Differences` columns are not taken directly from the Nextclade mutation
columns. Nextclade first produces translated amino-acid FASTA files. For each
protein, subtype, and reference family present, `mutation_finder.py` then:

1. loads `sequence_references/<family>/<subtype>/<protein>.fasta`;
2. globally aligns the reference and sample protein with Biopython
   `PairwiseAligner`, gap-open score −10 and gap-extension score −1;
3. uses the first highest-scoring alignment;
4. records substitutions as `referenceAA + 1-based position + sampleAA` and
   separates entries with semicolons; and
5. removes insertion and deletion tokens from the report-facing `Differences`
   value, writing `No mutations found` when no substitution remains.

The families are `human` for the standard human reference,
`human_vaccine` for HA/NA vaccine references, and `inhibition_human` for M2,
NA, and PA resistance-reference comparisons. Because the alignment's match and
mismatch defaults are not set explicitly, results can vary with the installed
Biopython version. A difference only means “different from this reference”; it
does not establish antigenic, clinical, or resistance significance.

### Resistance lookup and derived review codes

`table_lookup.py` filters the configured Excel database by subtype and protein
(`M2`, `NA1`, or `PA`). It compares each sample substitution with database
entries using only the numeric position plus alternate amino acid. For example,
`R143G` and `K143G` both match suffix `143G`; the reference amino acid is
ignored. Matches become `{M2,NA,PA} inhibtion mutations`; otherwise the value is
`No matching mutations found`.

`report.py` then applies a string-based summary:

| Source | Drug fields                                                                        | No database match | At least one match |
| ------ | ---------------------------------------------------------------------------------- | ----------------- | ------------------ |
| M2     | `DR_Res_Adamantine`                                                                | `AANI`            | `Review`           |
| NA     | `DR_Res_Oseltamivir`, `DR_Res_Zanamivir`, `DR_Res_Peramivir`, `DR_Res_Laninamivir` | `AANI`            | `Review`           |
| PA     | `DR_Res_Baloxavir`                                                                 | `AANS`            | `Review`           |

If the lookup source is absent, the classification is `NA`. `DR_M2_Mut`,
`DR_NA_Mut`, and `DR_PA_Mut` contain the matched mutation text when review is
required, `No Mutations` when no lookup match was found, and `NA` when the
source is absent.

These are workflow review codes, not susceptible/resistant phenotypes. The NA
classification applies the same mutation list to four drugs and does not model
mutation-specific evidence strength, mixtures, background effects, or assay
results. Reproducibility requires archiving the exact Excel lookup database.

### Reassortment screen

Each segment consensus is searched with BLASTN against the configured local
reassortment database. The preferred subject-header forms are
`ORIGIN|SUBTYPE|STRAIN|ACCESSION` and
`ORIGIN|SUBTYPE|STRAIN|SEGMENT|ACCESSION`. Legacy headers are supported; the
tracked `bin/reassortment_reference_metadata.csv` file fills origin and subtype
from accession for the current database. For every segment, the row with the
highest percent identity is selected, with bit score and alignment length used
as tie-breakers:

- `{segment}` contains `ORIGIN|SUBTYPE|STRAIN(percent_identity%)` at or above
  80%, `TooLow(identity%):ORIGIN|SUBTYPE|STRAIN` below 80%, or `Missing`.
- `Reassortment` remains a compact compatibility field: `No` for one complete
  reference profile, `Yes` for multiple complete profiles, and `Unknown` for
  missing, low-identity, or unannotated calls.
- `Conclusion` distinguishes a consistent human profile, possible
  within-subtype reassortment, human subtype discordance, mixed origins,
  non-human-only profiles, and inconclusive incomplete results.
- `Origins`, `Subtypes`, and `ReferenceStrains` list the distinct accepted
  values in compact semicolon-separated form.

This is a database-similarity screen. It does not infer phylogenetic trees,
ancestral segment exchange, or statistical support; results depend strongly on
database composition and annotation quality. Origin and subtype are reference
metadata, not inferred from the sample sequence itself.

### Seasonal HA subclade nomenclature

`Subclade_Nomenclature_*` is a custom rule-based call for supported seasonal HA
profiles: A/H3N2, A/H1N1pdm, and B/Victoria. The workflow downloads
machine-readable clade/subclade YAML definitions from the
`influenza-clade-nomenclature` repositories and obtains a profile-specific NCBI
reference sequence. The caller globally aligns HA nucleotide sequence to the
reference (match +2, mismatch −1, gap −5), translates profile features, and
compares observed nucleotide and amino-acid states with hierarchical defining
mutation rules.

| Suffix after `Subclade_Nomenclature_`                                                         | Meaning                                                                          |
| --------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `Profile`; `Source`                                                                           | Selected virus profile and rule repository                                       |
| `Clade`; `Clade_Long`; `Subclade`; `Lineage_Path`                                             | Assigned labels and parent path                                                  |
| `Key_Mutations`; `Lineage_Additive_Mutations`; `Lineage_Key_Mutations`; `Clade_Key_Mutations` | Observed rule-defining mutations at different hierarchy levels                   |
| `Closest_Subclade`; `Closest_Subclade_Missing_Mutations`; `Unique_Mutations`                  | Best incomplete match, missing defining states, and observed non-lineage changes |
| `Subclade_Match_Fraction`; `Clade_Match_Fraction`                                             | Matched rules divided by evaluated rules, rounded to three decimals              |

Candidates are ranked by matched-rule count, match fraction, hierarchy depth,
rule-set size, and name. Although the internal caller marks a non-exact result
as unassigned, the report currently writes the best candidate into both
`Subclade` and `Closest_Subclade`; missing defining mutations and the match
fraction must therefore be checked before treating the label as an exact call.
These genetic labels facilitate surveillance and do not necessarily
represent distinct phenotypes. Because the workflow downloads the current
`main` branch, rule revisions can change results unless the downloaded rules
are archived with the run.

### Reference-virus characterisation

`Characterisation_*` interprets the seasonal HA subclade result using the local
NH 2025/2026 H1, H3, B/Victoria, and B/Yamagata characterisation tables. Only a
guideline row with `reporting_category=yes` may become the primary
`Characterisation_Reference_Virus`. A sample matching a non-reporting descendant
inherits its nearest reporting ancestor and reports the mutations on the path
from that ancestor as `Characterisation_Guideline_Extra_Mutations`.

| Field                                                                                         | Meaning                                                                                                                                          |
| --------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `Characterisation_Profile`                                                                    | Selected H1, H3, B/Victoria, or B/Yamagata guideline profile                                                                                     |
| `Characterisation_Status`                                                                     | Exact reporting category, derived category with extra mutations, incomplete-call review, unassigned, or no reporting categories in the guideline |
| `Characterisation_Reference_Virus`; `Characterisation_Reference_Role`                         | Primary reportable reference-virus category and its guideline role                                                                               |
| `Characterisation_Reporting_Category_Subclade`                                                | Subclade of the primary reporting category                                                                                                       |
| `Characterisation_Closest_Guideline_Reference`; `Characterisation_Closest_Guideline_Subclade` | Most specific matching guideline row, including non-reporting rows                                                                               |
| `Characterisation_Guideline_Extra_Mutations`                                                  | Guideline mutations between the primary category and closest descendant                                                                          |
| `Characterisation_Sample_Extra_Mutations`                                                     | Additional amino-acid changes from `Subclade_Nomenclature_Unique_Mutations`; nucleotide-only changes are excluded                                |
| `Characterisation_All_Extra_Mutations`                                                        | De-duplicated combination of guideline and sample-specific additional mutations                                                                  |
| `Characterisation_Result`                                                                     | Compact display text such as `A/Victoria/4897/2022-like + R45K`                                                                                  |
| `Characterisation_Subclade_Match_Fraction`; `Characterisation_Missing_Subclade_Mutations`     | Evidence copied from the seasonal subclade call                                                                                                  |
| `Characterisation_Guideline_Source`                                                           | Characterisation CSV filename used for the result                                                                                                |

When defining mutations are missing or the subclade match fraction is below
one, the result is prefixed with `Provisional:` and must be reviewed. The
current B/Yamagata table contains no reporting-category rows, so the analysis
does not invent a primary reference-virus category; it can report the closest
non-reporting reference only when a Yamagata clade is available. These calls
describe genetic similarity to surveillance categories. They do not establish
antigenic phenotype, vaccine effectiveness, or clinical significance.

### Final NGS QC summary

`report_QC_calculation.py` evaluates each of the eight genome segments in the
order HA, NA, MP, NP, NS, PA, PB1, PB2:

- `LC` (low coverage): coverage is absent, non-numeric, or below 80%;
- `FS` (frame shift): any protein field for that segment has a non-missing value
  other than exactly `No frameShifts`;
- `MS` (mixed sites): the sum of the segment's protein-specific Nextclade mixed
  site counts is greater than 3.

Issues are written as `segment:codes` and joined with `|`, for example
`HA:MS|PB1:FS,LC`. `GISAID_Comment` is `Review` when `NGS_QC_Sum` is non-empty
and otherwise empty. Missing coverage is intentionally treated as a QC problem,
but missing frame-shift and mixed-site data are not; an empty summary therefore
does not prove that all QC inputs were present.

## Provenance required to reproduce a value

Retain these artifacts with the report:

- input samplesheet and raw/staged upstream CSVs;
- `params.json`, especially all database/reference paths and quality thresholds;
- `software_versions.yml` and container image digests;
- the exact subtype, resistance, reassortment, and sequence-reference databases;
- downloaded Nextclade datasets; and
- downloaded influenza-clade-nomenclature YAML rules and reference FASTA files.

`Release Version` alone is insufficient because several external datasets are
retrieved or supplied independently of the pipeline release.

## Software and reference sources

- Local implementations: [`report.py`](../bin/report.py),
  [`report_QC_calculation.py`](../bin/report_QC_calculation.py),
  [`coverage_finder.py`](../bin/coverage_finder.py),
  [`sequence_quality.py`](../bin/sequence_quality.py),
  [`mutation_finder.py`](../bin/mutation_finder.py),
  [`table_lookup.py`](../bin/table_lookup.py),
  [`detect_reassortment.py`](../bin/detect_reassortment.py), and
  [`subclade_nomenclature.py`](../bin/subclade_nomenclature.py).
- [IRMA documentation](https://wonder.cdc.gov/amd/flu/irma/irma.html),
  [IRMA output guide](https://wonder.cdc.gov/amd/flu/irma/output.html), and the
  [IRMA method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5011931/).
- [Nextclade tabular output definitions](https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html),
  [quality-control method](https://docs.nextstrain.org/projects/nextclade/en/stable/user/algorithm/06-quality-control.html),
  and [dataset documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).
- [NCBI BLAST+ command-line manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
  and [tabular output documentation](https://www.ncbi.nlm.nih.gov/books/NBK569862/).
- [Biopython pairwise-alignment documentation](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html).
- Seasonal influenza nomenclature definitions for
  [A/H3N2 HA](https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA),
  [A/H1N1pdm HA](https://github.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA),
  and [B/Victoria HA](https://github.com/influenza-clade-nomenclature/seasonal_B-Vic_HA),
  plus the [nomenclature method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12904685/).

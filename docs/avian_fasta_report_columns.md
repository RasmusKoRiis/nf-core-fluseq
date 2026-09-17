# Avian FASTA report column dictionary

## Scope

This document describes the sample-level CSV produced by the `avian-fasta`
workflow in `reporthuman/<runid>.csv`. It is based on the current workflow,
report merger, and analysis scripts rather than on the apparent meaning of a
column name. The example report supplied for this documentation contains 141
columns; every column in that example is covered below.

The report is a screening and annotation product. Mutation, adaptation,
antiviral, genotype, and reassortment fields require expert review and are not
phenotypic test results by themselves.

Several historical spellings are part of the CSV interface and are therefore
retained, including `mamailian` and `inhibtion`.

## How the report is assembled

The workflow splits the input multi-FASTA into individual records, identifies
the sample and segment, and assigns a deterministic internal sample identifier.
The analysis branches then produce many small CSV files. `reportfasta.py`:

1. reads every non-empty staged CSV it can parse;
2. normalises its sample-key column to `Sample`;
3. concatenates the rows;
4. groups by `Sample`; and
5. retains the first non-empty value that is not `NA`, `NaN`, or `None` for
   each column.

The report process then left-joins these results to `id_map.tsv`, ensuring one
row for every `SampleID`, adds run metadata, calculates final QC fields, and
writes the final CSV.

This has four important consequences:

- The output schema is partly dynamic. A column produced only by a failed or
  inapplicable analysis branch may be absent or contain `NA`.
- `NA` means unavailable or not produced. It does not automatically
  mean a negative biological result.
- Because process errors are ignored globally, a final report can be produced
  even when one or more upstream analyses failed. The Nextflow log and
  `pipeline_info/` must be reviewed with the CSV.

## Value and mutation notation

| Value or token | Meaning |
| --- | --- |
| `NA` | Result unavailable, source analysis absent, or value not applicable. |
| blank cell | Usually no final QC flag. In some source columns it can also mean that an upstream table was empty. |
| `Review` | A workflow rule found something requiring review. |
| `No matching mutations found` | The calculated substitutions did not match the configured lookup table. |
| `No Mutations` | Report-level summary used when the lookup completed without a matching listed mutation. |
| `A123T` | Reference amino acid `A`, 1-based reference position `123`, sample amino acid `T`. |
| `A123A` | Sample and reference agree at position 123. These matching tokens occur in full amino-acid lists. |
| `A123X` | The translated sample has an unknown amino acid at that aligned position. |
| `A123-` | Deletion/gap in the sample relative to the amino-acid reference. |
| `ins123...` / `del123...` | Internal insertion/deletion notation. These tokens are removed from report-facing `Differences` columns. |

Mutation coordinates are relative to the named reference sequence, not a
universal influenza coordinate system.

## Sample and run identity columns

| Column | What it contains | How it is calculated |
| --- | --- | --- |
| `SampleID` | Internal alphanumeric identifier used to join all results for one sample. | The selected original-name core is uppercased, all non-alphanumeric characters are removed, and the first eight hexadecimal characters of the SHA-1 digest of that normalised name are appended. Example: `AREDFOXSCOTLAND023888202326B3E7E5`. |
| `OriginalName` | Human-readable sample or strain name recovered from the original FASTA header. | The workflow compares the text before and after the first `|` and normally selects the side containing more `/` characters; ties select the left side. The wrapper's earlier header conversion can therefore affect this field. |
| `RunID` | Operational run identifier. | Copied from `--runid`; the avian wrapper supplies its `--run` value. |
| `Instrument ID` | Instrument label recorded for the run. | Copied from `--seq_instrument`. It is metadata and is not inferred from the FASTA. |
| `Date` | Date on which the final report process executed. | Generated using the container's current date. The CSV uses ISO `YYYY-MM-DD`; spreadsheet software may display it as local `DD.MM.YYYY`. It is not necessarily the collection or sequencing date. |
| `Release Version` | Operational method/release label. | Copied from `--release_version`. The wrapper currently supplies a literal release label, so this field is not guaranteed to be the Git commit or branch actually executed. |
| `Sample` | Intermediate sample key found in upstream CSV files. | Normally the internal UID. The report process removes it when it mirrors `SampleID` for all rows, but it can remain in legacy or mixed-key reports. Use `SampleID` as the authoritative join key. |

## Subtype and sequence-result columns

| Column | What it contains | How it is calculated |
| --- | --- | --- |
| `Subtype` | Combined HA and NA subtype call, for example `H5N1`. | HA and NA nucleotide sequences are searched separately with BLASTN against the configured HA and NA databases. Hits are ordered by bit score and then percent identity. The top HA hit is accepted when alignment length is at least `0.2 × 1701 = 340.2` nt; the top NA hit is accepted at `0.2 × 1410 = 282` nt. The second underscore-delimited token of each accepted database identifier is used, and the HA and NA tokens are concatenated. There is no explicit minimum percent-identity threshold. |
| `Sekvens_Resultat` | Display-ready subtype result. | Derived from `Subtype`, but only reported when both `Coverage-HA` and `Coverage-NA` are numeric and at least 30%. `H3N2` becomes `A/H3N2`, `H1N1` becomes `A/H1N1`, `VIC`/`VICVIC` becomes `B/Victoria`, and `YAM`/`YAMYAM` becomes `B/Yamagata`. Other values such as `H5N1` are retained unchanged. Otherwise the value is `NA`. |

`Subtype` and `Subtypes` are different fields. `Subtype` is the direct HA/NA
classification above; plural `Subtypes` is a reassortment-reference summary
described later.

## Consensus coverage columns

The eight example columns are:

- `Coverage-NP`
- `Coverage-HA`
- `Coverage-NA`
- `Coverage-PB2`
- `Coverage-PB1`
- `Coverage-M`
- `Coverage-PA`
- `Coverage-NS`

Each is calculated independently from that segment's consensus sequence:

```text
coverage (%) = 100 × (sequence length − number of N/n characters)
                   / selected denominator
```


Only `N` is counted as uncovered. Other IUPAC ambiguity codes count as covered,
and terminal sequence that is absent from the FASTA cannot be detected when the
observed length is used as denominator. Values are rounded to two decimals in
the report.

The coverage-filter module passes a segment downstream only when its value is
strictly greater than `--seq_quality_threshold` (default 80). The final QC rule
flags `LC` only below 80.

## Antiviral lookup and drug-review columns

Translated M2, NA, and PA proteins are globally aligned to references in the
`inhibition` reference family. The resulting substitutions are compared with
the configured inhibition Excel table. The lookup filters by segment and
subtype and matches using the numeric position plus alternate amino acid. The
leading reference amino acid is ignored; for example, `R143G` and `K143G` both
match a table entry ending in `143G`.

### Lookup result columns

| Column | What it contains |
| --- | --- |
| `M2 inhibtion mutations` | Sample M2 substitutions matching the configured inhibition table, separated by semicolons; otherwise `No matching mutations found`. |
| `NA inhibtion mutations` | Sample NA substitutions matching the inhibition table. Upstream `NA1 inhibtion mutations` is normalised to this column name. |
| `PA inhibtion mutations` | Sample PA substitutions matching the inhibition table. |

### Derived review-code columns

| Column | Source and rule |
| --- | --- |
| `DR_Res_Adamantine` | From `M2 inhibtion mutations`: `NA` if source unavailable, `AANI` for `No matching mutations found`, otherwise `Review`. |
| `DR_Res_Oseltamivir` | From `NA inhibtion mutations` using the same rule, with no-match code `AANI`. |
| `DR_Res_Zanamivir` | Same source and calculation as Oseltamivir. |
| `DR_Res_Peramivir` | Same source and calculation as Oseltamivir. |
| `DR_Res_Laninamivir` | Same source and calculation as Oseltamivir. |
| `DR_Res_Baloxavir` | From `PA inhibtion mutations`: `NA` if unavailable, `AANS` for no table match, otherwise `Review`. |
| `DR_M2_Mut` | Matching M2 mutation text when Adamantine is `Review`; `No Mutations` after a completed no-match lookup; otherwise `NA`. |
| `DR_NA_Mut` | Matching NA mutation text when any of the four NA drug fields is `Review`; `No Mutations` after a completed no-match lookup; otherwise `NA`. |
| `DR_PA_Mut` | Matching PA mutation text when Baloxavir is `Review`; `No Mutations` after a completed no-match lookup; otherwise `NA`. |

`AANI`, `AANS`, and `Review` are literal workflow reporting codes. The pipeline
does not calculate a drug-specific MIC, inhibition value, susceptible/resistant
phenotype, evidence grade, or clinical recommendation. In particular, all four
NA-inhibitor columns are derived from the same undifferentiated NA mutation
list.

## FluMut marker and effect columns

The `Flumut_*` and `Effect_*` columns come from FluMut's marker table. At run
time the pipeline attempts to update FluMutDB and falls back to the database in
the container if the update is unavailable. FluMut evaluates the combined
segment FASTA and returns marker mutations and literature-derived effect text.

| Protein | Mutation column | Effect column |
| --- | --- | --- |
| HA1 | `Flumut_HA1` | `Effect_HA1` |
| HA2 | `Flumut_HA2` | `Effect_HA2` |
| M1 | `Flumut_M1` | `Effect_M1` |
| M2 | `Flumut_M2` | `Effect_M2` |
| NA | `Flumut_NA` | `Effect_NA` |
| NP | `Flumut_NP` | `Effect_NP` |
| NS1 | `Flumut_NS1` | `Effect_NS1` |
| PA | `Flumut_PA` | `Effect_PA` |
| PB1 | `Flumut_PB1` | `Effect_PB1` |
| PB2 | `Flumut_PB2` | `Effect_PB2` |

For each FluMut marker row, the converter removes the marker prefix from
`Mutations in your sample`, assigns it to the matching protein column, removes
duplicate mutation strings within that protein, and joins values with `;`.
Effect strings are joined in the corresponding order. Repeated effect text is
possible when different mutations have the same literature annotation.


## GenIn2 genotype column

| Column | What it contains | How it is calculated |
| --- | --- | --- |
| `Genotype_Genin2` | Avian-influenza genotype label emitted by GenIn2, for example `EA-2021-AB`. | The combined, segment-labelled sample FASTA is passed directly to `genin2`. The workflow does not reinterpret the returned genotype; it renames GenIn2's `Genotype` field and keeps the first row per sample. The current invocation uses GenIn2's default minimum per-segment sequence coverage of 0.7. Genotype definitions and reference data are supplied by the pinned GenIn2 container. |

## Reassortment similarity-screen columns

The segment-named fields are **not sequence fields** and are not the direct
HA/NA subtype call. They report each segment's best reassortment-database match:

- `HA`
- `NA`
- `MP` (matrix segment)
- `NP`
- `NS`
- `PA`
- `PB1`
- `PB2`

All combined segment consensuses are searched with BLASTN against the configured
reassortment database. For each segment, the hit with the highest percent
identity is selected; bit score and alignment length break ties. The result is:

| Value | Meaning |
| --- | --- |
| `ORIGIN|SUBTYPE|STRAIN(percentage%)` | Best hit has at least 80% identity. An unavailable annotation can appear as `UNKNOWN` within the formatted value. |
| `TooLow(percentage%):ORIGIN|SUBTYPE|STRAIN` | Best hit is below 80% identity. |
| `Missing` | No segment hit was available. |

There is no explicit minimum query-coverage rule in this screen, so a short
high-identity alignment can be selected. The percentage is BLAST alignment
identity, not whole-segment completeness.

### Overall reassortment fields

| Column | What it contains | How it is calculated |
| --- | --- | --- |
| `Reassortment` | Compact compatibility result: `No`, `Yes`, or `Unknown`. | `Unknown` when any expected segment is missing, below 80%, or lacks reference metadata. Otherwise `No` when all accepted segments match one `(origin, subtype, strain)` profile, and `Yes` when more than one profile is represented. |
| `Conclusion` | Detailed rule-based interpretation of the segment matches. | Reports incomplete evidence, mixed origins, subtype discordance, non-human origin, or multiple reference strains. Examples include `INCONCLUSIVE`, `ALERT`, `REVIEW`, `FLAG`, and `CONSISTENT` messages. |
| `Origins` | Distinct origins among accepted segment hits. | Unique, sorted database annotations joined with `;`, excluding missing/low hits. |
| `Subtypes` | Distinct reference subtypes among accepted segment hits. | Unique, sorted database annotations joined with `;`. This can differ from singular `Subtype`. |
| `ReferenceStrains` | Distinct reference strain names among accepted segment hits. | Unique, sorted accepted reference annotations joined with `;`. |

This is a database-similarity screen, not a phylogenetic reassortment analysis.
It does not infer an ancestral reassortment event or statistical support.

## Reference-difference columns

The mutation workflow translates coverage-passing segment sequences and
globally aligns each sample protein to configured amino-acid references using
Biopython `PairwiseAligner`, gap-open score -10 and gap-extension score -1. It
uses the first highest-scoring alignment.

Three reference families can contribute to the example report:

| Name in column | Reference family | Purpose |
| --- | --- | --- |
| `mamailian` | `sequence_references/mamailian/<subtype>/` | Differences from configured mammalian-characterisation references. Historical misspelling retained. |
| `human_vaccine` | `sequence_references/human_vaccine/<subtype>/` | Differences from configured vaccine reference proteins. In the current workflow this branch is run for HA/NA products of H5N1 and H5N5. |
| `inhibition` | `sequence_references/inhibition/<subtype>/` | Full reference-relative protein lists used before the resistance-table lookup, principally for M2, NA, and PA. |

### Difference columns

Each `Differences` value contains only reference-relative amino-acid
substitutions separated by `;`. Alignment insertion/deletion tokens are removed
from this field. `No mutations found` means no remaining substitution relative
to that particular reference; it does not mean identity to every influenza
reference.

| Columns in the example | Interpretation |
| --- | --- |
| `HA Differences human_vaccine` | HA-product substitutions relative to the configured H5 vaccine reference. |
| `HA2 Differences human_vaccine` | HA2 substitutions relative to the configured vaccine HA2 reference. |
| `NA1 Differences human_vaccine` | NA1 substitutions relative to the configured vaccine NA reference. |
| `HA Differences mamailian`; `HA2 Differences mamailian`; `M Differences mamailian`; `NA1 Differences mamailian`; `NP Differences mamailian`; `NS Differences mamailian`; `PA Differences mamailian`; `PB1 Differences mamailian`; `PB2 Differences mamailian`; `SIG Differences mamailian` | Substitutions relative to the corresponding configured mammalian-characterisation protein reference. `SIG` is the HA signal peptide. |

### `_1`, `_2`, and `_3` mammalian-difference columns

The following columns are storage/display chunks of their unsuffixed source:

- `HA Differences mamailian_1`, `_2`, `_3`
- `HA2 Differences mamailian_1`, `_2`, `_3`
- `M Differences mamailian_1`, `_2`, `_3`
- `NA1 Differences mamailian_1`, `_2`, `_3`
- `NP Differences mamailian_1`, `_2`, `_3`
- `NS Differences mamailian_1`, `_2`, `_3`
- `PA Differences mamailian_1`, `_2`, `_3`
- `PB1 Differences mamailian_1`, `_2`, `_3`
- `PB2 Differences mamailian_1`, `_2`, `_3`
- `SIG Differences mamailian_1`, `_2`, `_3`

The report process splits semicolon-delimited tokens into up to three chunks,
aiming for at most 140 characters per chunk. If the list still exceeds three
chunks, all remaining tokens are appended to `_3`, so `_3` can exceed 140
characters. These columns are not three references, replicates, confidence
levels, or mutation classes. The complete value remains in the unsuffixed
column.

### Full amino-acid list columns

The example contains:

- `HA human_vaccine full amino acid list`
- `HA mamailian full amino acid list`
- `HA2 human_vaccine full amino acid list`
- `HA2 mamailian full amino acid list`
- `M mamailian full amino acid list`
- `NA1 human_vaccine full amino acid list`
- `NA1 inhibition full amino acid list`
- `NA1 mamailian full amino acid list`
- `NP mamailian full amino acid list`
- `NS mamailian full amino acid list`
- `PA inhibition full amino acid list`
- `PA mamailian full amino acid list`
- `PB1 mamailian full amino acid list`
- `PB2 mamailian full amino acid list`
- `SIG mamailian full amino acid list`

These fields list every aligned reference position, including matches such as
`A123A`, substitutions such as `A123T`, unknown residues such as `A123X`, and
sample gaps such as `A123-`. Insertions relative to the reference are not added
to `All_Positions`, so they are not represented in these full-list fields.

Large runs of `X` indicate unknown translated residues, commonly caused by
ambiguous nucleotides, incomplete coding sequence, or an unsuitable reading
frame/alignment. Large terminal runs of `-` indicate that the sample alignment
does not cover the corresponding reference region. These patterns require
review; they should not be read as hundreds of established biological
substitutions.

### Reference-name columns

| Column | What it contains | Limitation |
| --- | --- | --- |
| `Mutation reference` | Description of a reference FASTA record from the `mamailian` comparison, with its terminal segment suffix removed. | Many proteins can contribute a field with this same column name. The final merger keeps the first non-empty value, so this single field does not reliably document the reference used for every `mamailian` difference column. |
| `Vaccine mutation reference` | Description of a reference FASTA record from the `human_vaccine` comparison. | HA, HA2, and NA results share one final column name; the merger keeps the first non-empty value. Consult the controlled reference directory and reference manifest for complete provenance. |

## Nextclade HA columns

For H5 samples, the `NEXTCLADE` branch analyses HA with the community
`iav-h5/ha/all-clades` dataset. Other H5 segments are translated through a
separate local-dataset branch for mutation comparisons, but they do not
normally add native Nextclade summary columns to this report. This is why the
example contains Nextclade columns for HA1 and HA2 only.

Only segments passing the coverage-filter threshold reach these branches.

| Column | Native source and calculation |
| --- | --- |
| `Nextclade Mixed Sites HA1` | Nextclade `qc.mixedSites.totalMixedSites` for the HA segment. |
| `Nextclade Mixed Sites HA2` | The same HA segment-level value copied into the HA2 column. It is not independently counted for HA2. |
| `Nextclade QC HA1` | Nextclade `qc.overallStatus` for HA, copied to HA1. Values such as `good` are defined by the selected dataset's QC rules. |
| `Nextclade QC HA2` | The same HA segment-level QC status copied to HA2. |
| `aaDeletions HA1` | Native HA `aaDeletions`, copied to the HA1 report field. |
| `aaDeletions HA2` | The same native segment-level deletion field copied to HA2. |
| `aaInsertions HA1` | Native HA `aaInsertions`, copied to HA1. |
| `aaInsertions HA2` | The same native segment-level insertion field copied to HA2. |
| `frameShifts HA1` | Native HA `frameShifts`, copied to HA1; missing values become `No frameShifts`. |
| `frameShifts HA2` | The same native segment-level frameshift field copied to HA2. |
| `legacy-clade` | Nextclade dataset's legacy HA clade assignment. |
| `subclade` | Nextclade dataset's HA subclade assignment. |
| `glycosylation` | Dataset-provided HA glycosylation annotation. Commas are converted to semicolons during report conversion. |
| `HA_glycosylation_1` | First storage chunk of up to 22 comma-delimited glycosylation entries. |
| `HA_glycosylation_2` | Second storage chunk, populated only when more than 22 entries exist. |
| `HA_glycosylation_3` | Remaining entries after the first 44. |

The HA1/HA2 duplication is important for final QC: the two mixed-site columns
are summed even though they hold the same segment-level count. Consequently,
two native HA mixed sites appear as `2 + 2 = 4` and trigger the current
`MS > 3` rule.

Nextclade results depend on the exact dataset and reference. The process first
attempts to download the current official/community dataset and falls back to
the configured local dataset root when downloading is unavailable.

## Final QC columns

### `NGS_QC_Sum`

Compact segment-wise QC flags in this fixed order:

```text
HA | NA | MP | NP | NS | PA | PB1 | PB2
```

For each segment, the script can add:

| Code | Trigger |
| --- | --- |
| `LC` | Coverage is absent, non-numeric, or below 80%. Matrix coverage is read from `Coverage-M` or `Coverage-MP`. |
| `FS` | Any available segment/protein frameshift field is not exactly `No frameShifts` after case normalisation. |
| `MS` | Sum of available Nextclade mixed-site fields for that segment is greater than 3. |

Multiple issue codes are alphabetically sorted within a segment, and flagged
segments are joined with `|`; for example:

```text
HA:MS|PB1:FS,LC|NP:LC
```

The QC summary does not inspect `Conclusion`, `Reassortment`, FluMut effects,
antiviral review fields, GenIn2 genotype, `X` residues, the singular subtype,
or Nextclade's overall status text. Therefore a blank QC summary does not mean
that every other report field is unremarkable.

### `GISAID_Comment`

Set to `Review` when `NGS_QC_Sum` is non-empty. It is written as a blank cell
when there are no `LC`, `FS`, or `MS` flags. This is a submission-review hint,
not a complete GISAID validation.

## Reading the supplied mock row

The example row illustrates several independent analyses:

- `Subtype=H5N1` is the direct HA/NA classification. Because `Coverage-HA`
  and `Coverage-NA` are both above 30%, `Sekvens_Resultat` is also `H5N1`.
- All eight coverage values are above 80%, so none triggers `LC`.
- Both HA mixed-site fields are zero and both frameshift fields say
  `No frameShifts`, so the shown fields do not trigger `MS` or `FS`.
  Consequently `NGS_QC_Sum` and `GISAID_Comment` can be blank.
- `HA=Missing` belongs to the independent reassortment-database screen. It
  does **not** mean the sample lacks an HA consensus: the row has 91.44% HA
  non-`N` coverage and an HA-derived H5 subtype call.
- `Conclusion=INCONCLUSIVE - missing segments: HA; observed subtype
  discordance: H1N1,H3N2` summarises the reassortment reference matches, not
  the direct H5N1 subtype call. The other segments' best database matches were
  annotated as human H1N1 and H3N2 references.
- `Reassortment=Unknown` follows automatically because an expected segment was
  missing from that screen.
- `DR_Res_Oseltamivir`, `DR_Res_Zanamivir`, `DR_Res_Peramivir`, and
  `DR_Res_Laninamivir` are all `Review` because they are all derived from the
  same matched NA mutation list containing `S247N`.
- The long HA difference/full-list fields contain many `X` and `-` tokens.
  Those symbols describe unknown amino acids and alignment gaps relative to
  the configured references and warrant alignment/reference review; they are
  not automatically validated biological mutations.
- A blank `NGS_QC_Sum` does not clear the reassortment `INCONCLUSIVE` result or
  the antiviral `Review` fields because those are outside the final QC formula.

## Primary implementation files

- Workflow wiring: [`workflows/avian-fasta.nf`](../workflows/avian-fasta.nf)
- Final report process: [`modules/local/reportavianfasta/main.nf`](../modules/local/reportavianfasta/main.nf)
- CSV merger and DR summaries: [`bin/reportfasta.py`](../bin/reportfasta.py)
- Final QC: [`bin/report_QC_calculation.py`](../bin/report_QC_calculation.py)
- Coverage: [`bin/coverage_finder.py`](../bin/coverage_finder.py)
- Reference differences: [`bin/mutation_finder.py`](../bin/mutation_finder.py)
- Inhibition lookup: [`bin/table_lookup.py`](../bin/table_lookup.py)
- Nextclade conversion: [`bin/nextclade_converter.py`](../bin/nextclade_converter.py)
- FluMut conversion: [`bin/flumut_conversion.py`](../bin/flumut_conversion.py)
- Reassortment screen: [`bin/detect_reassortment.py`](../bin/detect_reassortment.py)
- GenIn2 reduction: [`bin/slim_genin2_report.py`](../bin/slim_genin2_report.py)

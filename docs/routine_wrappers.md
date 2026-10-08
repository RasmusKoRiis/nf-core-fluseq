# Routine wrapper contract

Wrapper copies are maintained outside this repository in
`/home/rasmuskopperud.riis/Coding/flu-wrappers` and
`/home/rasmuskopperud.riis/Coding/ngs_scripts/fluseq`. These are separate files and
currently differ in reference validation and download behaviour. The
`ngs_scripts/fluseq/README.md` run examples use
`/home/ngs/ngs_scripts/fluseq/` on the server. Check the path of the wrapper you
actually run when applying updates.

The launchers are:

- `fluseq_wrapper.sh` — human FASTQ;
- `fluseq_fasta_wrapper.sh` — human FASTA;
- `avianseq_wrapper.sh` — avian FASTQ;
- `avianseq_fasta_wrapper.sh` — avian FASTA.

The infrastructure changes preserve their historical parameter names and the `<outdir>/reporthuman/*.csv` upload location.

## Selecting pipeline code

The human wrappers and avian FASTQ wrapper accept `-b <branch-or-tag>`. The avian FASTA wrapper uses `-g <branch-or-tag>` because `-b` already denotes its source option. Use `infrastructure` for acceptance testing. After merge, use a release tag for routine production.

## Seqera/Tower credential

Credentials are no longer stored in the wrapper source. The wrappers use an inherited `TOWER_ACCESS_TOKEN`, or source a private file from:

```text
$HOME/.config/fluseq/tower.env
```

The file should contain an exported credential and be readable only by the service account:

```bash
export TOWER_ACCESS_TOKEN='<token>'
```

Set `FLUSEQ_TOWER_ENV` to use a different file. Runs continue with a warning when the credential is absent.

## File-safety behavior

- A wrapper refuses to start if `$HOME/<run-id>` already exists, preventing stale and new pipeline outputs from mixing.
- If `$HOME/out_fluseq/<run-id>` already exists after a successful analysis, it is renamed with a `.previous.<timestamp>` suffix before the new result is moved into place.
- Required databases, reference directories, sample sheets, and analysis inputs are checked before Nextflow starts.
- Temporary downloaded input is kept separate from published output.

The wrappers still obtain the controlled seasonal data from the existing SMB locations and invoke `-profile docker,server`.
They also default `NXF_SYNTAX_PARSER` to `v1`, which keeps the current DSL2 pipeline compatible with Nextflow 26.04 and newer while allowing an explicit environment override.

## EPI validation in the ngs_scripts human FASTQ wrapper

The `ngs_scripts/fluseq/fluseq_wrapper.sh` version requires every reference
header to contain a strain name, an EPI identifier, and a terminal segment suffix:

```text
>A/Example/1/2025|EPI_ISL_123456_HA1
```

The corresponding semicolon-separated reference table must contain the named
columns `Subtype;Reference;Type;GISAID_EPI`. The example header would match:

```text
H1N1;A/Example/1/2025;human;EPI_ISL_123456
```

These are illustrative values. Use the verified strain and identifier for the
actual reference. Both must match the table. The wrapper checks every FASTA
record and verifies the segment suffix against the filename, allowing established
aliases. Headers without the pipe-delimited EPI identifier fail with
`Missing EPI identifier`.

Routine runs download both `references/human/` and `references/human_vaccine/`
from the selected seasonal SMB directory into the server's
`sequence_references/` tree, then validate both types against the table. Both
types need table entries and populated FASTA directories with EPI headers.

Correct reference files in the seasonal SMB source as well as the deployed
bundle: the human reference download can overwrite edits made only to the
server's local copy. This wrapper provides an offline check that does not start
downloads or Nextflow:

```bash
bash /home/ngs/ngs_scripts/fluseq/fluseq_wrapper.sh --check-references \
  /mnt/tempdata/influensa_db/flu_seq_db/sequence_references \
  /mnt/tempdata/influensa_db/flu_seq_db/reference_table.csv
```

The separate `flu-wrappers` copies described below do not implement this EPI
check. Updates to those copies do not update `ngs_scripts/fluseq` automatically.

## Human mutation reference updates in the flu-wrappers copies

Both human wrappers in `Coding/flu-wrappers` recursively download `references/human/` and
`references/human_vaccine/` from the selected `Sesongfiler/${SEASON}` directory on
SMB. The subtype directories and FASTA filenames are preserved under
`/mnt/tempdata/influensa_db/flu_seq_db/sequence_references/human/` and
`/mnt/tempdata/influensa_db/flu_seq_db/sequence_references/human_vaccine/`.

Both reference types are validated against the downloaded `reference_table.csv`
before Nextflow starts. This semicolon-separated table must contain rows for
both `human` and `human_vaccine`, with columns in this order:
`Subtype;Reference;Type;GISAID`. Within each listed subtype directory, the first
FASTA header in every `.fasta` file must identify the same strain as the table.
Use `>strain_name_SEGMENT`, such as `>A/Victoria/2570/2019_HA1`; the wrapper removes
the final underscore suffix and normalizes slashes, spaces, and underscores
before comparing the strain names.

In these copies, the fourth column is read as GISAID metadata but is not
validated and may be empty. Their legacy header parser does not separately
recognize EPI identifiers; appending EPI metadata to a header can cause a
strain-name mismatch. This differs from the strict `ngs_scripts` FASTQ wrapper
above.

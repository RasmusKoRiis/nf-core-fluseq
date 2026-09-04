# Routine wrapper contract

The operational launchers are maintained outside this repository in `/home/rasmuskopperud.riis/Coding/flu-wrappers`:

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

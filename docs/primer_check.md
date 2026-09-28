# Primer checking

Primer checking is enabled by default for the human FASTQ and FASTA workflows
and runs PCR checks for influenza. Routine wrappers supply the database path;
direct pipeline runs must supply `--primer_check_pcr /path/to/database`.
Use `--primer_check false` to disable it.

The shared primer-checker integration also supports PCR plus NGS checks for
SARS-CoV-2/RSV.

See the [shared integration guide](https://github.com/RasmusKoRiis/primer-checker/blob/main/docs/PIPELINE_INTEGRATION.md)
for database layouts, wrapper overrides, latest-image deployment, ignored errors,
output interpretation and synthetic tests. The default container must be
published before first production use. Primer tasks use `errorStrategy 'ignore'`
and publish CSV/HTML under `primer_check/`; inspect `task_status.csv` for failures.

The local standalone harness is `tests/primer_check/main.nf`, with its own
`nextflow.config`. It can exercise the module using synthetic manifests without
running the sequencing pipeline or installing nf-core plugins. The canonical
harness and automated tests are maintained in primer-checker.

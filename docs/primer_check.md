# Primer checking

The optional primer-checker integration runs PCR checks for influenza and PCR
plus NGS checks for SARS-CoV-2/RSV. Routine wrappers enable it; direct pipeline
runs use `--primer_check true --primer_check_pcr /path/to/database`.

See the [shared integration guide](https://github.com/RasmusKoRiis/primer-checker/blob/main/docs/PIPELINE_INTEGRATION.md)
for database layouts, wrapper overrides, latest-image deployment, ignored errors,
output interpretation and synthetic tests. The default container must be
published before first production use. Primer tasks use `errorStrategy 'ignore'`
and publish CSV/HTML under `primer_check/`; inspect `task_status.csv` for failures.

The local standalone harness is `tests/primer_check/main.nf`, with its own
`nextflow.config`. It can exercise the module using synthetic manifests without
running the sequencing pipeline or installing nf-core plugins. The canonical
harness and automated tests are maintained in primer-checker.

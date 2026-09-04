# nf-core/fluseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Infrastructure - 2026-09-03

### `Added`

- Schema-driven parameter validation for all four analysis modes.
- Routine and canonical FASTQ input parsing tests.
- Reference SHA-256 manifests and controlled-reference documentation.
- CI coverage for Nextflow 24.10.2 and current stable, Python tests, nf-test, and formatting.

### `Changed`

- Reorganized reusable FASTA and reporting components into local modules and shared utilities.
- Replaced implicit large reference defaults with explicit run inputs.
- Restricted publishing to documented user-facing outputs while retaining `reporthuman/` wrapper compatibility.
- Pinned custom containers and runtime-downloaded subclade rules.

### `Fixed`

- Legacy routine parameter aliases, server work-directory configuration, executable pipeline utilities, and failing version probes.
- Duplicate lifecycle handlers, stale template tests/workflows, and runtime Nextclade latest-dataset downloads.

### `Deprecated`

- `--samplesDir`, `--seq_quality_thershold`, `--mamalian_mutation_db`, and `--inhibtion_mutation_db`; these remain supported for routine wrappers.

## v1.0dev - [date]

Initial release of nf-core/fluseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

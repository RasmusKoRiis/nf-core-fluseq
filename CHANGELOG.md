# nf-core/fluseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Infrastructure - 2026-09-11

### `Added`

- Schema-driven parameter validation for all four analysis modes.
- Routine and canonical FASTQ input parsing tests.
- Reference SHA-256 manifests and controlled-reference documentation.
- CI coverage for Nextflow 24.10.2 and current stable, Python tests, nf-test, and formatting.
- Metadata for all 35 local modules, with offline schema and channel-contract validation.

### `Changed`

- Reorganized reusable FASTA and reporting components into local modules and shared utilities.
- Replaced implicit large reference defaults with explicit run inputs.
- Restricted publishing to documented user-facing outputs while retaining `reporthuman/` wrapper compatibility.
- Pinned custom containers and runtime-downloaded subclade rules.

### `Fixed`

- Legacy routine parameter aliases, server work-directory configuration, executable pipeline utilities, and failing version probes.
- Duplicate lifecycle handlers, stale template tests/workflows, and runtime Nextclade latest-dataset downloads.
- Hidden input-check exceptions, Nextclade feature-to-reference mapping, HUMANMUTATION staging, and REPORTHUMAN CSV conversion.
- Global ignored process failures; required process errors now produce a failed workflow status after submitted tasks finish.
- An ignore rule that excluded `docs/testing.md`, and overly broad advisory lint exceptions that concealed technical migration gaps.

### `Deprecated`

- `--samplesDir`, `--seq_quality_thershold`, `--mamalian_mutation_db`, and `--inhibtion_mutation_db`; these remain supported for routine wrappers.

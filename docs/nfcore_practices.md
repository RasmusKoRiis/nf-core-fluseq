# Use of nf-core practices

This repository is maintained under the personal `RasmusKoRiis` GitHub account. It is not an official nf-core pipeline and has not undergone nf-core community review. Keeping the personal repository does not prevent improving its technical quality.

## Current status

The repository is **not yet fully aligned with current nf-core standards**. Passing selected CI checks or an advisory lint run with exclusions does not establish full compliance or biological acceptance.

Implemented checks include parameter-schema validation, synthetic input parsing and reporting tests, reference provenance, Nextflow configuration/help checks, and selected nf-tests. All 35 local processes now have `meta.yml` describing their tools, inputs, outputs, and maintainers. Processes that previously shared a module file have separate directories. An offline Python test validates every local metadata file against a vendored copy of the nf-core module schema and checks it against the process declarations.

The following work remains:

| Area                        | Remaining work                                                                                                                                                                                                                                                                    |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Software environments       | Only `cat_fastq` has `environment.yml`; the other 34 local modules need dependency review and tested environments where Conda support is intended. Most currently use containers. `seqkitfasta` also lacks an explicit runtime container and reports a hard-coded Alpine version. |
| Module execution tests      | Dedicated nf-tests currently cover `cat_fastq` and `surveillance_summary`. Other selected processes have Python-driven Nextflow regression tests, but coverage is incomplete.                                                                                                     |
| Stub runs                   | Only `cat_fastq` has a stub block. A complete pipeline stub run is not supported.                                                                                                                                                                                                 |
| Public pipeline tests       | `conf/test.config`, `conf/test_full.config`, and a pipeline-level `tests/default.nf.test` are missing. These need synthetic or redistributable fixtures and reference bundles. Private data is not a reason to omit public tests.                                                 |
| Current template and syntax | The workflow still uses legacy classes in `lib/`, nf-schema 2.4 compatibility settings, and the v1 syntax parser. Migration requires regression testing.                                                                                                                          |
| Module conventions          | Some local modules still use legacy names, lack execution guards or configurable arguments, or need more accurate software-version capture. Metadata documents the existing interfaces; it does not certify those implementations.                                                |
| Biological acceptance       | A full representative run with the intended references, plus review of mutation calls and final reports, remains unverified here.                                                                                                                                                 |

See the [nf-core module guidelines](https://nf-co.re/docs/guidelines/components/modules) for the conventions being assessed, and [testing.md](testing.md) for reproducible checks and private acceptance instructions.

## Lint scope

Run the advisory audit with nf-core tools 4.1.0:

```bash
NXF_SYNTAX_PARSER=v1 nf-core pipelines lint --dir .
```

`.nf-core.yml` lists exceptions for personal repository branding, nf-core institutional configuration, AWS automation, RO-Crate README synchronization, deliberately customized template files, and specific Python tests that construct literal Nextflow braces. The file/configuration checks are enabled with individual exceptions; technical migration findings remain visible and the advisory audit is expected to report failures until those gaps are resolved.

The audit on 2026-09-11 reported nine failures: the three missing pipeline test files listed above, four legacy files in `lib/`, the old `validation.failUnrecognisedParams` setting, and the absent `test` profile. Separately, the local regression suite passed 80 tests, the existing nf-test suites passed six tests, and configuration/help checks passed. These passing checks do not cover the outstanding items in the table.

Do not disable a whole technical check to obtain a zero-failure result. Review each finding, fix actionable defects, and document any justified exception. The legacy utility files and missing public test profiles remain actionable migration findings, not ownership exceptions.

## Module metadata

Each `modules/local/<module>/meta.yml` sits beside its `main.nf`. Inputs preserve process order and tuple structure. Outputs are keyed by the process's `emit` names; patterns describe staged or emitted files. Legacy spelling is retained where it is part of the current interface. Optional output descriptions explain when a channel may have no result.

Metadata and environment files serve different purposes: `meta.yml` documents a component, while `environment.yml` defines packages for Conda execution. An environment file must contain verified dependencies, not guessed versions copied from another module. New modules must include metadata and meaningful execution tests; changes to existing interfaces must update their metadata in the same commit.

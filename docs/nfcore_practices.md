# Use of nf-core practices

This repository is maintained under the personal `RasmusKoRiis` GitHub account. It is not an official nf-core pipeline, is not reviewed or released by the nf-core community, and must not be presented with an `nf-core/<pipeline>` repository identity.

The project still uses nf-core's technical conventions where they improve maintenance: DSL2 component directories, standard process labels, pinned containers, software-version channels, schema-driven parameters, nf-test, GitHub Actions, stable output publishing, and reference provenance. `nf-core pipelines lint` is used as an advisory technical audit.

Some lint checks assume an official repository generated from the newest template. They are not applicable here:

- `pipeline_name_conventions`, `manifest.name`, `manifest.homePage`, and expected nf-core logos require nf-core organization ownership and branding.
- AWS test workflows, branch automation, PR-comment actions, RO-Crate synchronization, and nf-core institutional configs require nf-core release infrastructure that this personal repository does not use.
- The legacy utility classes in `lib/` remain required while the workflow uses the nf-schema 2.4 compatibility line and Nextflow's v1 syntax parser.
- Public `test` and `test_full` profiles cannot contain operational influenza sequences or controlled databases. CI uses synthetic component fixtures, while the end-to-end acceptance procedure in [testing.md](testing.md) keeps private data and results outside Git.
- `template_strings` mistakes the literal Nextflow interpolation expressions constructed by Python unit tests for unrendered pipeline-template text.
- The nf-core template-version badge is omitted because this repository is not synchronized from a current nf-core pipeline template. The README identifies the nf-core tools version used for advisory audits instead.

The corresponding whole-test exceptions are listed explicitly in `.nf-core.yml`: `files_exist`, `files_unchanged`, `nextflow_config`, `multiqc_config`, `template_strings`, the AWS checks, pipeline naming, and RO-Crate README synchronization. Personal branding is retained in the manifest and MultiQC report instead of inserting false `nf-core` organization URLs merely to satisfy those checks.

Applicable lint findings remain defects to address. A change should not silence an actionable finding merely to reduce the count. The current migration priorities are adding accurate `meta.yml` and meaningful nf-test coverage to local components, then migrating the workflow to Nextflow's strict syntax parser with biological acceptance testing.

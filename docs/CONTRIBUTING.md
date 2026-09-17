# Contributing

Changes are accepted through pull requests to the personal `RasmusKoRiis/nf-core-fluseq` repository. This project follows relevant nf-core component and pipeline practices, but it is not an official community-owned nf-core pipeline.

Keep each change focused and describe its effect on the four analysis modes. Add new pipeline parameters to both `nextflow.config` and `nextflow_schema.json`, document new stable outputs in `docs/output.md`, and pin any new container image to an immutable version or digest.

Before opening a pull request, run the Python tests, Nextflow smoke checks, nf-test suite, Black, Prettier, and shell syntax checks described in [testing.md](testing.md). Changes to biological logic or report interpretation also require a private acceptance run with controlled references and representative data. Never commit operational sequence data, sample sheets, credentials, reference databases, Nextflow work directories, or acceptance results.

Local modules and subworkflows should use the nf-core directory layout, declare typed inputs and outputs, emit `versions.yml`, use standard process labels, include `meta.yml`, and provide an nf-test test when a compact synthetic fixture can exercise meaningful behavior.

Report security-sensitive problems privately to the repository maintainer instead of attaching operational data to a public issue.

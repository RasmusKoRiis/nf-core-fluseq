# nf-core-fluseq

[![Infrastructure CI](https://github.com/RasmusKoRiis/nf-core-fluseq/actions/workflows/ci.yml/badge.svg)](https://github.com/RasmusKoRiis/nf-core-fluseq/actions/workflows/ci.yml)
[![nf-test](https://github.com/RasmusKoRiis/nf-core-fluseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/RasmusKoRiis/nf-core-fluseq/actions/workflows/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core tools audit](https://img.shields.io/badge/nf--core_tools-4.1.0-24B064)](https://github.com/nf-core/tools/releases/tag/4.1.0)

A Nextflow DSL2 pipeline for analysing human and avian influenza A/B data from Nanopore FASTQ files or assembled FASTA sequences.

This is a personally maintained pipeline that follows relevant nf-core practices. It is not an official community-owned nf-core pipeline and is not hosted in the `nf-core` GitHub organization. The scope and intentional lint exceptions are documented in [use of nf-core practices](docs/nfcore_practices.md).

## Analysis modes

| `--file` value | Input                            | Main purpose                                               |
| -------------- | -------------------------------- | ---------------------------------------------------------- |
| `human-fastq`  | Sample sheet and FASTQ directory | Consensus generation and human-influenza reporting         |
| `human-fasta`  | Multi-FASTA                      | Human-influenza analysis from existing consensus sequences |
| `avian-fastq`  | Sample sheet and FASTQ directory | Consensus generation, avian genotyping and reporting       |
| `avian-fasta`  | Multi-FASTA                      | Avian analysis from existing consensus sequences           |

The pipeline requires Nextflow 24.10.2 or newer and is normally run with Docker. The FHI routine environment uses `-profile docker,server`. Nextflow 26.04 and newer must currently use `NXF_SYNTAX_PARSER=v1`; the routine wrappers set this compatibility mode automatically.

## Quick start

Inspect the complete parameter help:

```bash
nextflow run RasmusKoRiis/nf-core-fluseq --helpFull
```

For infrastructure-branch evaluation, select the branch explicitly:

```bash
nextflow run RasmusKoRiis/nf-core-fluseq \
  -r infrastructure \
  -profile docker,server \
  -params-file run-parameters.yml
```

Use a tagged release instead of a moving branch for production once the infrastructure work is merged and released.

Example human FASTQ parameters:

```yaml
file: human-fastq
input: /data/run/samplesheet.csv
samples_dir: /data/run/fastq_pass
outdir: /results/INF001
runid: INF001
ha_database: /references/human_HA.fasta
na_database: /references/human_NA.fasta
inhibition_mutation_db: /references/Inhibtion_Mutations_of_Intrest_2324.xlsx
reassortment_database: /references/reassortment_database.fasta
sequence_references: /references/sequence_references
nextclade_dataset: /references/nextclade_datasets
```

## Routine-wrapper compatibility

The operational wrappers in `/home/rasmuskopperud.riis/Coding/flu-wrappers` remain supported. Their historical parameter names are translated to the canonical names with a deprecation warning, and all four report processes still publish CSV files to `<outdir>/reporthuman/`.

| Historical parameter      | Canonical parameter        |
| ------------------------- | -------------------------- |
| `--samplesDir`            | `--samples_dir`            |
| `--seq_quality_thershold` | `--seq_quality_threshold`  |
| `--mamalian_mutation_db`  | `--mammalian_mutation_db`  |
| `--inhibtion_mutation_db` | `--inhibition_mutation_db` |

See [usage](docs/usage.md), [reference-data requirements](docs/reference_data.md), [routine-wrapper contract](docs/routine_wrappers.md), [outputs](docs/output.md), and [testing and private acceptance](docs/testing.md).

## Reproducibility and tests

- Runtime containers are versioned or digest-pinned.
- Nextclade uses the supplied local datasets rather than downloading the latest dataset during each sample task.
- Every run publishes `pipeline_info/reference_manifest.tsv` with SHA-256 hashes of staged reference files.
- CI parses the complete workflow and validates routine and canonical FASTQ sheets on Nextflow 24.10.2 and current stable.
- Python unit tests and nf-test coverage for all four surveillance-summary modes are included.
- A guarded local acceptance script can use private operational data without placing it or its results in the repository.

The full biological workflow still requires controlled influenza databases and representative run data; those are intentionally not embedded in this source repository.

## Credits

nf-core-fluseq was written by Rasmus Kopperud Riis.

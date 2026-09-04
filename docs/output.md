# nf-core-fluseq: Output

Only stable, user-facing files are published from the Nextflow work directory. Paths below are relative to `--outdir`; not every directory is produced in every mode.

| Directory        | Contents                                                                                             |
| ---------------- | ---------------------------------------------------------------------------------------------------- |
| `reporthuman/`   | Final human or avian CSV report. This path is a compatibility contract used by all routine wrappers. |
| `report/`        | Focused human FASTA drug-resistance report when `--drug_resistance_only` is enabled.                 |
| `irma/`          | IRMA consensus and run outputs for FASTQ modes.                                                      |
| `subtyping/`     | HA/NA subtype calls and supporting results.                                                          |
| `genotyping/`    | Avian or full FASTA genotype results.                                                                |
| `reassortment/`  | Reassortment analysis outputs.                                                                       |
| `nextclade/`     | Nextclade CSV and sequence outputs.                                                                  |
| `subclade/`      | Seasonal subclade nomenclature results.                                                              |
| `flumut/`        | FluMut markers, mutation, literature, and converted files.                                           |
| `genin2/`        | GenIn2 results and the slim report used by avian reporting.                                          |
| `surveillance/`  | Segment QC, reassortment, resistance, and sample-level surveillance summaries.                       |
| `fastqc/`        | Raw-read FastQC results for FASTQ modes.                                                             |
| `multiqc/`       | MultiQC report and data for FASTQ modes.                                                             |
| `pipeline_info/` | Execution metadata, software versions, parameters, and reference checksums.                          |

## Provenance

`pipeline_info/reference_manifest.tsv` contains a SHA-256 digest for every file staged from the run's controlled reference inputs. Together with the parameter dump and pipeline revision, this identifies the data and code used for the run.

Nextflow also writes timestamped execution reports, timelines, traces, and DAGs under `pipeline_info/`. These files are useful for run acceptance, performance tuning, and incident investigation.

## Publishing behavior

The default `--publish_dir_mode` is `copy`. `move` is intentionally disallowed because moving task outputs can invalidate Nextflow resume behavior. Intermediate files remain in the configured work directory and are not copied into results unless explicitly listed in the output contract.

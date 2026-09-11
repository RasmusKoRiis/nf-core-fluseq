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

## Final reports

All four full analysis modes publish their final CSV below `reporthuman/` because routine wrappers depend on that location. The filename is derived from `--runid` for the human and FASTA reports; the avian FASTQ merge currently publishes `fluseq_merged_report.csv`. Human FASTQ reporting also emits a filtered FASTA in the task work directory, but it is not part of the stable published-output contract.

The final report combines sample identity, subtype, segment coverage, IRMA statistics, Nextclade assignments, subclade nomenclature, reassortment evidence, and mutation annotations when those results apply. See [report_data_dictionary.md](report_data_dictionary.md) for the column-level contract and missing-value rules.

With `--file human-fasta --drug_resistance_only`, `report/<runid>_drug_resistance_report.csv` is the focused antiviral-resistance result.

## Analysis evidence

- `irma/` contains per-sample consensus FASTA, BAM/BAI, VCF, read-count tables, allele tables, figures, amended consensus files, and IRMA secondary outputs for FASTQ modes.
- `subtyping/` contains HA/NA BLAST hit tables, subtype calls, status records, and diagnostic error text when classification is incomplete.
- `genotyping/` and `reassortment/` contain per-sample CSV/TSV classifications and the segment-level evidence used by avian and full FASTA reporting.
- `nextclade/` contains raw Nextclade CSV output, translated coding sequences, filtered mutation tables, summaries, and normalized amino-acid mutation CSV files.
- `subclade/` contains seasonal nomenclature calls. These depend on the pinned rules and Nextclade results used for the run.
- `flumut/` and `genin2/` contain tool-native annotations plus the converted or trimmed CSV files consumed by final reporting.
- `surveillance/` contains `segment_qc.tsv`, `reassortment_summary.tsv`, `resistance_summary.tsv`, and `sample_summary.json`.
- `fastqc/` and `multiqc/` contain raw-read QC reports for FASTQ modes. MultiQC summarizes technical evidence; it does not replace review of failed tasks or the final influenza report.

Missing mode-specific directories can be expected when the corresponding analysis is not selected. A missing required final report is a failed acceptance condition.

## Provenance

`pipeline_info/reference_manifest.tsv` contains a SHA-256 digest for every file staged from the run's controlled reference inputs. Together with the parameter dump and pipeline revision, this identifies the data and code used for the run.

Nextflow also writes timestamped execution reports, timelines, traces, and DAGs under `pipeline_info/`. These files are useful for run acceptance, performance tuning, and incident investigation.

## Publishing behavior

The default `--publish_dir_mode` is `copy`. `move` is intentionally disallowed because moving task outputs can invalidate Nextflow resume behavior. Intermediate files remain in the configured work directory and are not copied into results unless explicitly listed in the output contract.

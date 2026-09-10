# nf-core-fluseq: Usage

## Select one analysis mode

Exactly one workflow is selected with `--file`: `human-fastq`, `human-fasta`, `avian-fastq`, or `avian-fasta`. The default is `human-fastq`.

FASTQ modes require `--input` and `--samples_dir`. FASTA modes require `--fasta`. Run `nextflow run RasmusKoRiis/nf-core-fluseq --helpFull` for the generated parameter reference.

## FASTQ sample sheets

The routine barcode layout is supported directly:

```csv
PCR-PlatePosition,SequenceID,Barcode,KonsCt
A1,SAMPLE01,barcode01,24.56
```

`Barcode` is resolved below `--samples_dir`, and all `.fastq.gz` or `.fq.gz` files in that barcode directory are assigned to the sample.

A canonical explicit-path layout is also supported:

```csv
sample,fastq_1,fastq_2
SAMPLE01,/data/run/barcode01/chunk.fastq.gz,
```

Sample identifiers must be unique and non-empty. Whitespace is converted to underscores. Every resolved FASTQ must exist, and a sample with no reads stops the run before analysis begins.

## Controlled reference inputs

| Parameter                        | Human FASTQ | Human FASTA                              | Avian FASTQ/FASTA |
| -------------------------------- | ----------- | ---------------------------------------- | ----------------- |
| `--ha_database`, `--na_database` | Required    | Required                                 | Required          |
| `--sequence_references`          | Required    | Required                                 | Required          |
| `--nextclade_dataset`            | Required    | Required                                 | Required          |
| `--inhibition_mutation_db`       | Required    | Required                                 | Required          |
| `--reassortment_database`        | Required    | Required unless `--drug_resistance_only` | Required          |
| `--genotype_database`            | Not used    | Required unless `--drug_resistance_only` | Required          |
| `--mammalian_mutation_db`        | Not used    | Not used                                 | Required          |

The large and season-controlled reference sets have no hidden repository defaults. Supply them explicitly or through a version-controlled parameters file. See [reference_data.md](reference_data.md) for layout and update guidance.

## Example commands

Human FASTQ:

```bash
nextflow run RasmusKoRiis/nf-core-fluseq \
  -r <release-or-branch> \
  -profile docker,server \
  --file human-fastq \
  --input samplesheet.csv \
  --samples_dir fastq_pass \
  --outdir results \
  --sequence_references /references/sequence_references \
  --nextclade_dataset /references/nextclade_datasets \
  --reassortment_database /references/reassortment_database.fasta \
  --inhibition_mutation_db /references/Inhibtion_Mutations_of_Intrest_2324.xlsx
```

Avian FASTA additionally supplies `--fasta`, `--genotype_database`, and `--mammalian_mutation_db`, and changes `--file` to `avian-fasta`.

For routine operation, prefer a YAML parameters file and keep it with the run record. Pipeline parameters belong in `-params-file`; executor, storage, and resource settings belong in a Nextflow config supplied with `-c`.

## Profiles and work files

All modules inherit `errorStrategy = 'ignore'`: failed tasks are logged and skipped while the remaining tasks continue, without automatic retries. Outputs that depend on a failed task may be missing for that sample. Check the execution trace and log for failures before accepting results. Input validation and workflow-level errors still stop the run.

- `docker` enables Docker and pulls missing pinned images.
- `server` sets the routine work directory to `/mnt/tempdata/work_fluseq`, preserves the work directory for resume/debugging, and caps tasks at 16 CPUs, 256 GB RAM, and 20 hours.
- Override the server work location with `--server_work_dir` when necessary.
- Use `-resume` after an interrupted run. Do not delete the work directory until the result has been accepted and archived.

Nextflow 26.04 enables its new strict syntax parser by default. This pipeline still uses dynamic DSL2/Groovy constructs, so set `NXF_SYNTAX_PARSER=v1` with Nextflow 26.04 or newer. The routine wrappers and CI set it automatically. Migrating the scientific workflows to the strict parser should be handled as a separately validated change because it touches channel-construction logic throughout the pipeline.

No fake `test` profile is shipped. Infrastructure smoke tests live under `tests/input_check`, and module tests under `tests/modules`. End-to-end biological acceptance must use the controlled routine references and representative influenza data.

## Legacy wrapper names

The routine wrappers may continue to pass `--samplesDir`, `--mamalian_mutation_db`, and `--inhibtion_mutation_db`. These aliases are translated before schema validation. New scripts and parameter files should use the correctly spelled canonical names.

## Reproducible execution

Pin a release or immutable commit with `-r`. `nextflow pull` only updates the local cache; the selected `-r` value controls the code that is executed. Preserve these files with each accepted run:

- the parameters file and exact launch command;
- `pipeline_info/reference_manifest.tsv`;
- execution report, timeline, trace, and DAG;
- software-version output;
- the pipeline Git revision or release tag.

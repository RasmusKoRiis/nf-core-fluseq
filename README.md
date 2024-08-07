# nf-core/fluseq :sneezing_face: 

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

This pipeline processes FASTQ files from Nanopore sequencing of Influenza A and B, generating consensus sequences and analyzing them for mutations, sequencing statistics, and drug resistance effects. The main steps include:

- Alignment and consensus sequencing with IRMA.
- Consensus sequence analysis with Nextclade.
- Mutation calling for all segments.
- Generation of a comprehensive report in CSV format.
- Output of high-quality sequences in a multiple FASTA file.

## Compatibility

- **Operating System**: Linux
- **Dependencies**: Docker and Nextflow

## Usage

### Sample Sheet Preparation

Prepare a sample sheet (CSV or TSV) in the `assets` folder with the following format:

```
PCR-PlatePosition,SequenceID,Barcode,KonsCt
A1*,sampleID,barcodeID,ct-value*
```
*not compulsory

Each row lists a sample to be analyzed. Samples not listed in the sheet will be excluded from the analysis.

### Directory Structure

Ensure your directory structure is as follows:

```
./
  |-data
         |-barcode1
               |-XXXX_pass_barcode01_XXXX.fastq.gz
               |-YYYY_pass_barcode01_YYYY.fastq.gz
  |-nf-core-fluseq
               |-assets
                     |-samplesheet.csv
                     |-samplesheet.tsv
               |-...
```

### Running the Pipeline

Navigate to the `nf-core-fluseq` folder and execute the following command with default parameters:

```bash
nextflow run main.nf -profile docker --runid runid_name --outdir ../outdir_name
```

### Important Parameters

- `--input` (default: `assets/samplesheet.csv`): Path to the samplesheet.
- `--seq_quality_threshold` (default: 95): Coverage threshold for consensus sequences.
- `--samplesDir` (default: `../data`): Directory containing the FASTQ files.

All parameters are detailed in the `nextflow.config` file.

## Pipeline Output

The output includes:

- Consensus sequences.
- Mutation calls.
- Sequencing statistics (coverage, quality parameters).
- Drug resistance effects.
- A report in CSV format.
- A multiple FASTA file of sequences that passed quality filters.


## Credits

nf-core/fluseq was originally written by Rasmus Kopperud Riis.


# nf-core/fluseq :sneezing_face: 

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

This pipeline processes FASTQ files from Nanopore sequencing of Influenza A and B, generating consensus sequences and analyzing them for mutations, sequencing statistics, and drug resistance effects. The main steps include:

- Alignment and consensus sequencing with IRMA.
- Consensus sequence analysis with Nextclade.
- Mutation calling for all segments.
- Generation of a comprehensive report in CSV format.
- Output of sequences in a multiple FASTA file.

The pipeline consist of four different worflows listed bellow:

1) Human Influenza FASTQ analysis (human)
  Alignment of FASTQ and mutation analysis 
2) Human Influenza FASTA analysis (human-fasta) (under development)
   Mutation analysis 
3) Avian Influenza FASTQ analysis (avian)
  Alignment of FASTQ and mutation analysis 
4) Avian Influenza FASTA analysis (avian-fasta)
   Mutation analysis 

## Compatibility

- **Operating System**: Linux
- **Dependencies**: Docker and Nextflow

## Usage

### Sample Sheet Preparation

Prepare a sample sheet (CSV or TSV*) in the `assets` folder with the following format:
* TSV file is not not compulsory

```
PCR-PlatePosition,SequenceID,Barcode,KonsCt
A1*,sampleID,barcodeID,ct-value*
```
*not compulsory

Each row lists a sample to be analyzed. Samples not listed in the sheet will be excluded from the analysis.

### Directory Structure

#### For FASTQ-analysis
Ensure your directory structure is as follows:

```
./
  |-data
         |-barcode3
               |-XXXX_pass_barcode03_XXXX.fastq.gz
               |-YYYY_pass_barcode03_YYYY.fastq.gz
  |-nf-core-fluseq
               |-assets
                     |-samplesheet.csv
                     |-samplesheet.tsv
               |-...
```

### Running the Pipeline

Navigate to the `nf-core-fluseq` folder and execute the following command with default parameters:

#### Human Influenza FASTQ analysis

```bash
nextflow run main.nf -profile docker --runid runid_name --outdir ../outdir_name
```

#### Avian Influenza FASTQ analysis

```bash
nextflow run main.nf -profile docker --file avian-fastq  --genotype_database database* --runid runid_name --outdir ../outdir_name
```
* The database given as the genotyping database must be in the format given bellow:
 ``` 
>DatabaseNumber_|Subtype|ID|Segment|SegmentNumber|GISAIDID
aa..
>DatabaseNumber_|Subtype|ID|Segment|SegmentNumber|GISAIDID
a..
```

Example of header:
```
>21_|H5N8|chicken/norway|HA|4|EPI_ESL_7473825
```

The database number is used in the genotyping02.py script to identify genotypes. Either the offical database has to be obtained or this script has to be adjusted to a be compatible to a in-house genotyping database.

#### Avian Influenza FASTA analysis

```bash
nextflow run main.nf -profile docker --file avian-fasta  --genotype_database database --runid runid_name --outdir ../outdir_name
```

### Important Parameters

- `--input` (default: `assets/samplesheet.csv`): Path to the samplesheet.
- `--seq_quality_threshold` (default: 20): Coverage threshold for analysis of consensus sequences.
- `--samplesDir` (default: `../data`): Directory containing the FASTQ files in the structure given above.

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

fluseq was originally written by Rasmus Kopperud Riis.

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

#### Human Influenza FASTA drug-resistance analysis only

Use the wrapper when only the drug-resistance result is needed:

```bash
bash bin/fasta_drug_resistance_wrapper.sh runid_name ../outdir_name input.fasta
```

The wrapper runs the required FASTA parsing, HA/NA subtyping and amino-acid
translation steps, but skips genotyping, reassortment, coverage, clade analysis,
surveillance summaries and the full report. Results are written below the
requested output directory. The focused report is written to
`report/<runid>_drug_resistance_report.csv`, and the individual lookup CSVs are
available in `tablelookup/`.
Set `FLUSEQ_PROFILE` to use a profile other than Docker, and append `-resume` to
reuse completed work:

```bash
FLUSEQ_PROFILE=apptainer bash bin/fasta_drug_resistance_wrapper.sh runid_name ../outdir_name input.fasta -resume
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

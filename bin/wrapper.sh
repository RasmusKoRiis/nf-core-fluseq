#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <runid_name> <outdir_name> <data_dir>"
    exit 1
fi

# Assign arguments to variables
runid=$1
outdir=$2
data_dir=$3

# Define the repository URL
repo_url="https://github.com/RasmusKoRiis/nf-core-fluseq"

# Clone the repository if it doesn't exist
if [ ! -d "nf-core-fluseq" ]; then
    echo "Cloning the repository from $repo_url..."
    git clone "$repo_url"
else
    echo "Repository already exists. Skipping clone."
fi


# Check if the samplesheet.csv exists and convert to samplesheet.tsv
if [ -f "samplesheet.csv" ]; then
    echo "Converting samplesheet.csv to samplesheet.tsv..."
    # Replace commas with tabs to create a TSV copy
    sed 's/,/\t/g' samplesheet.csv > samplesheet.tsv
else
    echo "samplesheet.csv not found!"
    exit 1
fi

# Move both CSV and TSV files to the assets directory
if [ -d "nf-core-flu/assets" ]; then
    echo "Moving samplesheet.csv and samplesheet.tsv to the assets directory..."
    mv samplesheet.csv nf-core-flu/assets/
    mv samplesheet.tsv nf-core-flu/assets/
else
    echo "assets directory not found! Creating assets directory..."
    mkdir -p assets
    mv samplesheet.csv nf-core-flu/assets/
    mv samplesheet.tsv nf-core-flu/assets/
fi

# Change to the cloned repository directory
cd nf-core-fluseq || { echo "Failed to enter nf-core-fluseq directory"; exit 1; }

# Find the fastq_pass folder inside the data_dir
fastq_pass_path=$(find "$data_dir" -type d -name "fastq_pass" 2>/dev/null | head -n 1)

# Check if the fastq_pass folder was found
if [ -n "$fastq_pass_path" ]; then
    echo "Found fastq_pass folder at $fastq_pass_path. Copying and renaming it to 'data'..."
    cp -r "$fastq_pass_path" ../data
else
    echo "fastq_pass folder not found in $data_dir!"
    exit 1
fi


# Run the nextflow command with provided arguments
nextflow run main.nf -profile docker --runid "$runid" --outdir "$outdir"

#!/usr/bin/env bash

set -euo pipefail

if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 <runid_name> <outdir_name> <fasta_file> [additional Nextflow arguments]" >&2
    echo "Example: $0 run42 results/run42 samples.fasta -resume" >&2
    exit 1
fi

runid=$1
outdir=$2
fasta_file=$3
shift 3

if [[ ! -f "$fasta_file" ]]; then
    echo "FASTA file not found: $fasta_file" >&2
    exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Nextflow is not available on PATH." >&2
    exit 1
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "$script_dir/.." && pwd)
profile=${FLUSEQ_PROFILE:-docker}

exec nextflow run "$project_dir/main.nf" \
    -profile "$profile" \
    --file human-fasta \
    --drug_resistance_only \
    --fasta "$fasta_file" \
    --runid "$runid" \
    --outdir "$outdir" \
    "$@"

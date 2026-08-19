#!/usr/bin/env bash

set -euo pipefail

if [[ "$#" -lt 4 ]]; then
    echo "Usage: $0 <runid_name> <results_dir> <work_dir> <fasta_file> [additional Nextflow arguments]" >&2
    echo "Example: $0 run42 results work samples.fasta -resume" >&2
    exit 1
fi

runid=$1
requested_outdir=$2
requested_workdir=$3
fasta_file=$4
shift 4

if [[ ! -f "$fasta_file" ]]; then
    echo "FASTA file not found: $fasta_file" >&2
    exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Nextflow is not available on PATH." >&2
    exit 1
fi

launch_dir=$(pwd -P)
fasta_dir=$(cd -- "$(dirname -- "$fasta_file")" && pwd -P)
fasta_file="$fasta_dir/$(basename -- "$fasta_file")"

make_absolute_dir() {
    local requested_dir=$1
    if [[ "$requested_dir" != /* ]]; then
        requested_dir="$launch_dir/$requested_dir"
    fi
    mkdir -p -- "$requested_dir"
    cd -- "$requested_dir" && pwd -P
}

results_dir=$(make_absolute_dir "$requested_outdir")
work_dir=$(make_absolute_dir "$requested_workdir")

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
local_project=$(cd -- "$script_dir/.." && pwd)
repo_url=${FLUSEQ_REPO_URL:-https://github.com/RasmusKoRiis/nf-core-fluseq}

if [[ -f "$local_project/main.nf" && -f "$local_project/workflows/human-fasta.nf" ]]; then
    project_dir=$local_project
else
    project_dir=${FLUSEQ_REPO_DIR:-$launch_dir/nf-core-fluseq}
    if [[ ! -f "$project_dir/main.nf" ]]; then
        if [[ -e "$project_dir" ]]; then
            echo "Repository path exists but is not an nf-core-fluseq checkout: $project_dir" >&2
            exit 1
        fi
        if ! command -v git >/dev/null 2>&1; then
            echo "Git is required to download nf-core-fluseq." >&2
            exit 1
        fi
        echo "Downloading nf-core-fluseq into $project_dir ..."
        if [[ -n "${FLUSEQ_REPO_REF:-}" ]]; then
            git clone --branch "$FLUSEQ_REPO_REF" "$repo_url" "$project_dir"
        else
            git clone "$repo_url" "$project_dir"
        fi
    else
        echo "Using existing repository: $project_dir"
    fi
fi

profile=${FLUSEQ_PROFILE:-docker}

echo "Persistent work directory: $work_dir"
echo "Results directory: $results_dir"

exec nextflow run "$project_dir/main.nf" \
    -profile "$profile" \
    -work-dir "$work_dir" \
    --file human-fasta \
    --drug_resistance_only \
    --fasta "$fasta_file" \
    --runid "$runid" \
    --outdir "$results_dir" \
    "$@"

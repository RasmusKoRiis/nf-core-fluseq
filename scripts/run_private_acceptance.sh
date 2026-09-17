#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

: "${FLUSEQ_ACCEPTANCE_DATA:?Set FLUSEQ_ACCEPTANCE_DATA to a private run-data directory outside this repository}"
: "${FLUSEQ_ACCEPTANCE_PARAMS:?Set FLUSEQ_ACCEPTANCE_PARAMS to a private YAML parameters file containing controlled reference paths}"
: "${FLUSEQ_ACCEPTANCE_OUTDIR:?Set FLUSEQ_ACCEPTANCE_OUTDIR to a result directory outside this repository}"

data_dir=$(realpath "$FLUSEQ_ACCEPTANCE_DATA")
params_file=$(realpath "$FLUSEQ_ACCEPTANCE_PARAMS")
out_dir=$(realpath -m "$FLUSEQ_ACCEPTANCE_OUTDIR")
profile=${FLUSEQ_ACCEPTANCE_PROFILE:-docker}
work_dir=${FLUSEQ_ACCEPTANCE_WORKDIR:-$(mktemp -d /tmp/fluseq-acceptance-work.XXXXXX)}

for private_path in "$data_dir" "$params_file" "$out_dir" "$work_dir"; do
    case "$private_path" in
        "$repo_dir"|"$repo_dir"/*)
            echo "ERROR: private acceptance data, parameters, work, and results must remain outside $repo_dir" >&2
            exit 1
            ;;
    esac
done

mapfile -t sample_sheets < <(find "$data_dir" -maxdepth 2 -type f -name '*.csv' -print)
mapfile -t fastq_dirs < <(find "$data_dir" -type d -name fastq_pass -print)

if [[ ${#sample_sheets[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one CSV sample sheet within two levels of the private data root; found ${#sample_sheets[@]}" >&2
    exit 1
fi
if [[ ${#fastq_dirs[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one fastq_pass directory under the private data root; found ${#fastq_dirs[@]}" >&2
    exit 1
fi

mkdir -p "$out_dir" "$work_dir"

nextflow run "$repo_dir" \
    -profile "$profile" \
    -params-file "$params_file" \
    -work-dir "$work_dir" \
    --file human-fastq \
    --input "${sample_sheets[0]}" \
    --samples_dir "${fastq_dirs[0]}" \
    --outdir "$out_dir"

if ! find "$out_dir/reporthuman" -maxdepth 1 -type f -name '*.csv' -print -quit | grep -q .; then
    echo "ERROR: acceptance run completed without a final report CSV" >&2
    exit 1
fi

echo "Private acceptance test passed. Results remain outside the repository: $out_dir"

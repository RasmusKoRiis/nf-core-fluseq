process EMIT_FASTA_RECORD {
    tag { sample_id }
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

    input:
    tuple val(sample_id), val(header_id), val(file_stem), val(seq), val(orig_name)

    output:
    tuple val(sample_id), val(orig_name), path("${file_stem}.fasta"), emit: fasta
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    cat > ${file_stem}.fasta <<'EOF'
    >${header_id}
    ${seq}
    EOF

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        alpine: 3.20.3
    END_VERSIONS
    """
}

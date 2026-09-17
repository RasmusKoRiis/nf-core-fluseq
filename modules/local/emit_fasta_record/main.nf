process EMIT_FASTA_RECORD {
    tag { sample_id }
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    tuple val(sample_id), val(header_id), val(file_stem), val(seq), val(orig_name)

    output:
    tuple val(sample_id), val(orig_name), path("${file_stem}.fasta"), emit: fasta
    path 'versions.yml', emit: versions, optional: true

    script:
    """
    set -euo pipefail
    cat > ${file_stem}.fasta <<'EOF'
    >${header_id}
    ${seq}
    EOF

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        bash: 5.2.15
    END_VERSIONS
    """
}

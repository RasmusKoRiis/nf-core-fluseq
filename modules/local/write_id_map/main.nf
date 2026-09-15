process WRITE_ID_MAP {
    tag 'id_map'
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    val pairs

    output:
    path 'id_map.tsv', emit: id_map
    path 'versions.yml', emit: versions, optional: true

    script:
    def lines = pairs.collect { row ->
        def uid = (row[0] ?: '').toString().trim()
        def original = (row[1] ?: '').toString().trim()
        "${uid}\t${original}"
    }.join('\n')
    """
    set -euo pipefail
    printf 'SampleID\\tOriginalName\\n' > id_map.tsv
    cat >> id_map.tsv <<'EOF'
    ${lines}
    EOF

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        bash: 5.2.15
    END_VERSIONS
    """
}

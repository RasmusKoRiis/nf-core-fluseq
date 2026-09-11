process WRITE_ID_MAP {
    tag 'id_map'
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

    input:
    val pairs

    output:
    path 'id_map.tsv', emit: id_map
    path 'versions.yml', emit: versions

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
        alpine: 3.20.3
    END_VERSIONS
    """
}

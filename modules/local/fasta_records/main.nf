process EMIT_FASTA_RECORD {
    tag { sample_id }
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

    input:
    tuple val(sample_id), val(header_id), val(file_stem), val(seq), val(orig_name)

    output:
    tuple val(sample_id), val(orig_name), path("${file_stem}.fasta")
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

process REHEADER_TO_UID {
    tag { meta.id }
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

    input:
    tuple val(meta), path(ha), path(na)

    output:
    tuple val(meta), path("HA_${meta.id}.fa"), path("NA_${meta.id}.fa")
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${ha} > HA_${meta.id}.fa
    awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${na} > NA_${meta.id}.fa

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        alpine: 3.20.3
    END_VERSIONS
    """
}

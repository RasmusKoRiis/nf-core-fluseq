process REHEADER_TO_UID {
    tag { meta.id }
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

    input:
    tuple val(meta), path(ha), path(na)

    output:
    tuple val(meta), path("HA_${meta.id}.fa"), path("NA_${meta.id}.fa"), emit: fasta
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

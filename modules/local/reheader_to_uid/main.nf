process REHEADER_TO_UID {
    tag { meta.id }
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    tuple val(meta), path(ha), path(na)

    output:
    tuple val(meta), path("HA_${meta.id}.fa"), path("NA_${meta.id}.fa"), emit: fasta
    path 'versions.yml', emit: versions, optional: true

    script:
    """
    set -euo pipefail
    awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${ha} > HA_${meta.id}.fa
    awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${na} > NA_${meta.id}.fa

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        bash: 5.2.15
    END_VERSIONS
    """
}

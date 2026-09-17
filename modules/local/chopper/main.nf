process CHOPPER {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_filtered.fastq") , emit: chopperfastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    
    chopper -q 11 -l 400 --tailcrop 22 --headcrop 22  -i $fastq > ${meta.id}_filtered.fastq

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | head -n 1)
    END_VERSIONS

    """
}

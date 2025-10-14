process CHOPPER {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*fastq") , emit: chopperfastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    chopper -q 11 -l 400 --tailcrop 22 --headcrop 22  -i $fastq > ${meta.id}_filtered.fastq

    """
}

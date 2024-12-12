process AMPLIGONE {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)

    output:
    tuple val(meta), path("*fastq") , emit: primertrimmedfastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    ampligone \
        --input $fastq\
        --output ${meta.id}.fastq \
        --reference $primerdir/SARS-CoV-2.reference.fasta \
        --primers $primerdir/SARS-CoV-2.scheme.bed \
        --export-primers ${meta.id}_removed_coordinates.bed \
        --amplicon-type end-to-end
    """
}
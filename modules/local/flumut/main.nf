
process FLUMUT {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
   


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //conda "bioconda::irma=1.0.3"
    //container 'docker.io/rasmuskriis/cdc_irma_custom:1.0'
    container 'quay.io/biocontainers/flumu:0.6.3--pyhdfd78af_0'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0':
        //'biocontainers/irma:1.0.3--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
   

    output:
    tuple val(meta), path("${meta.id}_markers_output.tsv") , emit: markers
    tuple val(meta), path("${meta.id}_mutations_output.tsv"),  path("$meta.id/*.bai"), emit: mutations
    tuple val(meta), path("${meta.id}_literature_output.tsv") , emit: literature


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    flumut --update
    flumut -m ${meta.id}_markers_output.tsv -M ${meta.id}_mutations_output.tsv -l ${meta.id}_literature_output.tsv $fasta
    """
}

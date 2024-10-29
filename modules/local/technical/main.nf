
process TECHNICAL {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
   

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(irma_stat)


    output:
    tuple val(meta), path("*irmastat.csv"), emit: depth_files
    path("*irmastat.csv"), emit: depth_files_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #DEPTH CACULATION
    # Make output name
    output_name="${meta.id}_irmastat.csv"
    
    python /project-bin/sequence_quality.py \
        $irma_stat \
        \$output_name \
        ${meta.id} 



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

}

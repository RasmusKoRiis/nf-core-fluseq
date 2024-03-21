
process SEGMENTIFENTIFIER  {
    tag "$meta.id"
    label 'process_single'
    debug true

    conda "bioconda::blast=2.15.0"
    container 'docker.io/ncbi/blast:latest'

    input:
    tuple val(meta), path(fasta)
    path(reference)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta_segment
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
    for fasta_file in ${fasta}; do

    basename=\$(basename \$fasta_file)
    blastn -query \$fasta_file -subject $reference -outfmt 6   > "\${fasta_file}_${prefix}.tsv"
    segment=\$(awk -F '\t' '{split(\$2,a,"_"); print a[3]}' "\${fasta_file}_${prefix}.tsv" | head -n 1)
    mv \$fasta_file "\${basename}"_\${segment}.fasta
    
    done

   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

}

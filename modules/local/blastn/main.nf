
process SUBTYPEFINDER {
    tag "$meta.id"
    label 'process_single'
   

    conda "bioconda::blast=2.15.0"
    container 'docker.io/ncbi/blast:latest'

    input:
    tuple val(meta), path(ha_fasta), path(na_fasta),path(pa_fasta), path(pb1_fasta),path(pb2_fasta), path(np_fasta),path(ns_fasta), path(m_fasta)
    path(ha_database)
    path(na_database)

    output:
    tuple val(meta), path("*subtype.txt"), emit: subtype
    tuple val(meta), path("*subtype.csv"), emit: subtype_file
    path("*subtype.csv"), emit: subtype_report
    tuple val(meta), path("*.tsv"), emit: subtype_file_test
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
    blastn -query $ha_fasta -subject $ha_database -outfmt 6   > "ha_${prefix}.tsv"
    subtype_ha=\$(awk -F '\t' '{split(\$2,a,\"_\"); print a[2]}' "ha_${prefix}.tsv" | head -n 1)
    
    blastn -query $na_fasta -subject $na_database -outfmt 6   > "na_${prefix}.tsv"
    subtype_na=\$(awk -F '\t' '{split(\$2,a,\"_\"); print a[2]}' "na_${prefix}.tsv" | head -n 1)

    subtype=""$meta.id",\$subtype_ha\$subtype_na"
    subtype_txt="\$subtype_ha\$subtype_na"

    echo \$subtype > ""$meta.id"_subtype.csv"
    echo "Sample,Subtype" | cat - "$meta.id"_subtype.csv > temp && mv temp "$meta.id"_subtype.csv

    echo \$subtype_txt > ""$meta.id"_subtype.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

}

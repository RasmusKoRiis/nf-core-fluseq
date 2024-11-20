
process SEGMENTIFENTIFIER  {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
   
    conda "bioconda::blast=2.15.0"
    container 'docker.io/ncbi/blast:latest'

    input:
    tuple val(meta), path(fasta)
    path(reference)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta_segment
    tuple val(meta), path("*.fa"), emit: fasta_configuration
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

    basename=\$(basename \$fasta_file .fasta)
    blastn -query \$fasta_file -subject $reference -outfmt 6   > "\${fasta_file}_${prefix}.tsv"
    segment=\$(awk -F '\t' '{split(\$2,a,"_"); print a[3]}' "\${fasta_file}_${prefix}.tsv" | head -n 1)

    # Assign a prefix based on the segment type
    case \$segment in
        PB2)
            prefix_number="1"
            ;;
        PB1)
            prefix_number="2"
            ;;
        PA)
            prefix_number="3"
            ;;
        HA)
            prefix_number="4"
            ;;
        NP)
            prefix_number="5"
            ;;
        NA)
            prefix_number="6"
            ;;
        MP)
            prefix_number="7"
            ;;
        NS)
            prefix_number="8"
            ;;
        *)
            echo "Unknown segment type: \$segment"
            prefix_number="unknown"
            ;;
    esac

    # Version 1 of the fasta file with one naming convention
    new_name_1="\${prefix_number}_\${segment}_\${basename}"
    echo ">\${new_name_1}" > "\${segment}_\${basename}.fasta"
    tail -n +2 \$fasta_file >> "\${segment}_\${basename}.fasta"

    # Version 2 of the fasta file with another naming convention
    new_name_2="\${basename}_\${prefix_number}"
    echo ">\${new_name_2}" > "\${basename}_\${prefix_number}.fa"
    tail -n +2 \$fasta_file >> "\${basename}_\${prefix_number}.fa"

    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

}

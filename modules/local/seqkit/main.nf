
process FASTA_CONFIGURATION {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::seqkit=2.7.0"
    container "docker.io/nanozoo/seqkit:latest"

    input:
    tuple val(meta), path(fasta_files), path(subtype)
  
    output:
    tuple val(meta), path("*.fa"), path(subtype), emit: fasta
    tuple val(meta), path("*.fasta"), emit: fasta_merge
    


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
    for fasta_file in \$(ls ${fasta_files}); do
        id=\$(echo "\$fasta_file" | sed 's/_.*//')
        number=\$(echo "\$fasta_file" | sed 's/.*_//' | sed 's/\\.fa//' | xargs printf "%02d")
        subtype=\$(cat ${subtype} )

        segment=""
        case "\$number" in
            "04") segment="01-HA" ;;
            "06") segment="02-NA" ;;
            "07") segment="03-M" ;;
            "02") segment="04-PB1" ;;
            "01") segment="05-PB2" ;;
            "05") segment="06-NP" ;;
            "03") segment="07-PA" ;;
            "08") segment="08-NS" ;;
            *) segment="Unknown" ;;
        esac

        # Prepare the new header
        new_header=">\${id}|\${segment}-\${subtype}"

        # Rename the header in the fasta file
        # Use awk to process the file. If it's the first line (the header), replace it. Otherwise, print the line as it is.
        awk -v new_header="\$new_header" 'NR == 1 {print new_header; next} {print}' "\$fasta_file" > temp_file.fa

        # Replace the original file with the modified one
        mv temp_file.fa ${meta.id}_\${segment}-\${subtype}.fa

        cat ${meta.id}_\${segment}-\${subtype}.fa >> ${meta.id}.fasta      

        
    done


    """


}

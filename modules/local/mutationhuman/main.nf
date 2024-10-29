
process MUTATIONHUMAN  {
    tag "$meta.id"
    label 'process_single'
    //errorStrategy 'ignore'
   
    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(fasta), path(subtype)
    path(sequence_references)
   
    
    output:
    //tuple val(meta), path("*.txt"), emit: genotype
    //tuple val(meta), path("*.csv"), emit: genotype_file
    tuple val(meta), path("*human_mutation.csv"), path(subtype), emit: human_mutation, optional: true
    tuple val(meta), path("*inhibtion_mutation.csv"), path(subtype), emit: inhibtion_mutation, optional: true
    tuple val(meta), path("*vaccine_mutation.csv"), emit: vaccine_mutation, optional: true
    
    path("*human_mutation_report.csv"), emit: human_mutation_report, optional: true
    path("*inhibtion_mutation_report.csv"), emit: inhibtion_mutation_report, optional: true
    path("*vaccine_mutation_report.csv"), emit: vaccine_mutation_report, optional: true

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    //errorStrategy 'ignore'

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
    subtype_name=\$(cat ${subtype} )
  
    #Subtype fix for B (not optimal)
    if [[ "\$subtype_name" == "VIC" ]]; then
        subtype_name="VICVIC"
    fi
    
    for fasta_file in ${fasta}; do

        filename=\$(basename \$fasta_file)
        filename_no_ext=\${filename%.*}  

        # Extract 'A_XX_XX' assuming it's after the last '.' in the filename_no_ext
        segment_subtype=\$(echo "\${filename_no_ext}" | awk -F. '{print \$NF}')
  
        # Extract 'XX' as the segment assuming it's after the last '_' in segment_subtype
        segment=\$(echo "\${segment_subtype}" | awk -F_ '{print \$NF}')

        # Make output name
        output_name_human=${meta.id}_\${segment}"_human_mutation.csv"
        output_name_mamailian=${meta.id}_\${segment}"_mamailian_mutation.csv"
        output_name_inhibition=${meta.id}_\${segment}"_inhibtion_mutation.csv"
        output_name_vaccine=${meta.id}_\${segment}"_vaccine_mutation.csv"

        reference_file=${sequence_references}

        mamailian=mamailian
        inhibition=inhibition
        human=human
        inhibition_human=inhibition_human
        human_vaccine=human_vaccine

        # HUMAN MUTATIONS - DEFAULT

        if [[ "\${segment}" != *"NEP"* && "\${segment}" != *"PA-X"* && "\${segment}" != *"PB1-F2"* ]]; then 

            python /project-bin/mutation_finder.py \
                \$fasta_file \
                \$reference_file \
                \${segment} \
                \$subtype_name \
                \$output_name_human\
                \$human\

        else
            echo "Skipping execution for a segment containing 'NEP' or 'PA-X' in its name: \${segment}"
        fi
        

        # HUMAN MUTATIONS - INHIBTION

        if [[ "\${segment}" == *"NA"* || ( "\${segment}" == *"PA"* && "\${segment}" != *"PA-X"* ) || "\${segment}" == *"M2"* ]]; then

            python /project-bin/mutation_finder.py \
                \$fasta_file \
                \$reference_file \
                \${segment} \
                \$subtype_name \
                \$output_name_inhibition\
                \$inhibition_human

        else
            echo "Skipping execution for a segment not containing 'NA', 'M' or 'PA' in its name: \${segment}"
        fi

        # HUMAN MUTATIONS - VACCINE
        
        if [[ ("\${segment}" == *"HA"* || "\${segment}" == *"NA"*) ]]; then

        python /project-bin/mutation_finder.py \
            \$fasta_file \
            \$reference_file \
            \${segment} \
            \$subtype_name \
            \$output_name_vaccine\
            \$human_vaccine



        else
            echo "Skipping execution for a segmentcontaining 'HA' or'NA' in its name: \${segment}"
        fi
        
    done
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

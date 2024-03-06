
process MUTATION  {
    tag "$meta.id"
    label 'process_single'


    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(fasta), path(subtype)
    path(sequence_references)
   
    

    output:
    //tuple val(meta), path("*.txt"), emit: genotype
    //tuple val(meta), path("*.csv"), emit: genotype_file
    tuple val(meta), path("*.csv"), path(subtype), emit: mutation_list
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

    for fasta_file in ${fasta}; do

        filename=\$(basename \$fasta_file)
        echo "Filename: \${filename}"

        filename_no_ext=\${filename%.*}  
        echo "Filename no ext: \${filename_no_ext}"

        # Extract 'A_H5_HA' assuming it's after the last '.' in the filename_no_ext
        segment_subtype=\$(echo "\${filename_no_ext}" | awk -F. '{print \$NF}')
        echo "Segment subtype: \${segment_subtype}"

        # Extract 'HA' as the segment assuming it's after the last '_' in segment_subtype
        segment=\$(echo "\${segment_subtype}" | awk -F_ '{print \$NF}')

        # Make output name
        output_name=${meta.id}_\${segment}"_mutation.csv"

        reference_file=${sequence_references}
        
        echo "Segment: \${segment}"
        echo "Processing file: \${fasta_file}"
        echo "Reference: ${sequence_references}"
        echo "Subtype name: \${subtype_name}"
        echo "Output name: \${output_name}"

        python /project-bin/mutation_finder.py \
            \$fasta_file \
            \$reference_file \
            \${segment} \
            \$subtype_name \
            \$output_name\

        
    done
    




    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process TABLELOOKUP {
    tag "$meta.id"
    label 'process_single'
    
    

   

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(mamailian_mutation), path(subtype)
    tuple val(meta), path(inhibition_mutation)
    path(mutation_table)

    output:
    //tuple val(meta), path("*.txt"), emit: genotype
    //tuple val(meta), path("*.csv"), emit: genotype_file
    tuple val(meta), path("*.csv"), emit: _mamalianadpation_mutations
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
    
    for mutation_file in ${mamailian_mutation}; do

        filename=\$(basename \$mutation_file)
        echo "Filename: \${filename}"

        filename_no_ext=\${filename%.*}  
        echo "Filename no ext: \${filename_no_ext}"

        segment=\$(echo "\${filename_no_ext}" | awk -F_ '{print \$2}')

        # Make output name
        output_name=${meta.id}_\${segment}"_mamalianadpation.csv"

        # State mutation table
        segment_table=${mutation_table}/"mamailian_adpation/\${segment}.csv"
        
        echo "Segment: \${segment}"
        echo "Processing file: \${mutation_file}"
        echo "Subtype name: \${subtype_name}"
        echo "Output name: \${output_name}"

        type=mamailian

        if [ ! -f \$segment_table ]; then
            echo "Segment table does not exist: \${segment_table}"
   
        else
            echo "Segment table exists: \${segment_table}"
            python /project-bin/table_lookup.py \
            \$mutation_file\
            \$segment_table \
            \$output_name\
            \${segment} \
            ${meta.id} \
            \${type} 

        fi

       
    done

    for mutation_file in ${inhibition_mutation}; do

        filename=\$(basename \$mutation_file)
        echo "Filename: \${filename}"

        filename_no_ext=\${filename%.*}  
        echo "Filename no ext: \${filename_no_ext}"

        segment=\$(echo "\${filename_no_ext}" | awk -F_ '{print \$2}')

        # Make output name
        output_name=${meta.id}_\${segment}"_inhibtion.csv"

        # State mutation table
        segment_table=${mutation_table}/"inhibtion_adaptaion/\${segment}.csv"
        
        echo "Segment: \${segment}"
        echo "Processing file: \${mutation_file}"
        echo "Subtype name: \${subtype_name}"
        echo "Output name: \${output_name}"

        type=inhitbion

        if [ ! -f \$segment_table ]; then
            echo "Segment table does not exist: \${segment_table}"
   
        else
            echo "Segment table exists: \${segment_table}"
            python /project-bin/table_lookup.py \
            \$mutation_file\
            \$segment_table \
            \$output_name\
            \${segment} \
            ${meta.id} \
            \${type} 
        fi
       
    done

    



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

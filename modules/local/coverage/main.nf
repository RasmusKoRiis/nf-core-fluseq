
process COVERAGE {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
    debug true

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(sequences), path(subtype)
    val(seq_quality_thershold)
   

    output:
    tuple val(meta), path("*.csv"), emit: coverage
    path("*.csv"), emit: coverage_report
    tuple val(meta), path("*fasta"), path(subtype), emit: filtered_fasta
    tuple val(meta), path("*coverage.fa"), path(subtype), emit:  merged_filtered_fasta
    path "versions.yml", emit: versions

    path("*fa"), emit: filtered_fasta_report

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
    for fasta_file in ${sequences}; do
        filename=\$(basename \$fasta_file)
        echo "Processing \$filename"  
        filename_no_ext=\${filename%.*}  
        segment_subtype=\${filename_no_ext#*-} 
        segment=\${segment_subtype%-*}  
        subtype_name=\${segment_subtype#*-} 
        

        output_csv="${meta.id}_\${segment}_coverage.csv"
        python /project-bin/coverage_finder.py "\$fasta_file" \$output_csv ${meta.id} \${segment}

        # Append "Coverage" to the second column header of the CSV
        #awk -F, 'NR == 1 {\$2 = \$2 " Coverage"; print; next} {print}' OFS=, \$output_csv  > temp.csv && mv temp.csv \$output_csv 
        # Append "Coverage" to the second column header of the CSV and rename the column
        awk -F, 'NR == 1 {
            split(\$2, a, "-"); 
            \$2 = "Coverage-" a[2]; 
            print; 
            next
        } {print}' OFS=, \$output_csv > temp.csv && mv temp.csv \$output_csv

        txt_filename="${meta.id}_\${segment}_coverage.txt"

        # Use awk to check for numbers above XX and capture any such number
        number_above_XX=\$(awk -F, '{for(i=1; i<=NF; i++) if(\$i+0 > ${seq_quality_thershold}) {print \$i; exit}}' "\$txt_filename")
    
        if [ ! -z "\$number_above_XX" ]; then
            echo "Found a number above XX: \$number_above_XX. Renaming \$fasta_file"
            # Rename the FASTA file to indicate it has passed
            mv "\$fasta_file" "\${fasta_file%.*}.fasta"
        else
            echo "No number above XX found in \$fasta_file"
            # Optionally, handle files that don't meet the criteria
        fi
    done



    cat *.fasta > "${meta.id}_coverage.fa"






    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

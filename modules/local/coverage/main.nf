
process COVERAGE {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
  

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(sequences), path(subtype)
    val(seq_quality_thershold)
   

    output:
    tuple val(meta), path("*.csv"), emit: coverage
    path("*.csv"), emit: coverage_report
    tuple val(meta), path("*fasta"), path(subtype), path("*csv") emit: filtered_fasta
    tuple val(meta), path("*coverage.fa"), path(subtype), emit:  merged_filtered_fasta
    path "versions.yml", emit: versions

    path("*fa"), emit: filtered_fasta_report

    when:
    task.ext.when == null || task.ext.when

    //errorStrategy 'ignore'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    for fasta_file in ${sequences}; do
        filename=\$(basename \$fasta_file)
        echo "Processing \$filename"  
        filename_no_ext=\${filename%.*}  
        segment_subtype=\${filename_no_ext#*-} 
        segment=\${segment_subtype%-*}  
        subtype_name=\${segment_subtype#*-} 
        
        echo "Segment: \$segment"

        output_csv="${meta.id}_\${segment}_coverage.csv"
        python /project-bin/coverage_finder.py "\$fasta_file" \$output_csv ${meta.id} \${segment}

        awk -F, 'NR == 1 { 
            sub(/.*-/, "", \$2); 
            \$2 = "Coverage-" \$2; 
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
        fi
    done

    cat *.fasta > "${meta.id}_coverage.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

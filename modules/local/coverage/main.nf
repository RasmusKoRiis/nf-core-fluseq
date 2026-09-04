
process COVERAGE {
    tag "$meta.id"
    label 'process_single'
  

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    tuple val(meta), path(sequences), path(subtype)
    val(seq_quality_thershold)
   

    output:
    tuple val(meta), path("coverage_reports/*.csv"), emit: coverage
    path("coverage_reports/*.csv"), emit: coverage_report
    tuple val(meta), path("passed/*.fasta"), path(subtype), path("coverage_reports/*.csv"), emit: filtered_fasta
    tuple val(meta), path("passed/${meta.id}_coverage.fa"), path(subtype), emit:  merged_filtered_fasta
    path "versions.yml", emit: versions

    path("passed/*.fasta"), emit: filtered_fasta_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail
    mkdir -p coverage_reports passed
    passed_files=()

    for fasta_file in ${sequences}; do
        filename=\$(basename \$fasta_file)
        echo "Processing \$filename"  
        filename_no_ext=\${filename%.*}  
        segment_subtype=\${filename_no_ext#*-} 
        segment=\${segment_subtype%-*}  
        subtype_name=\${segment_subtype#*-} 
        
        echo "Segment: \$segment"

        output_csv="${meta.id}_\${segment}_coverage.csv"
        coverage_finder.py "\$fasta_file" "\$output_csv" ${meta.id} "\${segment}"

        awk -F, 'NR == 1 {
            sub(/.*-/, "", \$2); 
            \$2 = "Coverage-" \$2; 
            print; 
            next 
        } {print}' OFS=, \$output_csv > temp.csv && mv temp.csv \$output_csv

        txt_filename="${meta.id}_\${segment}_coverage.txt"
        if [ ! -s "\$output_csv" ] || [ ! -s "\$txt_filename" ]; then
            echo "Coverage calculation did not produce expected output for \$fasta_file" >&2
            exit 1
        fi

        # Use awk to check for numbers above XX and capture any such number
        number_above_XX=\$(awk -F, '{for(i=1; i<=NF; i++) if(\$i+0 > ${seq_quality_thershold}) {print \$i; exit}}' "\$txt_filename")
    
        if [ ! -z "\$number_above_XX" ]; then
            echo "Found a number above XX: \$number_above_XX. Renaming \$fasta_file"
            # Copy only passing inputs into a dedicated output directory. This
            # prevents staged inputs and failed segments from being published.
            passed_fasta="passed/${meta.id}_\${segment}.fasta"
            cp "\$fasta_file" "\$passed_fasta"
            passed_files+=("\$passed_fasta")
        else
            echo "No number above XX found in \$fasta_file"
        fi
    done

    if [ "\${#passed_files[@]}" -eq 0 ]; then
        echo "No FASTA segments passed the coverage threshold for ${meta.id}" >&2
        exit 1
    fi
    cat "\${passed_files[@]}" > "passed/${meta.id}_coverage.fa"
    mv "\$output_csv" coverage_reports/ 2>/dev/null || true
    for report in ${meta.id}_*_coverage.csv; do
        [ -e "\$report" ] && mv "\$report" coverage_reports/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

process AMINOACIDTRANSLATION {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'ignore'

 
    container 'docker.io/rasmuskriis/nextclade-python@sha256:86ee1b9a00da7af2c113aaf937da3554cc72c1f954d79970027941eb2cf7ce52'


    input:
    tuple val(meta), path(fasta), path(subtype)
    path(dataset)
    

    output:
    path("*nextclade_mutations.csv"), emit: nextclade_csv
    tuple val(meta), path("*_nextclade_lookup_mutations.csv"), path(subtype), emit: mutation_lookup_csv
    tuple val(meta), path("*translation*fasta"), path(subtype), emit: aminoacid_sequence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def fasta_files = (fasta instanceof List ? fasta : [fasta]).collect { "'${it.toString().replace("'", "'\\''")}'" }.join(' ')
    """
    set -euo pipefail
    # COVERAGE renames files, but retains the segment and subtype in the header.
    subtype_name=\$(tr -d '\\r\\n' < "${subtype}")
    for fasta_file in ${fasta_files}; do
        header=\$(head -n 1 "\$fasta_file" | tr -d '\\r')
        segment_subtype=\${header#*|}
        segment=\${segment_subtype%-"\$subtype_name"}
        segment_name=\${segment##*-}
        if [[ "\$header" != '>'*'|'* || -z "\$subtype_name" || "\$segment" == "\$segment_subtype" ]]; then
            echo "Invalid FASTA header for ${meta.id} in \$fasta_file: \$header (expected subtype \$subtype_name)" >&2
            exit 1
        fi

        dataset_sample="${dataset}/\${subtype_name}_\${segment_name}"
        output_dir="${meta.id}_\${segment}_nextclade_output"
        
        nextclade run \
            --input-dataset "\$dataset_sample" \
            --output-all "\$output_dir" \
            "\$fasta_file"
        
        # Only move files, and keep outputs from different segments distinct.
        for output_file in "\$output_dir"/*; do
            [[ -f "\$output_file" ]] || continue
            output_name=\$(basename "\$output_file")
            mv "\$output_file" "${meta.id}_\${segment_name}_\${output_name}"
        done

        csv_conversion_nextclade.py \
                "${meta.id}_\${segment_name}_nextclade.csv" \
                "${meta.id}"
   
    done



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS

    """
}

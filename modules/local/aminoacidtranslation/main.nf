process AMINOACIDTRANSLATION {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'ignore'
 
    container 'docker.io/nextstrain/nextclade:latest'


    input:
    tuple val(meta), path(fasta), path(subtype)
    path(dataset)
    

    output:
    tuple val(meta), path("*.csv"), emit: nextclade_csv
    tuple val(meta), path("*translation*fasta"), path(subtype), emit: aminoacid_sequence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    for fasta_file in ${fasta}; do

        filename=\$(basename \$fasta_file)
        filename_no_ext=\${filename%.*}  
        segment_subtype=\${filename_no_ext#*-} 
        segment=\${segment_subtype%-*}  
        subtype_name=\${segment_subtype#*-} 
        
        dataset_sample=${dataset}/"\${subtype_name}_\${segment}"
        
        nextclade run \
            --input-dataset "\$dataset_sample" \
            --output-all=${meta.id}_\${segment}_nextclade_output/ \
            \$fasta_file
        
        for file in ${meta.id}_\${segment}_nextclade_output/*; do
            cat "\$file"
            basename=\$(basename \$file)
            mv "\$file" ./${meta.id}_\${basename}
        done

        subtype=\$(cat ${subtype})
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS

    """
}

process AMINOACIDTRANSLATION {
    tag "${meta.id}"
    label 'process_single'

 
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
    """
    set -euo pipefail
    for fasta_file in ${fasta}; do

    filename=\$(basename "\$fasta_file")
    filename_no_ext="\${filename%.*}"
    segment_subtype="\${filename_no_ext#*-}"       # INFA_01-HA-H5N1
    segment="\${segment_subtype%-*}"               # INFA_01-HA
    subtype_name="\${segment_subtype##*-}"         # H5N1
    segment_name="\${segment##*-}"                 # HA
    
        
        dataset_sample=${dataset}/"\${subtype_name}_\${segment_name}"
        
        nextclade run \
            --input-dataset "\$dataset_sample" \
            --output-all=${meta.id}_\${segment}_nextclade_output/ \
            \$fasta_file
        
        for file in ${meta.id}_\${segment}_nextclade_output/*; do
            cat "\$file"
            basename=\$(basename \$file)
            mv "\$file" ./${meta.id}_\${basename}
        done

        echo "1"

        mv ${meta.id}_nextclade.csv ${meta.id}_\${segment_name}_nextclade.csv

        subtype=\$(cat ${subtype})

        #NEXTCLADE CONVERSION
        echo ${fasta}
        csv_conversion_nextclade.py \
                ${meta.id}_\${segment_name}_nextclade.csv \
                ${meta.id} 
   
    done



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS

    """
}

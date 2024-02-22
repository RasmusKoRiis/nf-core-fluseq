process AMINOACIDTRANSLATION {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/nextstrain/nextclade:latest'
    //container logic as needed

    input:
    tuple val(meta), path(fasta), path(na_fasta)
    path(dataset)
    tuple val(meta), path(ha_fasta), path(subtype)

    output:
    tuple val(meta), path("*.csv"), emit: nextclade_csv
    tuple val(meta), path("*fasta"), emit: aminoacid_sequence
    tuple val(meta), path("${meta.id}_nextclade_output"), emit: nextclade_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """

    dataset=${dataset}/h5_ha
    nextclade run \
        --input-dataset "\$dataset" \
        --output-all=${meta.id}_nextclade_output/ \
        $fasta
    
    mv ${meta.id}_nextclade_output/nextclade.csv ${meta.id}_nextclade.csv

    for file in ${meta.id}_nextclade_output/nextclade.cds_translation.*.fasta;
    do
        filename=\$(basename -- "\$file")
        new_filename="${meta.id}_\${filename#nextclade.cds_translation.}"
        mv "\$file" "\${new_filename}"
    done

    subtype=\$(cat ${subtype})


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS

    """
}

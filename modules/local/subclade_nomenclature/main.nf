process SUBCLADE_NOMENCLATURE {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'
    input:
    tuple val(meta), path(fasta), path(subtype), path(coverage_csv)
    path rules_dir
    path caller_script
    path characterisation_script
    path characterisation_guidelines

    output:
    tuple val(meta), path("${meta.id}_subclade_nomenclature.csv"), emit: calls
    path("${meta.id}_subclade_nomenclature.csv"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${caller_script} \\
        --sample-id ${meta.id} \\
        --subtype-file ${subtype} \\
        --rules-dir ${rules_dir} \\
        --output ${meta.id}_subclade_nomenclature_raw.csv \\
        ${fasta}

    python ${characterisation_script} \\
        --input ${meta.id}_subclade_nomenclature_raw.csv \\
        --subtype-file ${subtype} \\
        --guidelines-dir ${characterisation_guidelines} \\
        --output ${meta.id}_subclade_nomenclature.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

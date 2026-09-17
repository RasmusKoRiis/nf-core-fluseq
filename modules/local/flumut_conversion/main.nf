process FLUMUT_CONVERSION {
    label 'process_single'

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    tuple val(meta), path(tsv)

    output:
    path("${meta.id}_flumut_report.csv"), emit: flumut_report
    path "versions.yml", emit: versions
   
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """  
    set -euo pipefail
    flumut_conversion.py $tsv $meta.id

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """

}

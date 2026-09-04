process DRUG_RESISTANCE_REPORT {
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path subtype
    path resistance
    path id_map
    val  runid

    output:
    path("${runid}_drug_resistance_report.csv"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    drug_resistance_report.py \
        --id-map ${id_map} \
        --subtype ${subtype} \
        --resistance ${resistance} \
        --output ${runid}_drug_resistance_report.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

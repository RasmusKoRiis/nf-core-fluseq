process DRUG_RESISTANCE_REPORT {
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    path subtype
    path resistance
    path id_map
    val  runid

    output:
    path("${runid}_drug_resistance_report.csv"), emit: report

    script:
    """
    python /project-bin/drug_resistance_report.py \
        --id-map ${id_map} \
        --subtype ${subtype} \
        --resistance ${resistance} \
        --output ${runid}_drug_resistance_report.csv
    """
}

process REPORT_QC_HTML {
    tag "${report_csv.simpleName}"
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path report_csv
    path report_script
    path report_template

    output:
    path '*_qc.html', emit: html
    path '*_qc.json', emit: summary
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python "${report_script}" "${report_csv}" --template "${report_template}"

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        report_qc_html: \$(python "${report_script}" --version)
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}

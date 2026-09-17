process SLIM_GENIN2_REPORT {
    tag { genin2_csv.simpleName }
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path genin2_csv

    output:
    path "${genin2_csv.simpleName}_trimmed.csv", emit: genin2_report_trimmed
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    slim_genin2_report.py ${genin2_csv} ${genin2_csv.simpleName}_trimmed.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}

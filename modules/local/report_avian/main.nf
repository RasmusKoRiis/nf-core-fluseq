process REPORT_AVIAN {

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    /*
     * Five lists of CSV paths + run ID.
     * Each list arrives as a Bash array; we don’t have to touch them.
     */
    input:
    path subtype
    path genotype
    path coverage
    path mammalian_mutations
    path nextclade
    val  runid

    /*
     * One merged file + (optionally) keep all csvs for provenance.
     */
    output:
    path "fluseq_merged_report.csv", emit: report
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    # Merge every CSV in the current directory
    reportavian.py

    # Rename to a fixed, pipeline-wide name
    mv merged_report.csv fluseq_merged_report.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

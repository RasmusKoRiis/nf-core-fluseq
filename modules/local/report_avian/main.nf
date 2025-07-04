process REPORT_AVIAN {

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions "-v ${baseDir}/bin:/project-bin"   // reportavian.py lives here

    /*
     * Six lists of CSV paths – exactly what you pass with .collect()
     * Each list arrives as a Bash array; we don’t have to touch them.
     */
    input:
    path subtype
    path genotype
    path coverage
    path mammalian_mutations
    path lookup_inhibition
    path lookup_mammalian

    /*
     * One merged file + (optionally) keep all csvs for provenance.
     */
    output:
    path "fluseq_merged_report.csv", emit: report
    path "*.csv",                    emit: all   // keep originals + merged

    script:
    """
    # Merge every CSV in the current directory
    python /project-bin/reportavian.py

    # Rename to a fixed, pipeline-wide name
    mv merged_report.csv fluseq_merged_report.csv
    """
}

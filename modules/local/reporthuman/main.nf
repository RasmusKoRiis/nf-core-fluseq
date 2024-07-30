process REPORTHUMAN {
    label 'process_single'
    debug = true


    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    path(subtype)
    path(coverage)
    path(mutation_human)
    path(mutation_inhibtion)
    path(lookup)
    path(nextclade_summary_ha)
    path(nextclade_sample)
    path(mutation_vaccine)
    val runid
    path(filtered_fasta)
    


    output:

    path("${runid}.csv"), emit: report
    path("${runid}.fasta"), emit: filtered_fasta


    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    python /project-bin/report.py

    #Add constant parameters to the report
    awk -v runid=${runid} -v OFS=',' '{ if (NR == 1) { print \$0, "RunID" } else { print \$0, runid } }' merged_report.csv > ${runid}.csv

    #Merge all filtered fasta files
    cat ${filtered_fasta} > ${runid}.fasta
    


    """

}

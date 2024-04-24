process REPORT {
    label 'process_single'
    debug = true


    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    path(subtype)
    path(genotype)
    path(coverage)
    path(mutation_mamilian)
    path(mutation_inhibtion)
    path(lookup)

    output:

    path("*_report.csv"), emit: report
    path("*.csv"), emit: all


    when:
    task.ext.when == null || task.ext.when

    //errorStrategy 'ignore'

    script:
    """  
    python /project-bin/report.py 
    """

}

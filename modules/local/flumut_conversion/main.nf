process FLUMUT_CONVERSION {
    label 'process_single'
    errorStrategy 'ignore'

    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    tuple val(meta), path(tsv)

    output:
    path("*_report.csv"), emit: flumut_report
   
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """  
    python /project-bin/flumut_conversion.py $tsv $meta.id
    """

}

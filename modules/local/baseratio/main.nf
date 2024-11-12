process BASERATIO {
    label 'process_single'



    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    path(depth)

    
    output:

    path("*.csv"), emit: report



    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    python /project-bin/depth_analysis_merge.py 

    
    """

}

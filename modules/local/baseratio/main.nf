process BASERATIO {
    label 'process_single'
   



    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path(depth)

    
    output:

    path("depth_report.csv"), emit: report
    path "versions.yml", emit: versions



    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    set -euo pipefail
    depth_analysis_merge.py

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS

    """

}

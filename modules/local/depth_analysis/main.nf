
process DEPTH_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    debug false
    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //conda "bioconda::irma=1.0.3"
    //container 'docker.io/rasmuskriis/cdc_irma_custom:1.0'
    container 'docker.io/rasmuskriis/nextclade-python@sha256:86ee1b9a00da7af2c113aaf937da3554cc72c1f954d79970027941eb2cf7ce52'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0':
        //'biocontainers/irma:1.0.3--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
   

    output:
    path("${meta.id}_long.csv") , emit: depth_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    depth_analysis.py \
                ${meta.id}   

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

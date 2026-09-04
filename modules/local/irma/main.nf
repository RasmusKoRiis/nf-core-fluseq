
process IRMA {
    tag "$meta.id"
    label 'process_medium'
   


    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //conda "bioconda::irma=1.0.3"
    //container 'docker.io/rasmuskriis/cdc_irma_custom:1.0'
    // v1.3.4+ fixes FASTQ deduplication/inflation for tab-delimited ONT headers.
    container 'docker.io/cdcgov/irma:v1.3.5'
    // The hardened image defaults to UID 65532, which cannot write host-owned work directories.
    containerOptions = '-u $(id -u):$(id -g)'


    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/irma:1.0.3--pl5321hdfd78af_0':
        //'biocontainers/irma:1.0.3--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
   

    output:
    tuple val(meta), path("$meta.id/*.fasta") , emit: fasta
    tuple val(meta), path("$meta.id/*.bam"),  path("$meta.id/*.bai"), emit: bam
    tuple val(meta), path("$meta.id/*.vcf") , emit: vcf
    tuple val(meta), path("$meta.id/tables/READ_COUNTS.txt") , emit: read_count
    tuple val(meta), path("$meta.id/figures/*.pdf") , emit: figures
    tuple val(meta), path("$meta.id/amended_consensus/*.fa") , emit: amended_consensus
    tuple val(meta), path("$meta.id/secondary") , emit: secondary
    tuple val(meta), path("$meta.id/tables/*txt") , emit: alleles
    path "versions.yml", emit: versions
  
  

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    IRMA FLU-minion $fastq ${meta.id}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        irma: \$(IRMA --version 2>&1 | head -n 1 || true)
    END_VERSIONS

    """
}

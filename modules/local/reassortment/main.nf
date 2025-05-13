process REASSORTMENT {
    tag "$meta.id"
    label 'process_single'
    //errorStrategy 'ignore'

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(sequences)                         // Multi-segment FASTA with headers like >ID_segment
    path(reassortment_database)                              // Multi-lineage reference DB

    output:
    tuple val(meta), path("${meta.id}_reassortment_summary.csv"), emit: genotype_report
    path("versions.yml"), emit: versions

    script:
    def prefix = meta.id

    """
    blastn -query ${sequences} -subject ${reassortment_database} -outfmt 6 -max_target_seqs 5 -out ${prefix}_blast.tsv

    python /project-bin/detect_reassortment.py \\
        --blast ${prefix}_blast.tsv \\
        --output ${prefix}_reassortment_summary.csv \\
        --sample ${prefix}

    echo "${task.process}:" > versions.yml
    echo "    blast: \$(blastn -version 2>&1 | head -n 1)" >> versions.yml
    """
}

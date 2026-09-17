process REASSORTMENT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    stageInMode 'copy'          
    input:
    tuple val(meta), path(sequences)                         // Multi-segment FASTA with headers like >ID_segment
    file reassortment_database                           // Multi-lineage reference DB

    output:
    tuple val(meta), path("${meta.id}_reassortment_summary.csv"), emit: genotype
    path("${meta.id}_reassortment_summary.csv"), emit: genotype_report
    path("versions.yml"), emit: versions

    script:
    def prefix = meta.id

    """
    set -euo pipefail
    blastn -query ${sequences} -subject ${reassortment_database} -outfmt 6 -max_target_seqs 5 -out ${prefix}_blast.tsv

    detect_reassortment.py \\
        --blast ${prefix}_blast.tsv \\
        --output ${prefix}_reassortment_summary.csv \\
        --sample ${prefix} \\
        --metadata "\$(dirname "\$(command -v detect_reassortment.py)")/reassortment_reference_metadata.csv"

    echo "${task.process}:" > versions.yml
    echo "    blast: \$(blastn -version 2>&1 | head -n 1)" >> versions.yml
    """
}

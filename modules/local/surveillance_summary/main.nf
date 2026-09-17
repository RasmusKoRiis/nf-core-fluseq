process SURVEILLANCE_SUMMARY {
    tag "${mode}"
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    val mode
    path fasta
    path coverage
    path subtype
    path subtype_hits
    path reassortment
    path resistance
    path resistance_databases
    path summary_script

    output:
    path "segment_qc.tsv", emit: segment_qc
    path "reassortment_summary.tsv", emit: reassortment_summary
    path "resistance_summary.tsv", emit: resistance_summary
    path "sample_summary.json", emit: sample_summary
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    python ${summary_script} \
        --mode '${mode}' \
        --fasta ${fasta} \
        --coverage ${coverage} \
        --subtype ${subtype} \
        --subtype-hits ${subtype_hits} \
        --reassortment ${reassortment} \
        --resistance ${resistance} \
        --resistance-databases ${resistance_databases}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        surveillance_summary: "1.0.0"
        python: \$(python --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}

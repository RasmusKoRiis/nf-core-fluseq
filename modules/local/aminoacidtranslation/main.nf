process AMINOACIDTRANSLATION {
    tag "${meta.id}"
    label 'process_single'

    //container "YOUR-TOOL-HERE"
    //container logic as needed

    input:
    tuple val(meta), path(fasta), path(ha_fasta)

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    Echo "hello"

    """
}

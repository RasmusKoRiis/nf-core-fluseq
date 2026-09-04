#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = null
params.samples_dir = null

include { INPUT_CHECK } from '../../subworkflows/local/input_check'

workflow {
    if (!params.input) {
        error 'Provide --input to the input-check harness.'
    }

    INPUT_CHECK(Channel.fromPath(params.input, checkIfExists: true))

    INPUT_CHECK.out.reads
        .map { meta, reads ->
            assert meta.id == 'SAMPLE01'
            assert reads.size() == 1
            assert reads[0].exists()
            "${meta.id}\t${reads[0]}"
        }
        .collectFile(name: 'validated_input.tsv', newLine: true)
        .view { path -> "Validated input contract: ${path}" }
}

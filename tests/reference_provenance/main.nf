#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.references = null

include { REFERENCE_PROVENANCE } from '../../modules/local/reference_provenance/main'

workflow {
    if (!params.references) {
        error 'Provide --references to the reference-provenance harness.'
    }

    REFERENCE_PROVENANCE(Channel.value([file(params.references, checkIfExists: true)]))

    REFERENCE_PROVENANCE.out.manifest
        .map { manifest ->
            def lines = manifest.readLines()
            assert lines[0] == 'path\tsha256'
            assert lines.size() == 9
            assert lines.drop(1).every { it ==~ /.+\t[0-9a-f]{64}/ }
            "Validated reference manifest: ${manifest}"
        }
        .view()
}

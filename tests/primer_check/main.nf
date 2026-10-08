include { PRIMER_CHECK_RUN } from '../../subworkflows/local/primer_check/main'

process AFTER_PRIMER_CHECK {
    publishDir params.outdir, mode: 'copy'
    input:
    val reports
    output:
    path 'continued.txt'
    script:
    """
    echo 'Pipeline continued after optional primer checks' > continued.txt
    """
}

workflow {
    def manifest = new groovy.json.JsonSlurper().parse(file(params.test_manifest).toFile())
    samples = Channel.fromList(manifest).map { sample ->
        tuple([id: sample.sample_id, subtype: sample.subtype ?: ''],
              sample.fasta.collect { file(it, checkIfExists: true) },
              sample.subtype_file ? file(sample.subtype_file, checkIfExists: true) : [])
    }
    PRIMER_CHECK_RUN(Channel.value(file(params.test_fasta, checkIfExists: true)), samples,
        [virus: params.test_virus, run_id: 'SYNTHETIC',
         assays: params.test_virus == 'Influenza' ? ['pcr'] : ['pcr', 'ngs'],
         ngs_dir: params.primer_check_ngs_dir, ngs_scheme: 'TEST', offline: true])
    AFTER_PRIMER_CHECK(PRIMER_CHECK_RUN.out.csv.toList())
}

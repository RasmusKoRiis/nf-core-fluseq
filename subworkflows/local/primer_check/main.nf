// Maintained in primer-checker/integrations/nextflow; use sync_nextflow_modules.py.
include { PRIMER_CHECK } from '../../../modules/local/primer_check/main'

workflow PRIMER_CHECK_RUN {
    take:
    samples // tuple(meta, consensus file/list, subtype file or [])
    settings // virus, run_id, ngs_dir, ngs_scheme, offline, assays

    main:
    // Missing optional QC resources must fail inside the ignored process, not
    // during workflow initialisation or channel construction.
    def existing = { value ->
        if (!value) return []
        def resolved = file(value)
        resolved.exists() ? resolved : []
    }
    def pcr = existing(params.primer_check_pcr)
    def ngs = existing(settings.ngs_dir)
    def source = existing(params.primer_check_source)

    ch_inputs = samples.toList().map { entries ->
        def records = []
        def sequences = []
        def subtypes = []
        entries.sort { a, b -> a[0].id <=> b[0].id }.each { meta, fastas, subtype ->
            def files = fastas instanceof List ? fastas : [fastas]
            def indices = files.collect { fasta ->
                sequences.add(fasta)
                sequences.size() - 1
            }
            def subtypeIndex = null
            if (subtype) {
                subtypes.add(subtype)
                subtypeIndex = subtypes.size() - 1
            }
            records.add([id: meta.id, subtype: meta.subtype ?: '',
                         fasta_indices: indices, subtype_index: subtypeIndex])
        }
        tuple(records, sequences, subtypes)
    }.combine(Channel.fromList(settings.assays))

    PRIMER_CHECK(ch_inputs, pcr, ngs, source, settings)

    // Also emit a fresh status when an ignored failure produces no reports.
    // On -resume, older published reports may still exist in the output folder.
    ch_task_status = PRIMER_CHECK.out.csv.toList().map { reports ->
        def completed = reports.collect { it.getFileName().toString() }
        (['Assay_Type,Status'] + settings.assays.collect { assay ->
            "${assay},${completed.contains(assay + '_primer_report.csv') ? 'complete' : 'failed'}"
        }).join('\n') + '\n'
    }.collectFile(name: 'task_status.csv', storeDir: "${params.outdir}/primer_check", cache: false)

    emit:
    csv = PRIMER_CHECK.out.csv
    html = PRIMER_CHECK.out.html
    status = PRIMER_CHECK.out.status
    provenance = PRIMER_CHECK.out.provenance
    task_status = ch_task_status
}

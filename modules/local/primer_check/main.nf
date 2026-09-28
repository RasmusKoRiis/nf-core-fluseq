// Maintained in primer-checker/integrations/nextflow; use sync_nextflow_modules.py.
process PRIMER_CHECK {
    tag "${settings.virus}:${assay}"
    label 'process_single'
    errorStrategy 'ignore'
    cache false
    container params.primer_check_container
    containerOptions { workflow.containerEngine == 'docker' && !settings.offline && !checker_source ? '--pull=always' : '' }

    input:
    tuple val(records), path(consensus, stageAs: 'consensus??/*'), path(subtypes, stageAs: 'subtypes??/*'), val(assay)
    path pcr_database, stageAs: 'pcr_database'
    path ngs_directory, stageAs: 'ngs_primers'
    path checker_source, stageAs: 'checker_source'
    val settings

    output:
    path "${assay}_primer_report.csv", emit: csv
    path "${assay}_primer_report.html", emit: html
    path "${assay}_primer_report.status.csv", emit: status
    path "${assay}_primer_report.provenance.json", emit: provenance

    script:
    def asList = { value -> value instanceof List ? value : [value] }
    def sequences = asList(consensus)
    def types = asList(subtypes)
    def manifest = records.collect { record ->
        [sample_id: record.id, subtype: record.subtype,
         fasta: record.fasta_indices.collect { sequences[it].toString() },
         subtype_file: record.subtype_index == null ? null : types[record.subtype_index].toString()]
    }
    def encoded = groovy.json.JsonOutput.toJson(manifest).bytes.encodeBase64().toString()
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    def source = checker_source ? checker_source.toString() : '/opt/primer-checker'
    def arguments = ['--manifest', 'input_manifest.json', '--virus', settings.virus,
                     '--assay-type', assay, '--run-id', settings.run_id,
                     '--ngs-scheme', settings.ngs_scheme ?: '',
                     '--output-prefix', "${assay}_primer_report"]
    if (pcr_database) arguments.addAll(['--pcr-db', pcr_database.toString()])
    if (ngs_directory) arguments.addAll(['--ngs-dir', ngs_directory.toString()])
    def command = arguments.collect(quote).join(' ')
    """
    python3 -c 'import base64; from pathlib import Path; Path("input_manifest.json").write_bytes(base64.b64decode("${encoded}"))'
    python3 ${quote(source + '/scripts/run_pipeline_primer_check.py')} ${command}
    """

    stub:
    """
    printf 'Primer_Name,Hit_Status\n' > ${assay}_primer_report.csv
    printf '<html><body>Primer check stub</body></html>\n' > ${assay}_primer_report.html
    printf 'Sample_ID,Status\n' > ${assay}_primer_report.status.csv
    printf '{}\n' > ${assay}_primer_report.provenance.json
    """
}

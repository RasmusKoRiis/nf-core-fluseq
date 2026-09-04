// Validate either the routine barcode sheet or an explicit FASTQ sheet and
// return one tuple per sample: [meta, [fastq files]].

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    reads = samplesheet
        .splitCsv(header: true, sep: ',', strip: true)
        .toList()
        .flatMap { rows -> create_fastq_channels(rows, params.samples_dir) }

    emit:
    reads
}

def create_fastq_channels(List rows, samplesDir) {
    if (!rows) {
        error 'The input sample sheet contains no samples.'
    }

    def parsed = []
    rows.eachWithIndex { row, index ->
        def sample = (row.SequenceID ?: row.sample_id ?: row.sample)?.toString()?.trim()
        if (!sample) {
            error "Missing sample identifier on sample-sheet line ${index + 2}. Use SequenceID, sample_id, or sample."
        }
        sample = sample.replaceAll(/\s+/, '_')

        def reads = []
        def barcode = (row.Barcode ?: row.barcode)?.toString()?.trim()
        if (barcode) {
            if (!samplesDir) {
                error "Sample ${sample} uses a barcode but --samples_dir was not provided."
            }
            ['*.fastq.gz', '*.fq.gz'].each { suffix ->
                reads.addAll(files("${samplesDir}/${barcode}/${suffix}").findAll { it.exists() })
            }
        } else {
            [row.fastq_1, row.fastq_2].findAll { it?.toString()?.trim() }.each { readPath ->
                def value = readPath.toString().trim()
                def resolved = value.startsWith('/') || value ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/
                    ? file(value)
                    : samplesDir
                        ? file("${samplesDir}/${value}")
                        : file(value)
                if (!resolved.exists()) {
                    error "FASTQ file for sample ${sample} does not exist: ${resolved}"
                }
                reads.add(resolved)
            }
        }

        reads = reads.unique().sort { it.toString() }
        if (!reads) {
            def location = barcode ? "${samplesDir}/${barcode}" : 'the explicit FASTQ columns'
            error "No .fastq.gz or .fq.gz files were found for sample ${sample} in ${location}."
        }

        parsed.add(tuple([id: sample, single_end: reads.size() == 1], reads))
    }

    def duplicates = parsed.groupBy { it[0].id }.findAll { id, values -> values.size() > 1 }.keySet()
    if (duplicates) {
        error "Sample identifiers must be unique. Duplicates: ${duplicates.sort().join(', ')}"
    }
    parsed
}

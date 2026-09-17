//
// This file holds several functions specific to the main.nf workflow in the nf-core/fluseq pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    private static final List<String> MODES = [
        'human-fastq',
        'human-fasta',
        'avian-fastq',
        'avian-fasta'
    ]

    // Preserve the command-line contract used by the operational wrappers while
    // exposing correctly spelled parameter names to new callers.
    public static void applyCompatibilityAliases(params, log) {
        applyAlias(params, log, 'samples_dir', 'samplesDir')
        applyAlias(params, log, 'seq_quality_threshold', 'seq_quality_thershold')
        applyAlias(params, log, 'mammalian_mutation_db', 'mamalian_mutation_db')
        applyAlias(params, log, 'inhibition_mutation_db', 'inhibtion_mutation_db')
    }

    private static void applyAlias(params, log, String canonical, String legacy) {
        if (params[legacy] != null && params[legacy].toString() != '') {
            log.warn "Parameter --${legacy} is deprecated; use --${canonical}. The legacy name remains supported for routine wrappers."
        }
    }

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        if (!MODES.contains(params.file?.toString())) {
            Nextflow.error("Invalid --file mode '${params.file}'. Choose one of: ${MODES.join(', ')}")
        }

        def required = ['ha_database', 'na_database', 'sequence_references', 'nextclade_dataset']

        if (params.file.endsWith('-fastq')) {
            required.addAll(['input', 'samples_dir'])
        } else {
            required.add('fasta')
        }

        if (params.file.startsWith('avian')) {
            required.addAll([
                'genotype_database',
                'reassortment_database',
                'mammalian_mutation_db',
                'inhibition_mutation_db'
            ])
        } else {
            required.add('inhibition_mutation_db')
            def resistanceOnly = params.file == 'human-fasta' && params.drug_resistance_only.toString().toBoolean()
            if (!resistanceOnly) {
                required.add('reassortment_database')
            }
            if (params.file == 'human-fasta' && !resistanceOnly) {
                required.add('genotype_database')
            }
        }

        def missing = required.unique().findAll { key ->
            def value = resolvedParameter(params, key)
            value == null || value.toString().trim() == ''
        }
        if (missing) {
            Nextflow.error("Missing required parameter(s) for --file ${params.file}: ${missing.collect { '--' + it }.join(', ')}")
        }
    }

    private static Object resolvedParameter(params, String canonical) {
        def aliases = [
            samples_dir: 'samplesDir',
            seq_quality_threshold: 'seq_quality_thershold',
            mammalian_mutation_db: 'mamalian_mutation_db',
            inhibition_mutation_db: 'inhibtion_mutation_db'
        ]
        def legacy = aliases[canonical]
        if (legacy && params[legacy] != null && params[legacy].toString().trim() != '') {
            return params[legacy]
        }
        params[canonical]
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    
}

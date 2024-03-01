/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowFluseq.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { IRMA                        } from '../modules/local/irma/main'
include { AMINOACIDTRANSLATION        } from '../modules/local/aminoacidtranslation/main'
include { SUBTYPEFINDER               } from '../modules/local/blastn/main'
include { GENOTYPING                  } from '../modules/local/genotyping/main'
include { FASTA_CONFIGURATION         } from '../modules/local/seqkit/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


// Function to parse the sample sheet
def parseSampleSheet(sampleSheetPath) {
        return Channel
            .fromPath(sampleSheetPath)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row ->
                def sampleId = row.sample_id
                def files = file("${params.samplesDir}/${row.barcode}/*.fastq.gz")
                // Creating a metadata map
                def meta = [ id: sampleId, single_end: true ]
                return tuple(meta, files)
            }
}

workflow AVIAN {

    //
    // INPUT PARSE
    //

    ch_sample_information = parseSampleSheet(params.input)
    
    ch_versions = Channel.empty()

    ch_sample_information
        .map { meta, files ->
            tuple(meta, files.toList())
        }
    .set { read_input }

    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    //INPUT_CHECK (
    //    file(params.input)
    //)
    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: CAT_FASTQ
    //

    CAT_FASTQ (
        read_input
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    
    //
    // MODULE: IRMA
    //

    IRMA (
        CAT_FASTQ.out.reads
    )

    /// SUBTYPE CHANNEL
    IRMA.out.fasta
    .map { meta, files -> 
        def ha_files = files.findAll { it.getName().contains('_HA') }
        def na_files = files.findAll { it.getName().contains('_NA') }
        return (ha_files && na_files) ? tuple(meta, ha_files, na_files) : null
    }
    .filter { it != null }
    .set { IRMA_ha_na_fasta }

    /// GENOTYPING CHANNEL
    IRMA.out.amended_consensus
    .map { meta, files -> 
        def amended_consensus_files = files.findAll { it.getName().contains('.fa') }
        return (amended_consensus_files) ? tuple(meta, amended_consensus_files) : null
    }
    .filter { it != null }
    .set { IRMA_amended_consensus_files }

    //
    // MODULE: SUBTYPE FINDER
    //
    def currentDir = System.getProperty('user.dir')
    def fullPathHA = "${currentDir}/${params.ha_database}"
    def fullPathNA = "${currentDir}/${params.na_database}"


    SUBTYPEFINDER (
        IRMA_ha_na_fasta, fullPathHA, fullPathNA
    )


    /// MUTATION CHANNELS

    IRMA.out.amended_consensus
    .join(SUBTYPEFINDER.out.subtype, by: [0]) // Assuming meta.id is the first element in the tuple
    .map { items ->
        def meta = items[0] // The common meta.id
        def fasta = items[1] // The fasta file from IRMA_out_amended_consensus
        def subtype = items[2] // The subtype file from SUBTYPEFINDER_out_subtype
        return [meta, fasta, subtype]
    }
    .set { fasta_subtype }

    fasta_subtype
    .map { meta, files, subtype -> 
        def pb2 = files.findAll { it.getName().contains('_1') }
        return (pb2) ? tuple(meta, pb2, subtype) : null
    }
    .filter { it != null }
    .set { pb2_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def pb1 = files.findAll { it.getName().contains('_2') }
        return (pb1) ? tuple(meta, pb1, subtype) : null
    }
    .filter { it != null }
    .set { pb1_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def pa = files.findAll { it.getName().contains('_3') }
        return (pa) ? tuple(meta, pa, subtype) : null
    }
    .filter { it != null }
    .set { pa_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def ha = files.findAll { it.getName().contains('_4') }
        return (ha) ? tuple(meta, ha, subtype) : null
    }
    .filter { it != null }
    .set { ha_fasta }

    fasta_subtype
    .map { meta, files, subtype -> 
        def np = files.findAll { it.getName().contains('_5') }
        return (np) ? tuple(meta, np, subtype) : null
    }
    .filter { it != null }
    .set { np_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def na = files.findAll { it.getName().contains('_6') }
        return (na) ? tuple(meta, na, subtype) : null
    }
    .filter { it != null }
    .set { na_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def m = files.findAll { it.getName().contains('_7') }
        return (m) ? tuple(meta, m, subtype) : null
    }
    .filter { it != null }
    .set { m_fasta }
    
    fasta_subtype
    .map { meta, files, subtype -> 
        def ns = files.findAll { it.getName().contains('_8') }
        return (ns) ? tuple(meta, ns, subtype) : null
    }
    .filter { it != null }
    .set { ns_fasta }



    //
    // MODULE: GENOTYPING
    //

    ch_genotype_database = params.genotype_database

    GENOTYPING (
        IRMA_amended_consensus_files, ch_genotype_database
    )


    //
    // MODULE: FASTA CONFIGURATION
    //

    
    FASTA_CONFIGURATION(
         fasta_subtype 
    )

    // FASTA_CONFIGURATION.out.fasta.view()



    //
    // MODULE: AMINO ACID TRANSLATION
    //

    
    def fullPath_nextclade_dataset = "${currentDir}/${params.nextclade_dataset}"

    AMINOACIDTRANSLATION (
        IRMA_ha_na_fasta, fullPath_nextclade_dataset, fasta_subtype
    )

    //AMINOACIDTRANSLATION.out.fasta.view()


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        CAT_FASTQ.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: REPORT
    //






    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFluseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowFluseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

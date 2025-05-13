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
include { NEXTCLADE                   } from '../modules/local/nextclade/main'
include { SUBTYPEFINDER               } from '../modules/local/blastn/main'
include { GENOTYPING                  } from '../modules/local/genotyping/main'
include { COVERAGE                    } from '../modules/local/coverage/main'
include { FASTA_CONFIGURATION         } from '../modules/local/seqkit/main'
include { MUTATIONHUMAN               } from '../modules/local/mutationhuman/main'
include { TABLELOOKUP                 } from '../modules/local/tablelookup/main'
include { REPORTHUMAN                 } from '../modules/local/reporthuman/main'
include { TECHNICAL                   } from '../modules/local/technical/main'
include { DEPTH_ANALYSIS              } from '../modules/local/depth_analysis/main'
include { BASERATIO                   } from '../modules/local/baseratio/main'
include { CHOPPER                     } from '../modules/local/chopper/main'
include { REASSORTMENT                } from '../modules/local/reassortment/main'






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
                def sampleId = row.SequenceID
                def files = file("${params.samplesDir}/${row.Barcode}/*.fastq.gz")
                // Creating a metadata map
                def meta = [ id: sampleId, single_end: true ]
                return tuple(meta, files)
            }
}

workflow HUMAN {

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
    // MODULE: CAT_FASTQ
    //Concatenate the fastq files into singel files for each sample as required by IRMA

    CAT_FASTQ (
        read_input
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    
    //
    // MODULE: Chopper
    //CQuality and length trim of FASTQ-sequences

    CHOPPER (
        CAT_FASTQ.out.reads
    )
   

    //
    // MODULE: IRMA
    //Run IRMA to get the consensus sequences

    IRMA (
        CHOPPER.out.chopperfastq
        
    )
    
    

    //
    // MODULE: DEPTH ANALYSIS
    //

   DEPTH_ANALYSIS (
     IRMA.out.bam
    )


    //
    // MODULE: MDEPTH REPORT
    //

    BASERATIO (
        DEPTH_ANALYSIS.out.depth_report.collect()
    )

    //ch_versions = ch_versions.mix(BASERATIO.out.versions.first())
    
    
    //
    // MODULE: IRMA STAT
    //Get technical data from IRMA

    TECHNICAL (
        IRMA.out.read_count
    )


    /// SUBTYPE CHANNEL
    /// This channel is used to store the fasta files for each segment
    IRMA.out.amended_consensus
    .map { meta, files -> 
        def ha_files = files.findAll { it.getName().contains('_4') }
        def na_files = files.findAll { it.getName().contains('_6') }
        def pa_files = files.findAll { it.getName().contains('_3') }
        def pb1_files = files.findAll { it.getName().contains('_2') }
        def pb2_files = files.findAll { it.getName().contains('_1') }
        def ns_files = files.findAll { it.getName().contains('_7') }
        def np_files = files.findAll { it.getName().contains('_5') }
        def m_files = files.findAll { it.getName().contains('_8') }
        return (ha_files && na_files && pa_files && pb1_files && pb2_files && ns_files && np_files && m_files) ? tuple(meta, ha_files, na_files, pa_files, pb1_files, pb2_files, ns_files, np_files, m_files) : null
    }
    .filter { it != null }
    .set { IRMA_all_fasta }

        IRMA.out.amended_consensus
    .map { meta, files -> 
        def ha_files = files.findAll { it.getName().contains('_4') }
        def na_files = files.findAll { it.getName().contains('_6') }
        return (ha_files && na_files) ? tuple(meta, ha_files, na_files) : null
    }
    .filter { it != null }
    .set { IRMA_ha_na_fasta }


    //
    // MODULE: SUBTYPE FINDER
    //Use BLAST on local database to find the subtype of the sequences

    SUBTYPEFINDER (
        IRMA_ha_na_fasta, Channel.value(file(params.ha_database)),Channel.value(file(params.na_database))
    )

    ch_versions = ch_versions.mix(SUBTYPEFINDER.out.versions.first())


    /// MUTATION CHANNELS
    /// This channel is used to store the fasta files for each segment togheter with subtype and ID

    IRMA.out.amended_consensus
    .join(SUBTYPEFINDER.out.subtype, by: [0]) // Assuming meta.id is the first element in the tuple
    .map { items ->
        def meta = items[0] // The common meta.id
        def fasta = items[1] // The fasta file from IRMA_out_amended_consensus
        def subtype = items[2] // The subtype file from SUBTYPEFINDER_out_subtype
        return [meta, fasta, subtype]
    }
    .set { fasta_subtype }


    //
    // MODULE: FASTA CONFIGURATION
    //Configure the fasta file names for reporting

    
    FASTA_CONFIGURATION (
         fasta_subtype 
    )

    //
    // MODULE: REASSORTMENT
    //Configure the fasta file names for reporting

    
    REASSORTMENT (
         FASTA_CONFIGURATION.out.fasta_flumut, Channel.value(file(params.reassortment_database))
    )


    //
    // MODULE: COVERAGE
    //Calculate the coverage of the sequences and filter out low quality sequences

    /// Coverage threshold from the params/user
    def seq_quality_thershold = params.seq_quality_thershold

    COVERAGE (
         FASTA_CONFIGURATION.out.fasta, seq_quality_thershold
    )

    ch_versions = ch_versions.mix(COVERAGE.out.versions.first())


    //
    // MODULE: NEXTCLADE
    //Translate the nucleotide sequences to amino acid sequences using Nextclade

    NEXTCLADE (
        COVERAGE.out.filtered_fasta
    )

    ch_versions = ch_versions.mix(NEXTCLADE.out.versions.first())

    //
    // MODULE: MUTATION
    //Compare translated amino acid sequences to references to find mutations

    //def fullPath_references_2 = "${currentDir}/${params.sequence_references}"
    
    MUTATIONHUMAN  (
       NEXTCLADE.out.aminoacid_sequence, Channel.value(file(params.sequence_references))
    )

    ch_versions = ch_versions.mix(MUTATIONHUMAN.out.versions.first())
 
    //
    // MODULE: TABLELOOKUP
    //Check if mutations are annotated in mammalian and inhibition databases

    TABLELOOKUP  (
        MUTATIONHUMAN.out.inhibtion_mutation, Channel.value(file(params.inhibtion_mutation_db))
        
    )

    //
    // MODULE: REPORT
    //Combine all data into a report

    def runid = params.runid
    def seq_instrument   = params.seq_instrument  
    def release_version = params.release_version

    REPORTHUMAN  (
        SUBTYPEFINDER.out.subtype_report.collect(), 
        COVERAGE.out.coverage_report.collect(), 
        MUTATIONHUMAN.out.human_mutation_report.collect(), 
        MUTATIONHUMAN.out.inhibtion_mutation_report.collect(), 
        TABLELOOKUP.out.lookup_report.collect(),
        NEXTCLADE.out.nextclade_summary_rapport.collect(),
        NEXTCLADE.out.nextclade_report.collect(),
        MUTATIONHUMAN.out.vaccine_mutation_report.collect(),
        runid,
        release_version,
        COVERAGE.out.filtered_fasta_report.collect(),
        TECHNICAL.out.depth_files_report.collect(),
        seq_instrument,
        Channel.value(file(params.input))

    )
    
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
    ch_multiqc_files = ch_multiqc_files.mix(NEXTCLADE.out.versions.first().ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(SUBTYPEFINDER.out.versions.first().ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(MUTATIONHUMAN.out.versions.first().ifEmpty(null))
    



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

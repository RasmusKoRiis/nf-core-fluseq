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
include { SEGMENTIFENTIFIER           } from '../modules/local/blastnfasta/main'
include { GENOTYPING                  } from '../modules/local/genotyping/main'
include { COVERAGE                    } from '../modules/local/coverage/main'
include { FASTA_CONFIGURATION         } from '../modules/local/seqkit/main'
include { MUTATION                    } from '../modules/local/mutation/main'
include { TABLELOOKUP                 } from '../modules/local/tablelookup/main'
include { REPORT                      } from '../modules/local/report/main'





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// Function to parse a multi-FASTA file and channel sequences individually


def parseFasta(String fastaFilePath) {
    Channel
        .fromPath(fastaFilePath)
        .splitFasta(record: [id: true, seqString: true])
        .map { record ->
            def sampleId = record.id.tokenize('|')[0] // Extract sampleId
            def fileId = record.id // Full ID from the record
            def sequenceFile = "${fileId}.fasta" // Define file name with uniqueId to ensure uniqueness

            // Write sequence to file
            new File(sequenceFile).withWriter { writer ->
                writer.writeLine(">${record.id}")
                writer.writeLine(record.seqString)
            }

            // Construct metadata for this sequence
            def meta = [id: sampleId, sampleId: fileId]

            // Return a tuple of sampleId, metadata, and the file reference
            // This allows for grouping by sampleId while retaining detailed info for each sequence
            return [sampleId, [meta, file(sequenceFile)]]
        }
        .groupTuple()
        .map { sampleId, seqInfos ->
            // For each sampleId, seqInfos contains a list of [meta, file] tuples for each sequence
            
            // If you want to return a structure where each sample's sequences are together:
            // Flatten the list of files for easier access
            def files = seqInfos.collect { it[1] } // Collect all file references

            // Optionally, merge or manipulate metadata here. Assuming keeping first meta or creating a summary
            // For simplicity, let's just take the first sequence's metadata as a sample representation
            def meta = seqInfos[0][0]

            // Return a structure combining metadata with all sequence file paths for the sample
            [meta, files]
        }
}


workflow AVIANFASTA {

    //
    // INPUT PARSE
    //

    ch_sample_information = parseFasta(params.fasta)
    
    ch_versions = Channel.empty()

    def currentDir = System.getProperty('user.dir')
    def fullPath_references = "${currentDir}/${params.sequence_references}/references_2324.fasta"

    //
    // FIND SEGMENTS
    //

    SEGMENTIFENTIFIER (
        ch_sample_information, fullPath_references
    )

    /// SUBTYPE CHANNEL
    // SUBTYPE CHANNEL
    SEGMENTIFENTIFIER.out.fasta_segment
    .map { meta, files -> 
        def ha_files = files.findAll { it.getName().contains('HA') }
        def na_files = files.findAll { it.getName().contains('NA') }
        def pa_files = files.findAll { it.getName().contains('PA') }
        def pb1_files = files.findAll { it.getName().contains('PB1') }
        def pb2_files = files.findAll { it.getName().contains('PB2') }
        def ns_files = files.findAll { it.getName().contains('NS') }
        def np_files = files.findAll { it.getName().contains('NP') }
        def m_files = files.findAll { it.getName().contains('M') }
        return tuple(meta, ha_files ?: [], na_files ?: [], pa_files ?: [], pb1_files ?: [], pb2_files ?: [], ns_files ?: [], np_files ?: [], m_files ?: [])
    }
    .set { subtype_channel }



    /// GENOTYPING CHANNEL
    SEGMENTIFENTIFIER.out.fasta_segment
    .map { meta, files -> 
        def amended_consensus_files = files.findAll { it.getName().contains('.fa') }
        return (amended_consensus_files) ? tuple(meta, amended_consensus_files) : null
    }
    .filter { it != null }
    .set { genotyping_channel }

    //
    // MODULE: SUBTYPE FINDER
    //
    def fullPathHA = "${currentDir}/${params.ha_database}"
    def fullPathNA = "${currentDir}/${params.na_database}"


    SUBTYPEFINDER (
        subtype_channel, fullPathHA, fullPathNA
    )

    /// MUTATION CHANNELS

    SEGMENTIFENTIFIER.out.fasta_segment
    .join(SUBTYPEFINDER.out.subtype, by: [0]) // Assuming meta.id is the first element in the tuple
    .map { items ->
        def meta = items[0] // The common meta.id
        def fasta = items[1] // The fasta file from IRMA_out_amended_consensus
        def subtype = items[2] // The subtype file from SUBTYPEFINDER_out_subtype
        return [meta, fasta, subtype]
    }
    .set { fasta_subtype }


    /// FASTACONFIG CHANNELS

    SEGMENTIFENTIFIER.out.fasta_configuration
    .join(SUBTYPEFINDER.out.subtype, by: [0]) // Assuming meta.id is the first element in the tuple
    .map { items ->
        def meta = items[0] // The common meta.id
        def fasta = items[1] // The fasta file from IRMA_out_amended_consensus
        def subtype = items[2] // The subtype file from SUBTYPEFINDER_out_subtype
        return [meta, fasta, subtype]
    }
    .set { fastaconfig_subtype }


    //
    // MODULE: GENOTYPING
    //

    ch_genotype_database = params.genotype_database

    GENOTYPING (
        genotyping_channel, ch_genotype_database
    )

    //
    // MODULE: FASTA CONFIGURATION
    //

    
    FASTA_CONFIGURATION (
         fastaconfig_subtype 
    )

    //
    // MODULE: COVERAGE
    //

    
    COVERAGE (
         FASTA_CONFIGURATION.out.fasta
    )

    //
    // MODULE: AMINO ACID TRANSLATION
    //

    def fullPath_nextclade_dataset = "${currentDir}/${params.nextclade_dataset}"

    AMINOACIDTRANSLATION (
        COVERAGE.out.filtered_fasta, fullPath_nextclade_dataset
    )

    //
    // MODULE: MUTATION
    //

    def fullPath_references_2 = "${currentDir}/${params.sequence_references}"
    
    MUTATION  (
        AMINOACIDTRANSLATION.out.aminoacid_sequence, fullPath_references_2
    )
 
    //
    // MODULE: TABLELOOKUP
    //

    def fullPath_tables = "${currentDir}/${params.mutation_tables}"
    def fullPath_mammalian_mutation = "${currentDir}/${params.mamalian_mutation_db}"
    def fullPath_inhibtion_mutation = "${currentDir}/${params.inhibtion_mutation_db }"

    


    
    TABLELOOKUP  (
        MUTATION.out.mamailian_mutation, MUTATION.out.inhibtion_mutation, fullPath_mammalian_mutation, fullPath_inhibtion_mutation 
    )

    //
    // MODULE: REPORT
    //
    REPORT  (
        SUBTYPEFINDER.out.subtype_report.collect(), 
        GENOTYPING.out.genotype_report.collect(), 
        COVERAGE.out.coverage_report.collect(), 
        MUTATION.out.mamailian_mutation_report.collect(), 
        MUTATION.out.inhibtion_mutation_report.collect(), 
        TABLELOOKUP.out.lookup_report.collect()
    )
    

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

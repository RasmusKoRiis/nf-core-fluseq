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
include { TABLELOOKUP_MAMMALIAN       } from '../modules/local/tablelookup_mammalian/main'
include { REPORT_AVIAN                } from '../modules/local/report_avian/main'
include { FLUMUT                      } from '../modules/local/flumut/main'
include { FLUMUT_CONVERSION           } from '../modules/local/flumut_conversion/main'
include { GENIN2                      } from '../modules/local/genin2/main'





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// Function to parse a multi-FASTA file and channel sequences individually


/* ---------- imports & globals ---------- */
import groovy.transform.Field

@Field Map SEGMENT_INDEX = [
    'PB2':1, 'PB1':2, 'PA':3, 'HA':4,
    'NP':5,  'NA':6,  'MP':7, 'NS':8
]

/* duplicate counter, initialised once for the whole run */
@Field Map dupCounter = [:].withDefault{ 0 }

/* helper to make IDs filesystem-safe */
String sanitise( String s ) {
    s.replaceAll('[^A-Za-z0-9._-]', '-')   // swap illegal chars → -
     .replaceAll('-{2,}', '-')             // squeeze repeats
     .replaceAll(/^[-.]+|[-.]+$/, '')      // trim ends
}

import groovy.transform.Field

/* map segment → number */
@Field Map SEGMENT_NUM = [
    'PB2':1, 'PB1':2, 'PA':3, 'HA':4,
    'NP':5,  'NA':6,  'MP':7, 'NS':8
]

/* remember which (sample, segment) we have already written */
@Field Set taken = new HashSet<>()          // e.g. "A-Cambodia-…|HA"

String clean(String s) {
    s.replaceAll('[^A-Za-z0-9._-]', '-')
     .replaceAll('-{2,}', '-')
     .replaceAll(/^[-.]+|[-.]+$/, '')
}

def parseFasta(String fastaPath) {

    Channel
        .fromPath(fastaPath)
        .splitFasta(record:[id:true, seqString:true])
        .filter { rec ->                          // keep only first of each segment
            def seg = rec.id.split(/\|+/)
                         .find { SEGMENT_NUM.containsKey(it) } ?: 'UNK'
            def sid = clean(rec.id.split(/\|/)[0])
            taken.add("${sid}|${seg}")            // returns false if seen before
        }
        .map { rec ->
            def parts      = rec.id.split(/\|+/)
            def sampleId   = clean(parts[0])
            def segmentTag = parts.find { SEGMENT_NUM.containsKey(it) } ?: 'UNK'
            def segNum     = SEGMENT_NUM.getOrDefault(segmentTag, 0)

            /* simple name: ID_#.fasta */
            def fn = "${sampleId}_${segNum}.fasta"

            new File(fn).withWriter { w ->
                w << ">${sampleId}|${segmentTag}\n${rec.seqString}\n"
            }

            def meta = [ id: sampleId ]        // ← only the ID
            [ sampleId, [ meta, file(fn) ] ]
        }
        .groupTuple()
        .map { sampleId, seqInfos ->
            [ seqInfos[0][0], seqInfos.collect { it[1] } ]
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

    // SUBTYPE CHANNEL
 /* ---------- SUBTYPE CHANNEL (one file per segment) ---------- */
def pick = { fset, num ->                              // helper
    fset.find { it.getName().endsWith("_${num}.fasta") }
}

SEGMENTIFENTIFIER.out.fasta_segment
    .map { meta, files ->

        // mandatory segments
        def ha = pick(files, 4)
        def na = pick(files, 6)

        // skip sample if either HA or NA is missing
        if (!ha || !na) return null

        // optional segments
        def pa  = pick(files, 3)
        def pb1 = pick(files, 2)
        def pb2 = pick(files, 1)
        def ns  = pick(files, 8)
        def np  = pick(files, 5)
        def mp  = pick(files, 7)

        /* emit exactly nine elements: meta + eight paths (some may be null) */
        tuple(meta, ha, na, pa, pb1, pb2, ns, np, mp)
    }
    .filter { it != null }     // remove samples lacking HA or NA
    .set    { subtype_channel }


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

    def seq_quality_thershold = params.seq_quality_thershold

    COVERAGE (
        FASTA_CONFIGURATION.out.fasta, seq_quality_thershold
    )


    //
    // MODULE: FLUMUT
    //

    FLUMUT (
        FASTA_CONFIGURATION.out.fasta_flumut
    )

        
    //
    // MODULE: FLUMUT COVERSION
    //

    FLUMUT_CONVERSION (
        FLUMUT.out.markers
    )

    //
    // MODULE: GENIN2
    //

    GENIN2 (
        FASTA_CONFIGURATION.out.fasta_genin
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

    def fullPath_mammalian_mutation = "${currentDir}/${params.mamalian_mutation_db}"
    def fullPath_inhibtion_mutation = "${currentDir}/${params.inhibtion_mutation_db }"

    TABLELOOKUP  (
        MUTATION.out.inhibtion_mutation, fullPath_inhibtion_mutation
    )

    TABLELOOKUP_MAMMALIAN  (
        MUTATION.out.mamailian_mutation, fullPath_mammalian_mutation
    )

    //
    // MODULE: REPORT
    //
    REPORT_AVIAN  (
        SUBTYPEFINDER.out.subtype_report.collect(), 
        GENOTYPING.out.genotype_report.collect(), 
        COVERAGE.out.coverage_report.collect(), 
        MUTATION.out.mamailian_mutation_report.collect(), 
        TABLELOOKUP.out.lookup_report.collect(),
        TABLELOOKUP_MAMMALIAN.out.lookup_report.collect()
    )
    

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

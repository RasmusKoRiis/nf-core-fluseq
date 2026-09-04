/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap } from 'plugin/nf-schema'

def summary_params = paramsSummaryMap(workflow)

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
include { COVERAGE                    } from '../modules/local/coverage/main'
include { FASTA_CONFIGURATION         } from '../modules/local/seqkit/main'
include { MUTATION                    } from '../modules/local/mutation/main'
include { TABLELOOKUP_MAMMALIAN                } from '../modules/local/tablelookup_mammalian/main'
include { REPORT_AVIAN                } from '../modules/local/report_avian/main'
include { FLUMUT                      } from '../modules/local/flumut/main'
include { FLUMUT_CONVERSION           } from '../modules/local/flumut_conversion/main'
include { GENIN2                      } from '../modules/local/genin2/main'
include { REASSORTMENT                } from '../modules/local/reassortment/main'
include { SURVEILLANCE_SUMMARY        } from '../modules/local/surveillance_summary/main'
include { REFERENCE_PROVENANCE        } from '../modules/local/reference_provenance/main'







/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


workflow AVIAN {

    //
    // INPUT PARSE
    //

    ch_versions = Channel.empty()

    REFERENCE_PROVENANCE(Channel.value([
        file(params.ha_database, checkIfExists: true),
        file(params.na_database, checkIfExists: true),
        file(params.genotype_database, checkIfExists: true),
        file(params.reassortment_database, checkIfExists: true),
        file(params.mamalian_mutation_db ?: params.mammalian_mutation_db, checkIfExists: true),
        file(params.inhibtion_mutation_db ?: params.inhibition_mutation_db, checkIfExists: true),
        file(params.sequence_references, checkIfExists: true),
        file(params.nextclade_dataset, checkIfExists: true)
    ]))

    INPUT_CHECK(Channel.fromPath(params.input, checkIfExists: true))
    INPUT_CHECK.out.reads.set { read_input }

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
        def pa_files = files.findAll { it.getName().contains('_PA') }
        def pb1_files = files.findAll { it.getName().contains('_PB1') }
        def pb2_files = files.findAll { it.getName().contains('_PB2') }
        def ns_files = files.findAll { it.getName().contains('_NS') }
        def np_files = files.findAll { it.getName().contains('_NP') }
        def m_files = files.findAll { it.getName().contains('_M') }
        return (ha_files && na_files && pa_files && pb1_files && pb2_files && ns_files && np_files && m_files) ? tuple(meta, ha_files, na_files, pa_files, pb1_files, pb2_files, ns_files, np_files, m_files) : null
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
    SUBTYPEFINDER (
        IRMA_ha_na_fasta,
        Channel.value(file(params.ha_database, checkIfExists: true)),
        Channel.value(file(params.na_database, checkIfExists: true))
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


    //
    // MODULE: GENOTYPING
    //

    GENOTYPING (
        IRMA_amended_consensus_files, Channel.value(file(params.genotype_database, checkIfExists: true))
    )


    //
    // MODULE: FASTA CONFIGURATION
    //

    
    FASTA_CONFIGURATION (
         fasta_subtype 
    )

    REASSORTMENT(
        FASTA_CONFIGURATION.out.fasta_flumut,
        Channel.value(file(params.reassortment_database, checkIfExists: true))
    )


    //
    // MODULE: FLUMUT
    //

    FLUMUT (
        FASTA_CONFIGURATION.out.fasta_flumut
    )

    //
    // MODULE: GENIN2
    //

    GENIN2 (
        FASTA_CONFIGURATION.out.fasta_genin
    )



    //
    // MODULE: COVERAGE
    //

    /// Coverage threshold from the params/user
    def seq_quality_threshold = params.seq_quality_thershold ?: params.seq_quality_threshold

    
    COVERAGE (
         FASTA_CONFIGURATION.out.fasta, seq_quality_threshold
    )


    //
    // MODULE: AMINO ACID TRANSLATION
    //

    AMINOACIDTRANSLATION (
        COVERAGE.out.filtered_fasta, Channel.value(file(params.nextclade_dataset, checkIfExists: true))
    )
   

    //
    // MODULE: MUTATION
    //

    MUTATION  (
        AMINOACIDTRANSLATION.out.aminoacid_sequence, Channel.value(file(params.sequence_references, checkIfExists: true))
    )

    //
    // MODULE: TABLELOOKUP
    //

    TABLELOOKUP_MAMMALIAN  (
        AMINOACIDTRANSLATION.out.mutation_lookup_csv, Channel.value(file(params.mamalian_mutation_db ?: params.mammalian_mutation_db, checkIfExists: true))
    )


    SURVEILLANCE_SUMMARY(
        params.file,
        FASTA_CONFIGURATION.out.fasta_flumut.map { meta, fasta -> fasta }.collect(),
        COVERAGE.out.coverage_report.collect(),
        SUBTYPEFINDER.out.subtype_report.collect(),
        SUBTYPEFINDER.out.subtype_hits
            .map { meta, ha_hits, na_hits -> [ha_hits, na_hits] }
            .flatten()
            .collect(),
        REASSORTMENT.out.genotype_report.collect(),
        TABLELOOKUP_MAMMALIAN.out.lookup_report.collect(),
        Channel.value([file(params.mamalian_mutation_db ?: params.mammalian_mutation_db, checkIfExists: true)]),
        Channel.value(file("$projectDir/bin/surveillance_summary.py", checkIfExists: true))
    )

    //
    // MODULE: REPORT
    //
    def runid = params.runid


    REPORT_AVIAN  (
        SUBTYPEFINDER.out.subtype_report.collect(), 
        GENOTYPING.out.genotype_report.collect(), 
        COVERAGE.out.coverage_report.collect(),
        TABLELOOKUP_MAMMALIAN.out.lookup_report.collect(),
        AMINOACIDTRANSLATION.out.nextclade_csv.collect(),
        runid 
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

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
===============================================================================
UID-FIRST FASTA STRATEGY (FULL PIPELINE → REPORTHUMANFASTA)
- Robust to mixed headers; stable UID per biological “core”
- Never call list ops (e.g., collect) on GroupTupleOp
- Build segments from our own per-UID FASTAs (avoid grouped module outputs)
- Maintain UID ↔ OriginalName mapping (id_map.tsv) for final report
===============================================================================
*/

/* ──────────────────────────────────────────────────────────────────────────
   MODULES (paths follow your original layout)
   ────────────────────────────────────────────────────────────────────────── */
include { AMINOACIDTRANSLATION } from '../modules/local/aminoacidtranslation/main'
include { SUBTYPEFINDER        } from '../modules/local/blastn/main'
include { SEGMENTIFENTIFIER    } from '../modules/local/blastnfasta/main'
include { GENOTYPING           } from '../modules/local/genotyping/main'
include { COVERAGE             } from '../modules/local/coverage/main'
include { FASTA_CONFIGURATIONFASTA  } from '../modules/local/seqkitfasta/main'
include { MUTATIONHUMAN        } from '../modules/local/mutationhuman/main'
include { REPORTAVIANFASTA     } from '../modules/local/reportavianfasta/main'
include { NEXTCLADE            } from '../modules/local/nextclade/main'
include { REASSORTMENT         } from '../modules/local/reassortment/main'
include { MUTATION                    } from '../modules/local/mutation/main'
include { TABLELOOKUP                 } from '../modules/local/tablelookup/main'
include { TABLELOOKUP_MAMMALIAN       } from '../modules/local/tablelookup_mammalian/main'
include { FLUMUT                      } from '../modules/local/flumut/main'
include { FLUMUT_CONVERSION           } from '../modules/local/flumut_conversion/main'
include { GENIN2                      } from '../modules/local/genin2/main'
include { SURVEILLANCE_SUMMARY         } from '../modules/local/surveillance_summary/main'
include { REFERENCE_PROVENANCE         } from '../modules/local/reference_provenance/main'
include { EMIT_FASTA_RECORD; WRITE_ID_MAP; REHEADER_TO_UID } from '../modules/local/fasta_records/main'
include { SLIM_GENIN2_REPORT } from '../modules/local/slim_genin2_report/main'

/* ──────────────────────────────────────────────────────────────────────────
   MAIN WORKFLOW
   ────────────────────────────────────────────────────────────────────────── */
workflow AVIANFASTA {

  if ( !params.fasta ) error "Missing required parameter: --fasta"

  // Refs/DBs
  def ref_fasta             = file("${params.sequence_references}/references_2324.fasta", checkIfExists: true)
  def ref_dir_all           = file(params.sequence_references, checkIfExists: true)
  def ha_db                 = file(params.ha_database, checkIfExists: true)
  def na_db                 = file(params.na_database, checkIfExists: true)
  def genotype_db           = file(params.genotype_database, checkIfExists: true)
  def nextclade_dataset_dir = file(params.nextclade_dataset, checkIfExists: true)
  def inhib_mut_db          = file(params.inhibition_mutation_db, checkIfExists: true)
  def reassortment_db       = file(params.reassortment_database, checkIfExists: true)

  REFERENCE_PROVENANCE(Channel.value([
    ref_dir_all,
    ha_db,
    na_db,
    genotype_db,
    nextclade_dataset_dir,
    inhib_mut_db,
    file(params.mammalian_mutation_db, checkIfExists: true),
    reassortment_db
  ]))

  /* 1) Split multi-FASTA → per-record UID FASTAs */
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .splitFasta(record: [ id: true, seqString: true ])
    .map { rec ->
      def hdrFull  = rec.id as String
      def core     = FastaUtils.coreFromHeader(hdrFull)
      def uid      = FastaUtils.uidFromCore(core)
      def fileStem = hdrFull.replaceAll(/[^A-Za-z0-9._-]+/, '_')
      tuple(uid, uid, fileStem, rec.seqString, core)
    }
    .set { ch_records }

  EMIT_FASTA_RECORD( ch_records )

  /* 2) UID ↔ OriginalName map (materialize once) */
  EMIT_FASTA_RECORD.out
    .map { uid, core, f -> tuple(uid?.toString()?.trim(), core?.toString()?.trim()) }
    .distinct()
    .toList()
    .set { ch_id_pairs_list }

  WRITE_ID_MAP( ch_id_pairs_list )
  WRITE_ID_MAP.out.id_map.set { ch_id_map_file }

  /* 3) Per-sample bundles (meta + files) — plain tuples only */
  EMIT_FASTA_RECORD.out
    .groupTuple() // -> uid, [core], [files]
    .map { uid, cores, files -> tuple([ id: uid, orig: (cores ? cores[0] : uid) ], files) }
    .set { ch_sample_info }

  /* 4) Optionally run SEGMENTIFENTIFIER (not consumed for grouping here) */
  SEGMENTIFENTIFIER( ch_sample_info, ref_fasta )


  /* 5) Build segment groups from our own files (avoid GroupTupleOp traps) */
  ch_sample_info
    .map { meta, files -> tuple(meta, files.findAll { f -> FastaUtils.filename(f) ==~ /(?i).*\.fa(sta)?$/ }) }
    .filter { meta, files -> files && files.size() > 0 }
    .set { ch_segments_grouped }


  /* 6) SUBTYPING: pick HA + NA by filename tokens */
  ch_segments_grouped
    .map { meta, files ->
      def pick = { String seg ->
        files.find { f ->
          def fname = FastaUtils.filename(f)
          def toks  = fname.replaceFirst(/\.[^.]+$/, '').split(/[^A-Za-z0-9]+/)*.toUpperCase()
          toks.contains(seg.toUpperCase())
        }
      }
      def ha = pick('HA'); def na = pick('NA')
      (ha && na) ? tuple(meta, ha, na) : null
    }
    .filter { it != null }
    .set { ch_subtype_pairs }

  REHEADER_TO_UID( ch_subtype_pairs )
  SUBTYPEFINDER( REHEADER_TO_UID.out, ha_db, na_db )

  /* 7) GENOTYPING: all FASTAs per sample */
  ch_segments_grouped
    .map { meta, files ->
      def fas = files.findAll { f -> FastaUtils.filename(f) ==~ /(?i).*\.fa(sta)?$/ }
      fas ? tuple(meta, fas) : null
    }
    .filter { it != null }
    .set { ch_genotyping }



  /* 8) FASTA CONFIGURATION input = (meta, files, subtype) */
    ch_segments_grouped
    .join(SUBTYPEFINDER.out.subtype, by: [0])   // -> meta, files, subtype_file
    .set { ch_segments_with_subtype }

    // Call the new robust module
    FASTA_CONFIGURATIONFASTA( ch_segments_with_subtype )
    GENOTYPING(FASTA_CONFIGURATIONFASTA.out.fasta_genotyping, genotype_db )

  /* 9) REASSORTMENT */
  REASSORTMENT( FASTA_CONFIGURATIONFASTA.out.fasta_flumut, Channel.value(reassortment_db) )

  /* 10) Coverage */
  COVERAGE( FASTA_CONFIGURATIONFASTA.out.fasta, params.seq_quality_threshold )

  /* 11) Nextclade */
  NEXTCLADE(
    COVERAGE.out.filtered_fasta,
    Channel.value(nextclade_dataset_dir)
  )

  /* 12) Mutation vs references */
  //MUTATIONHUMAN( NEXTCLADE.out.aminoacid_sequence, Channel.value(ref_dir_all) )

  /* 13) Table lookups */
  //TABLELOOKUP( MUTATIONHUMAN.out.inhibtion_mutation, Channel.value(inhib_mut_db) )

  //
  // MODULE: FLUMUT
  //

  FLUMUT (
       FASTA_CONFIGURATIONFASTA.out.fasta_flumut
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
        FASTA_CONFIGURATIONFASTA.out.fasta_genin
  )

  SLIM_GENIN2_REPORT( GENIN2.out.genin2_report )

  //
  // MODULE: AMINO ACID TRANSLATION
  //

  def fullPath_nextclade_dataset           = file(params.nextclade_dataset, checkIfExists: true)
  def fullPath_references_2                = file(params.sequence_references, checkIfExists: true)
  def fullPath_mammalian_mutation          = file(params.mammalian_mutation_db, checkIfExists: true)
  def fullPath_inhibition_mutation          = file(params.inhibition_mutation_db, checkIfExists: true)

  AMINOACIDTRANSLATION (
      COVERAGE.out.filtered_fasta, fullPath_nextclade_dataset
  )

  //
  // MODULE: MUTATION
  //


    
  MUTATION  (
      AMINOACIDTRANSLATION.out.aminoacid_sequence, fullPath_references_2
  )
 
  //
  // MODULE: TABLELOOKUP
  //

  def ch_full_mutation_lists = MUTATION.out.full_mutation_list
  def ch_full_inhib = ch_full_mutation_lists.filter { meta, f, subtype -> FastaUtils.filename(f) ==~ /.*_inhibtion_.*/ }
  def ch_full_mamm  = ch_full_mutation_lists.filter { meta, f, subtype -> FastaUtils.filename(f) ==~ /.*_mamailian_.*/ }


  TABLELOOKUP  (
      ch_full_inhib, fullPath_inhibition_mutation
  )

  TABLELOOKUP_MAMMALIAN  (
      ch_full_mamm, fullPath_mammalian_mutation
  )

  SURVEILLANCE_SUMMARY(
    params.file,
    FASTA_CONFIGURATIONFASTA.out.fasta_flumut.map { meta, fasta -> fasta }.collect(),
    COVERAGE.out.coverage_report.collect(),
    SUBTYPEFINDER.out.subtype_report.collect(),
    SUBTYPEFINDER.out.subtype_hits
      .map { meta, ha_hits, na_hits -> [ha_hits, na_hits] }
      .flatten()
      .collect(),
    REASSORTMENT.out.genotype_report.collect(),
    TABLELOOKUP.out.lookup_report
      .mix(TABLELOOKUP_MAMMALIAN.out.lookup_report)
      .mix(FLUMUT_CONVERSION.out.flumut_report)
      .collect(),
    Channel.value([
      file(params.inhibition_mutation_db, checkIfExists: true),
      file(params.mammalian_mutation_db, checkIfExists: true)
    ]),
    Channel.value(file("$projectDir/bin/surveillance_summary.py", checkIfExists: true))
  )

  /* 14) Report (materialize leaf streams only) */
  REPORTAVIANFASTA(
    SUBTYPEFINDER.out.subtype_report.collect(),
    COVERAGE.out.coverage_report.collect(),
    MUTATION.out.mamailian_mutation_report.collect(),
    MUTATION.out.vaccine_mutation_report.collect(),
    MUTATION.out.full_mutation_list_report.collect(),
    TABLELOOKUP.out.lookup_report.collect(),
    TABLELOOKUP_MAMMALIAN.out.lookup_report.collect(),
    NEXTCLADE.out.nextclade_summary_rapport.collect(),
    NEXTCLADE.out.nextclade_report.collect(),
    ch_id_map_file,                                // id_map.tsv (Path)
    params.runid,
    params.release_version,
    COVERAGE.out.filtered_fasta_report.collect(),
    params.seq_instrument,
    Channel.value(file(params.input ?: params.fasta)),
    SLIM_GENIN2_REPORT.out.genin2_report_trimmed.collect(),
    FLUMUT_CONVERSION.out.flumut_report.collect(),
    REASSORTMENT.out.genotype_report.collect()
  )

}

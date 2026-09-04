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
include { SUBTYPEFINDER        } from '../modules/local/blastn/main'
include { SEGMENTIFENTIFIER    } from '../modules/local/blastnfasta/main'
include { GENOTYPING           } from '../modules/local/genotyping/main'
include { COVERAGE             } from '../modules/local/coverage/main'
include { FASTA_CONFIGURATIONFASTA  } from '../modules/local/seqkitfasta/main'
include { MUTATIONHUMAN        } from '../modules/local/mutationhuman/main'
include { TABLELOOKUP          } from '../modules/local/tablelookup/main'
include { DRUG_RESISTANCE_REPORT } from '../modules/local/drug_resistance_report/main'
include { REPORTHUMANFASTA     } from '../modules/local/reporthumanfasta/main'
include { NEXTCLADE            } from '../modules/local/nextclade/main'
include { SUBCLADE_NOMENCLATURE; SUBCLADE_NOMENCLATURE_RULES } from '../modules/local/subclade_nomenclature/main'
include { REASSORTMENT         } from '../modules/local/reassortment/main'
include { SURVEILLANCE_SUMMARY } from '../modules/local/surveillance_summary/main'
include { REFERENCE_PROVENANCE } from '../modules/local/reference_provenance/main'
include { EMIT_FASTA_RECORD; WRITE_ID_MAP; REHEADER_TO_UID } from '../modules/local/fasta_records/main'

/* ──────────────────────────────────────────────────────────────────────────
   MAIN WORKFLOW
   ────────────────────────────────────────────────────────────────────────── */
workflow HUMANFASTA {

  if ( !params.fasta ) error "Missing required parameter: --fasta"

  def drugResistanceOnly = params.drug_resistance_only.toString().toBoolean()

  // Refs/DBs
  def ref_dir_all           = file(params.sequence_references, checkIfExists: true)
  def ha_db                 = file(params.ha_database, checkIfExists: true)
  def na_db                 = file(params.na_database, checkIfExists: true)
  def inhib_mut_db          = file(params.inhibition_mutation_db, checkIfExists: true)
  def ref_fasta             = drugResistanceOnly ? null : file("${params.sequence_references}/references_2324.fasta", checkIfExists: true)
  def genotype_db           = drugResistanceOnly ? null : file(params.genotype_database, checkIfExists: true)
  def reassortment_db       = drugResistanceOnly ? null : file(params.reassortment_database, checkIfExists: true)

  def provenance_refs = [ref_dir_all, ha_db, na_db, inhib_mut_db, file(params.nextclade_dataset, checkIfExists: true)]
  if (!drugResistanceOnly) {
    provenance_refs.addAll([genotype_db, reassortment_db])
  }
  REFERENCE_PROVENANCE(Channel.value(provenance_refs))

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

  /* 2) UID ↔ OriginalName map for both full and resistance-only reports */
  ch_id_pairs_list = EMIT_FASTA_RECORD.out
    .map { uid, core, f -> tuple(uid?.toString()?.trim(), core?.toString()?.trim()) }
    .distinct()
    .toList()

  WRITE_ID_MAP( ch_id_pairs_list )
  ch_id_map_file = WRITE_ID_MAP.out.id_map

  /* 3) Per-sample bundles (meta + files) — plain tuples only */
  EMIT_FASTA_RECORD.out
    .groupTuple() // -> uid, [core], [files]
    .map { uid, cores, files -> tuple([ id: uid, orig: (cores ? cores[0] : uid) ], files) }
    .set { ch_sample_info }

  /* 4) SEGMENTIFENTIFIER is informational and is not needed by resistance-only mode. */
  if ( !drugResistanceOnly ) {
    SEGMENTIFENTIFIER( ch_sample_info, ref_fasta )
  }

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

  /* 7) GENOTYPING: full workflow only */
  def ch_genotyping = null
  if ( !drugResistanceOnly ) {
    ch_genotyping = ch_segments_grouped
      .map { meta, files ->
        def fas = files.findAll { f -> FastaUtils.filename(f) ==~ /(?i).*\.fa(sta)?$/ }
        fas ? tuple(meta, fas) : null
      }
      .filter { it != null }

    GENOTYPING( ch_genotyping, genotype_db )
  }

  /* 8) FASTA CONFIGURATION input = (meta, files, subtype) */
    ch_segments_grouped
    .join(SUBTYPEFINDER.out.subtype, by: [0])   // -> meta, files, subtype_file
    .set { ch_segments_with_subtype }

    // Call the new robust module
    FASTA_CONFIGURATIONFASTA( ch_segments_with_subtype )

  /* 9-11) Select the smallest valid input path to Nextclade. */
  def ch_nextclade_input
  if ( drugResistanceOnly ) {
    // Drug resistance uses only NA, PA and M2. The M segment is required to
    // produce the M2 translation, while HA/NA were already used for subtyping.
    ch_nextclade_input = FASTA_CONFIGURATIONFASTA.out.fasta
      .map { meta, fasta, subtype ->
        def resistanceFasta = FastaUtils.flatten([fasta]).findAll { f ->
          def name = FastaUtils.filename(f)?.toUpperCase() ?: ''
          name.contains('-NA-') || name.contains('-PA-') || name.contains('-MP-')
        }
        resistanceFasta ? tuple(meta, resistanceFasta, subtype) : null
      }
      .filter { it != null }
  } else {
    REASSORTMENT( FASTA_CONFIGURATIONFASTA.out.fasta_flumut, Channel.value(reassortment_db) )

    COVERAGE( FASTA_CONFIGURATIONFASTA.out.fasta, params.seq_quality_threshold )

    ch_nextclade_input = COVERAGE.out.filtered_fasta
      .map { meta, fasta, subtype, coverage_csv -> tuple(meta, fasta, subtype) }

    SUBCLADE_NOMENCLATURE_RULES()
    ch_subclade_nomenclature_rules = SUBCLADE_NOMENCLATURE_RULES.out.rules_dir.first()
    ch_subclade_nomenclature_script = Channel.value(file("$projectDir/bin/subclade_nomenclature.py", checkIfExists: true))
    ch_characterisation_script = Channel.value(file("$projectDir/bin/characterise_reference_virus.py", checkIfExists: true))
    ch_characterisation_guidelines = Channel.value(file("$projectDir/assets/characterisation_guidelines", checkIfExists: true))

    SUBCLADE_NOMENCLATURE(
      COVERAGE.out.filtered_fasta,
      ch_subclade_nomenclature_rules,
      ch_subclade_nomenclature_script,
      ch_characterisation_script,
      ch_characterisation_guidelines
    )
  }

  NEXTCLADE(
    ch_nextclade_input,
    Channel.value(file(params.nextclade_dataset, checkIfExists: true))
  )

  /* 12) Mutation vs references */
  MUTATIONHUMAN( NEXTCLADE.out.aminoacid_sequence, Channel.value(ref_dir_all) )

  /* 13) Table lookups */
  TABLELOOKUP( MUTATIONHUMAN.out.inhibtion_mutation, Channel.value(inhib_mut_db) )

  if ( drugResistanceOnly ) {
    DRUG_RESISTANCE_REPORT(
      SUBTYPEFINDER.out.subtype_report.collect(),
      TABLELOOKUP.out.lookup_report.collect(),
      ch_id_map_file,
      params.runid
    )
  } else {
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
      TABLELOOKUP.out.lookup_report.collect(),
      Channel.value([file(params.inhibition_mutation_db, checkIfExists: true)]),
      Channel.value(file("$projectDir/bin/surveillance_summary.py", checkIfExists: true))
    )

    /* 14) Report (materialize leaf streams only) */
    REPORTHUMANFASTA(
      SUBTYPEFINDER.out.subtype_report.collect(),
      COVERAGE.out.coverage_report.collect(),
      MUTATIONHUMAN.out.human_mutation_report.collect(),
      MUTATIONHUMAN.out.inhibtion_mutation_report.collect(),
      TABLELOOKUP.out.lookup_report.collect(),
      NEXTCLADE.out.nextclade_summary_rapport.collect(),
      NEXTCLADE.out.nextclade_report.collect(),
      MUTATIONHUMAN.out.vaccine_mutation_report.collect(),
      ch_id_map_file,                                // id_map.tsv (Path)
      params.runid,
      params.release_version,
      COVERAGE.out.filtered_fasta_report.collect(),
      params.seq_instrument,
      Channel.value(file(params.input ?: params.fasta)),
      REASSORTMENT.out.genotype_report.collect(),
      SUBCLADE_NOMENCLATURE.out.report.collect()
    )
  }

  if ( drugResistanceOnly ) {
    log.info 'Running the human FASTA drug-resistance path.'
  }
}

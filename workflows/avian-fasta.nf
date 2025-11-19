nextflow.enable.dsl = 2

/*
===============================================================================
UID-FIRST FASTA STRATEGY (FULL PIPELINE → REPORTHUMANFASTA)
- Robust to mixed headers; stable UID per biological “core”
- Never call list ops (e.g., collect) on GroupTupleOp
- Build segments from our own per-UID FASTAs (avoid grouped module outputs)
- Maintain UID ↔ OriginalName mapping (id_map.tsv) for final report
===============================================================================
*/

params.fasta                      = params.fasta ?: ''
params.sequence_references        = params.sequence_references ?: './refs'
params.ha_database                = params.ha_database ?: ''
params.na_database                = params.na_database ?: ''
params.genotype_database          = params.genotype_database ?: ''
params.nextclade_dataset          = params.nextclade_dataset ?: ''
params.inhibtion_mutation_db      = params.inhibtion_mutation_db ?: ''
params.reassortment_database      = params.reassortment_database ?: ''
params.seq_quality_thershold      = params.seq_quality_thershold ?: 0.0
params.runid                      = params.runid ?: 'unknown'
params.seq_instrument             = 'unknown'
params.release_version            = params.release_version ?: 'unknown'
params.input                      = params.input ?: ''

/* ──────────────────────────────────────────────────────────────────────────
   MODULES (paths follow your original layout)
   ────────────────────────────────────────────────────────────────────────── */
include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'

include { IRMA                 } from '../modules/local/irma/main'
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

/* ──────────────────────────────────────────────────────────────────────────
   PROCESSES
   ────────────────────────────────────────────────────────────────────────── */

/* Single-record FASTA writer */
process EMIT_FASTA_RECORD {
  tag     { sample_id }
  cpus    1
  memory  '512 MB'
  scratch true

  input:
  tuple val(sample_id), val(header_id), val(file_stem), val(seq), val(orig_name)

  output:
  tuple val(sample_id), val(orig_name), path("${file_stem}.fasta")

  script:
  """
  cat > ${file_stem}.fasta <<'EOF'
  >${header_id}
  ${seq}
  EOF
  """
}

/* UID ↔ OriginalName mapping TSV */
process WRITE_ID_MAP {
  tag 'id_map'
  cpus 1
  memory '256 MB'

  input:
  val pairs  // List<[uid, orig]>

  output:
  path('id_map.tsv'), emit: id_map

  script:
  def lines = pairs.collect { row ->
    def uid  = (row[0] ?: '').toString().trim()
    def orig = (row[1] ?: '').toString().trim()
    "${uid}\t${orig}"
  }.join('\n')
  """
  printf 'SampleID\\tOriginalName\\n' > id_map.tsv
  cat >> id_map.tsv <<'EOF'
  ${lines}
  EOF
  """
}

/* Reheader HA/NA to UID */
process REHEADER_TO_UID {
  tag { meta.id }
  cpus 1
  memory '256 MB'

  input:
  tuple val(meta), path(ha), path(na)

  output:
  tuple val(meta), path("HA_${meta.id}.fa"), path("NA_${meta.id}.fa")

  script:
  """
  awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${ha} > HA_${meta.id}.fa
  awk 'BEGIN{h=">${meta.id}"} /^>/{print h; next} {print}' ${na} > NA_${meta.id}.fa
  """
}

/* ──────────────────────────────────────────────────────────────────────────
   GROOVY HELPERS
   ────────────────────────────────────────────────────────────────────────── */
def coreFromHeader = { String hdr ->
  // tolerate "01-HA|A/Poland/28/2024" or "A/Poland/28/2024|meta"
  def parts = hdr.split(/\|/, 2)
  def L = parts[0]; def R = parts.size() > 1 ? parts[1] : ''
  def cnt = { String x -> (x.findAll('/')?.size() ?: 0) }
  def best = (cnt(L) >= cnt(R) ? L : R).trim()
  if (best ==~ /^\d{2}-[A-Za-z0-9]+$/ && parts.size() > 1) return R.trim()
  best
}

def normCore  = { String s -> s?.toUpperCase()?.replaceAll(/[^A-Z0-9]/, '') ?: '' }

def shortHash = { String s ->
  java.security.MessageDigest.getInstance('SHA-1')
    .digest((s ?: '').bytes).encodeHex().toString().substring(0,8).toUpperCase()
}

def uidFromCore = { String core ->
  def base = normCore(core)
  ((base ?: 'S') + shortHash(base))
}

def flattenAny = { xs ->
  def out = []
  xs.each { x -> out.addAll( x instanceof List || x instanceof Object[] ? x : [x] ) }
  out
}

/* Path-safe filename */
def fnameOf = { obj ->
  if (obj instanceof java.nio.file.Path) return obj.getFileName().toString()
  if (obj instanceof File)               return obj.getName()
  return obj?.toString()
}

/* Extract a core-like key from filename tokens, e.g. "02-NA_A_City_123_2024.fa" -> "A/City/123/2024" */
def coreKeyFromFilename = { String name ->
  def base = name.replaceFirst(/\.[^.]+$/, '')
  def toks = base.split(/[^A-Za-z0-9]+/)
  int i    = toks.findIndexOf { it.equalsIgnoreCase('A') || it.equalsIgnoreCase('B') }
  def tail = (i >= 0 ? toks[i..-1] : toks)
  tail.join('/')
}

/* ──────────────────────────────────────────────────────────────────────────
   MAIN WORKFLOW
   ────────────────────────────────────────────────────────────────────────── */
workflow AVIANFASTA {

  if ( !params.fasta ) error "Missing required parameter: --fasta"

  // Refs/DBs
  def ref_fasta             = file("${params.sequence_references}/references_2324.fasta")
  def ref_dir_all           = file(params.sequence_references)
  def ha_db                 = file(params.ha_database)
  def na_db                 = file(params.na_database)
  def genotype_db           = file(params.genotype_database)
  def nextclade_dataset_dir = file(params.nextclade_dataset) // used inside module
  def inhib_mut_db          = file(params.inhibtion_mutation_db)
  def reassortment_db       = file(params.reassortment_database)

  /* 1) Split multi-FASTA → per-record UID FASTAs */
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .splitFasta(record: [ id: true, seqString: true ])
    .map { rec ->
      def hdrFull  = rec.id as String
      def core     = coreFromHeader(hdrFull)
      def uid      = uidFromCore(core)
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
    .map { meta, files -> tuple(meta, files.findAll { f -> fnameOf(f) ==~ /(?i).*\.fa(sta)?$/ }) }
    .filter { meta, files -> files && files.size() > 0 }
    .set { ch_segments_grouped }


  /* 6) SUBTYPING: pick HA + NA by filename tokens */
  ch_segments_grouped
    .map { meta, files ->
      def pick = { String seg ->
        files.find { f ->
          def fname = fnameOf(f)
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
      def fas = files.findAll { f -> fnameOf(f) ==~ /(?i).*\.fa(sta)?$/ }
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
  COVERAGE( FASTA_CONFIGURATIONFASTA.out.fasta, params.seq_quality_thershold )

  /* 11) Nextclade */
  NEXTCLADE( COVERAGE.out.filtered_fasta )

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

  GENIN2.out.genin2_report.view()


  //
  // MODULE: AMINO ACID TRANSLATION
  //

  def fullPath_nextclade_dataset           = file(params.nextclade_dataset)
  def fullPath_references_2                = file(params.sequence_references) 
  def fullPath_mammalian_mutation          = file(params.mamalian_mutation_db)
  def fullPath_inhibtion_mutation          = file(params.inhibtion_mutation_db)

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


  TABLELOOKUP  (
      MUTATION.out.inhibtion_mutation, fullPath_inhibtion_mutation
  )

  TABLELOOKUP_MAMMALIAN  (
      MUTATION.out.mamailian_mutation, fullPath_mammalian_mutation
  )

  /* 14) Report (materialize leaf streams only) */
  REPORTAVIANFASTA(
    SUBTYPEFINDER.out.subtype_report.collect(),
    COVERAGE.out.coverage_report.collect(),
    MUTATION.out.mamailian_mutation_report.collect(), 
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
    GENIN2.out.genin2_report.collect(),
    FLUMUT_CONVERSION.out.flumut_report.collect()
  )

  /* Bring-up logs */
  ch_subtype_pairs.count().subscribe   { n -> log.info "🔎 HA/NA subtype inputs: ${n}" }
  ch_genotyping.count().subscribe      { n -> log.info "🧬 Genotyping inputs: ${n}" }
  ch_id_map_file.view                  { f -> "📄 Wrote ID map: ${f}" }
}

/* Entrypoint */
workflow {
  AVIANFASTA()
}

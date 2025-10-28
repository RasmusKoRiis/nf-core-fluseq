process FASTA_CONFIGURATIONFASTA {

  tag     { meta.id }
  label   'fastaconfig'
  cpus    1
  memory  '1 GB'
  scratch true

  /*
    Inputs:
      - meta: [ id: UID, orig: core-like string ]    (not used inside bash)
      - segments: list of PATHs (.fa/.fasta) for this sample (staged into work dir)
      - subtype_txt: text file named like ${meta.id}_subtype.txt (staged into work dir)
  */
  input:
  tuple val(meta), path(segments), path(subtype_txt)

  output:
    tuple val(meta), path("*.fa"), path(subtype_txt), emit: fasta
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta_merge
    tuple val(meta), path("${meta.id}_flumut.fasta"), emit: fasta_flumut
    tuple val(meta), path("${meta.id}_genin.fasta"), emit: fasta_genin

  /*
    Strategy:
    - Do NOT reference ${meta.id}, ${segments}, ${subtype_txt} in the script.
    - Work only with staged files in CWD to avoid any Groovy interpolation issues.
  */
  script: '''
  set -euo pipefail

  # Find subtype file and derive META_ID
  subfile="$(ls *_subtype.txt 2>/dev/null | head -n1 || true)"
  if [ -z "$subfile" ]; then
    echo "No *_subtype.txt found in work dir" >&2
    exit 1
  fi
  META_ID="${subfile%_subtype.txt}"

  # Read subtype (first non-empty token)
  tr -d "\r" < "$subfile" | awk "NF{print; exit}" > subtype.tmp
  subtype="$(cat subtype.tmp)"

  # Map a segment token to numbered label; VIC/VICVIC swaps PB1/PB2 numbers
  map_segment_label() {
    seg="$1"   # HA NA M PB1 PB2 NP PA NS
    case "$subtype" in
      VIC|VICVIC)
        case "$seg" in
          HA)  echo "01-HA" ;;
          NA)  echo "02-NA" ;;
          M)   echo "03-M"  ;;
          PB1) echo "04-PB1" ;;  # VIC swap
          PB2) echo "05-PB2" ;;  # VIC swap
          NP)  echo "06-NP" ;;
          PA)  echo "07-PA" ;;
          NS)  echo "08-NS" ;;
          *)   echo "Unknown" ;;
        esac
        ;;
      *)
        case "$seg" in
          HA)  echo "01-HA" ;;
          NA)  echo "02-NA" ;;
          M)   echo "03-M"  ;;
          PB1) echo "04-PB1" ;;
          PB2) echo "05-PB2" ;;
          NP)  echo "06-NP" ;;
          PA)  echo "07-PA" ;;
          NS)  echo "08-NS" ;;
          *)   echo "Unknown" ;;
        esac
        ;;
    esac
  }

  # flumut/genin combined headers need bare segment token; change M -> MP
  seg_token_for_fg() {
    seg="$1"
    if [ "$seg" = "M" ]; then
      echo "MP"
    else
      echo "$seg"
    fi
  }

  # Detect segment token from filename (case-insensitive).
  # Returns: HA NA M PB1 PB2 NP PA NS ; empty string if none found.
  detect_seg_token() {
    fn="$1"
    u="$(printf "%s" "$fn" | tr '[:lower:]' '[:upper:]')"
    for t in HA NA M PB1 PB2 NP PA NS; do
      # require non-alnum boundaries to avoid partial matches
      pattern="(^|[^A-Z0-9])${t}([^A-Z0-9]|$)"
      if echo "$u" | grep -Eq "$pattern"; then
        echo "$t"
        return
      fi
    done
    echo ""
  }

  # Clean outputs
  : > "${META_ID}.fasta"
  : > "${META_ID}_flumut.fasta"
  : > "${META_ID}_genin.fasta"

  any_written=0

  # Iterate over staged fasta-like files, skip outputs & subtype file
  found_files=0
  for f in *.fa *.fasta; do
    [ -e "$f" ] || continue
    case "$f" in
      "${META_ID}.fasta"|"${META_ID}_flumut.fasta"|"${META_ID}_genin.fasta") continue ;;
      *_subtype.txt) continue ;;
    esac
    found_files=1

    base="$(basename "$f")"
    ext="${base##*.}"
    case "${ext,,}" in
      fa|fasta) ;;  # accept fasta-like
      *) continue ;;
    esac

    seg_tok="$(detect_seg_token "$base")"
    [ -z "$seg_tok" ] && continue

    seg_label="$(map_segment_label "$seg_tok")"
    [ "$seg_label" = "Unknown" ] && continue

    # Headers
    new_header=">${META_ID}|${seg_label}-${subtype}"
    seg_fg="$(seg_token_for_fg "$seg_tok")"
    flumut_header=">${META_ID}_${seg_fg}"
    genin_header=">${META_ID}_${seg_fg}"

    # Per-segment rewritten file (kept in work dir)
    out_seg="${META_ID}_${seg_label}-${subtype}.fa"

    # Rewrite first header line; keep sequence as-is
    awk -v H="$new_header" 'NR==1{print H; next} {print}' "$f" > "$out_seg"

    # Append to combined outputs with their respective headers
    awk -v H="$new_header"    'NR==1{print H; next} {print}' "$f" >> "${META_ID}.fasta"
    awk -v H="$flumut_header" 'NR==1{print H; next} {print}' "$f" >> "${META_ID}_flumut.fasta"
    awk -v H="$genin_header"  'NR==1{print H; next} {print}' "$f" >> "${META_ID}_genin.fasta"

    any_written=1
  done

  if [ "$found_files" -ne 1 ]; then
    echo "No input segment files (*.fa/*.fasta) staged for ${META_ID}" >&2
    exit 1
  fi

  # Sanity: ensure we produced at least one sequence
  if [ "$any_written" -ne 1 ] || [ ! -s "${META_ID}.fasta" ]; then
    echo "No valid segments found for ${META_ID} among staged files." >&2
    exit 1
  fi
  '''
}

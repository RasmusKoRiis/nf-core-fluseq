process SUBTYPEFINDER {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'
    debug   true

    conda "bioconda::blast=2.15.0"
    container 'biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    tuple val(meta), path(ha_fasta), path(na_fasta)
    path(ha_database)
    path(na_database)

    output:
    tuple val(meta), path("*subtype.txt"), emit: subtype
    tuple val(meta), path("*subtype.csv"), emit: subtype_file
    path("*subtype.csv"), emit: subtype_report
    path("*_subtype.status.tsv"), emit: subtype_status
    path("*_subtype.errors.txt"), optional: true, emit: subtype_errors
    tuple val(meta), path("*.tsv"), emit: subtype_file_test
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -uo pipefail

    # ---------------------------
    # Status + error capture files
    # ---------------------------
    status_file="${meta.id}_subtype.status.tsv"
    errors_file="${meta.id}_subtype.errors.txt"

    echo -e "sample\\tstep\\tok\\tmessage" > "\$status_file"
    : > "\$errors_file"

    log_ok()   { echo -e "${meta.id}\\t\$1\\t1\\t\$2" >> "\$status_file"; }
    log_fail() { echo -e "${meta.id}\\t\$1\\t0\\t\$2" >> "\$status_file"; echo "[\$1] \$2" >> "\$errors_file"; }

    # Always create these so Nextflow outputs exist even on failure
    : > "ha_${prefix}.tsv"
    : > "na_${prefix}.tsv"
    : > "ha_${prefix}.blast.err"
    : > "na_${prefix}.blast.err"
    : > "ha_${prefix}_sorted.tsv"
    : > "na_${prefix}_sorted.tsv"
    : > "ha_${prefix}_filtered.tsv"
    : > "na_${prefix}_filtered.tsv"

    subtype_ha=""
    subtype_na=""

    # ---------------------------
    # HA BLAST
    # ---------------------------
    if [[ ! -s "$ha_fasta" ]]; then
        log_fail "HA_input" "HA query fasta is missing or empty: $ha_fasta"
    else
        if blastn -task blastn -query "$ha_fasta" -subject "$ha_database" -outfmt 6 -max_target_seqs 2 \
            > "ha_${prefix}.tsv" 2> "ha_${prefix}.blast.err"; then
            log_ok "HA_blast" "blastn exited 0"
        else
            log_fail "HA_blast" "blastn failed (non-zero exit). See ha_${prefix}.blast.err"
        fi

        if [[ -s "ha_${prefix}.tsv" ]]; then
            sort -k12,12nr -k3,3nr "ha_${prefix}.tsv" | head -n 1 > "ha_${prefix}_sorted.tsv" || true
            head -n 1 "ha_${prefix}_sorted.tsv" > "ha_${prefix}_filtered.tsv" || true
            log_ok "HA_hits" "BLAST produced hits"
        else
            log_fail "HA_hits" "No hits produced for HA (or BLAST failed)."
        fi
    fi

    # HA subtype call (coverage threshold)
    if [[ -s "ha_${prefix}_filtered.tsv" ]]; then
        if awk -F '\\t' '{if(\$4 >= 0.2 * 1701) exit 0; else exit 1}' "ha_${prefix}_filtered.tsv"; then
            subtype_ha=\$(awk -F '\\t' '{split(\$2,a,"_"); print a[2]; exit}' "ha_${prefix}_filtered.tsv")
            log_ok "HA_call" "Subtype HA called: \${subtype_ha}"
        else
            log_fail "HA_call" "Top HA hit below coverage threshold (<20% of 1701bp)"
        fi
    else
        log_fail "HA_call" "No HA filtered hit line available"
    fi

    # ---------------------------
    # NA BLAST
    # ---------------------------
    if [[ ! -s "$na_fasta" ]]; then
        log_fail "NA_input" "NA query fasta is missing or empty: $na_fasta"
    else
        if blastn -task blastn -query "$na_fasta" -subject "$na_database" -outfmt 6 -max_target_seqs 2 \
            > "na_${prefix}.tsv" 2> "na_${prefix}.blast.err"; then
            log_ok "NA_blast" "blastn exited 0"
        else
            log_fail "NA_blast" "blastn failed (non-zero exit). See na_${prefix}.blast.err"
        fi

        if [[ -s "na_${prefix}.tsv" ]]; then
            sort -k12,12nr -k3,3nr "na_${prefix}.tsv" | head -n 1 > "na_${prefix}_sorted.tsv" || true
            head -n 1 "na_${prefix}_sorted.tsv" > "na_${prefix}_filtered.tsv" || true
            log_ok "NA_hits" "BLAST produced hits"
        else
            log_fail "NA_hits" "No hits produced for NA (or BLAST failed)."
        fi
    fi

    # NA subtype call (coverage threshold)
    if [[ -s "na_${prefix}_filtered.tsv" ]]; then
        if awk -F '\\t' '{if(\$4 >= 0.2 * 1410) exit 0; else exit 1}' "na_${prefix}_filtered.tsv"; then
            subtype_na=\$(awk -F '\\t' '{split(\$2,a,"_"); print a[2]; exit}' "na_${prefix}_filtered.tsv")
            log_ok "NA_call" "Subtype NA called: \${subtype_na}"
        else
            log_fail "NA_call" "Top NA hit below coverage threshold (<20% of 1410bp)"
        fi
    else
        log_fail "NA_call" "No NA filtered hit line available"
    fi

    # If no failures were logged, remove the error file (optional cleanliness)
    if ! grep -qP "\\t0\\t" "\$status_file"; then
        rm -f "\$errors_file"
        log_ok "overall" "No errors detected"
    else
        log_fail "overall" "One or more steps failed; see \$errors_file and *.blast.err"
    fi

    # ---------------------------
    # Final subtype outputs
    # ---------------------------
    subtype="${meta.id},\${subtype_ha}\${subtype_na}"
    subtype_txt="\${subtype_ha}\${subtype_na}"

    echo "\$subtype" > "${meta.id}_subtype.csv"
    echo "Sample,Subtype" | cat - "${meta.id}_subtype.csv" > temp && mv temp "${meta.id}_subtype.csv"

    echo "\$subtype_txt" > "${meta.id}_subtype.txt"

    # ---------------------------
    # Versions
    # ---------------------------
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    blastn: $(blastn -version 2>&1 | sed -E 's#^.*blastn: ([^ ]+).*$#\1#')
    END_VERSIONS
    """
}

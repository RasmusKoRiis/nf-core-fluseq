process NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/rasmuskriis/nextclade-python@sha256:86ee1b9a00da7af2c113aaf937da3554cc72c1f954d79970027941eb2cf7ce52'

    input:
    tuple val(meta), path(fasta), path(subtype)
    path datasets

    output:
    tuple val(meta), path("*nextclade.csv"), emit: nextclade_csv, optional: true
    // H5 non-HA segments are intentionally handled by AMINOACIDTRANSLATION.
    // Keep this channel optional so an all-skipped sample remains a valid task.
    tuple val(meta), path("*translation*fasta"), path(subtype), emit: aminoacid_sequence, optional: true
    tuple val(meta), path("*mutation.csv"), emit: nextclade_filtered, optional: true
    path("*summary.csv"), emit: nextclade_summary_rapport, optional: true
    path("*NC_mutation.csv"), emit: nextclade_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    processed=0
    # COVERAGE removes the subtype from filenames but preserves the FASTA header.
    subtype_name=\$(tr -d '\\r\\n' < "${subtype}")
    for fasta_file in ${fasta}; do
        header=\$(head -n 1 "\$fasta_file" | tr -d '\\r')
        segment_subtype=\${header#*|}
        segment=\${segment_subtype%-"\$subtype_name"}
        segment_name=\${segment##*-}
        if [[ "\$header" != '>'*'|'* || "\$segment" == "\$segment_subtype" ]]; then
            echo "Invalid Nextclade FASTA header for ${meta.id} in \$fasta_file: \$header (expected subtype \$subtype_name)" >&2
            exit 1
        fi

        case "\$subtype_name" in
            H1*) dataset_subtype='H1N1' ;;
            H3*) dataset_subtype='H3N2' ;;
            VIC*|BVIC*|B-VIC*) dataset_subtype='VIC' ;;
            H5*) dataset_subtype="\$subtype_name" ;;
            *)
                echo "Unsupported Nextclade subtype for ${meta.id}: \$subtype_name" >&2
                exit 1
                ;;
        esac

        dataset_segment="\$segment_name"
        [[ "\$dataset_segment" == 'MP' ]] && dataset_segment='M'

        # Preserve the established H5 behavior, which performs Nextclade only
        # for HA. Other avian segments are translated by AMINOACIDTRANSLATION.
        if [[ "\$dataset_subtype" == H5* && "\$dataset_segment" != 'HA' ]]; then
            continue
        fi

        dataset_dir=''
        for candidate in \
            "${datasets}/\${dataset_subtype}_\${dataset_segment}" \
            "${datasets}/B_VIC_\${dataset_segment}" \
            "${datasets}/BVIC_\${dataset_segment}"; do
            if [[ -d "\$candidate" && -f "\$candidate/pathogen.json" ]]; then
                dataset_dir="\$candidate"
                break
            fi
        done
        if [[ -z "\$dataset_dir" ]]; then
            echo "No local Nextclade dataset for \${dataset_subtype}_\${dataset_segment} under ${datasets}" >&2
            exit 1
        fi

        output_dir="${meta.id}_\${segment}_nextclade_output"
        nextclade run \
            --input-dataset "\$dataset_dir" \
            --output-all "\$output_dir" \
            "\$fasta_file"
        processed=1

        if compgen -G "\$output_dir/*" > /dev/null; then
            for output_file in "\$output_dir"/*; do
                output_name=\$(basename "\$output_file")
                if [[ "\$output_file" == *.csv ]]; then
                    mv "\$output_file" "${meta.id}_\${segment}_\${output_name}"
                else
                    mv "\$output_file" "${meta.id}_\${output_name}"
                fi
            done
        fi

        nextclade_csv="${meta.id}_\${segment}_nextclade.csv"
        if [[ -f "\$nextclade_csv" ]]; then
            nextclade_converter.py \
                "\$nextclade_csv" \
                ${meta.id} \
                "\$segment" \
                NC
        fi
    done

    if [ "\$processed" -eq 0 ]; then
        echo "No segments required Nextclade processing for ${meta.id}" >&2
    fi

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1 | sed 's/^.*nextclade //; s/ .*\$//')
        dataset_root: ${datasets}
    END_VERSIONS
    """
}

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
    dataset_records=''
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
        dataset_candidates=("${datasets}/\${dataset_subtype}_\${dataset_segment}")

        # Accept the established flat aliases only for Victoria datasets.
        if [[ "\$dataset_subtype" == 'VIC' ]]; then
            dataset_candidates+=(
                "${datasets}/B_VIC_\${dataset_segment}"
                "${datasets}/BVIC_\${dataset_segment}"
                "${datasets}/B-VIC_\${dataset_segment}"
            )
        fi

        # Also accept the hierarchy produced by the official Nextclade dataset
        # collection. HA and NA normally have an additional reference-accession
        # directory (for example vic/ha/KX058884/pathogen.json).
        official_subtype=''
        case "\$dataset_subtype" in
            H1N1) official_subtype='h1n1pdm' ;;
            H3N2) official_subtype='h3n2' ;;
            VIC) official_subtype='vic' ;;
        esac
        if [[ -n "\$official_subtype" ]]; then
            official_segment=\$(printf '%s' "\$dataset_segment" | tr '[:upper:]' '[:lower:]')
            [[ "\$official_segment" == 'm' ]] && official_segment='mp'
            dataset_candidates+=(
                "${datasets}/nextstrain/flu/\${official_subtype}/\${official_segment}"
                "${datasets}/data/nextstrain/flu/\${official_subtype}/\${official_segment}"
                "${datasets}/flu/\${official_subtype}/\${official_segment}"
                "${datasets}/\${official_subtype}/\${official_segment}"
            )
        fi

        # Keep the HA/NA reference choices used by the production workflow.
        # Internal segments use the official subtype/segment dataset shortcut.
        dataset_name=''
        case "\${dataset_subtype}:\${dataset_segment}" in
            H1N1:HA) dataset_name='nextstrain/flu/h1n1pdm/ha/california-7-2009' ;;
            H1N1:NA) dataset_name='nextstrain/flu/h1n1pdm/na/wisconsin-588-2019' ;;
            H3N2:HA) dataset_name='nextstrain/flu/h3n2/ha/CY163680' ;;
            H3N2:NA) dataset_name='nextstrain/flu/h3n2/na/EPI1857215' ;;
            VIC:HA) dataset_name='nextstrain/flu/vic/ha/KX058884' ;;
            VIC:NA) dataset_name='nextstrain/flu/vic/na/CY073894' ;;
            H1N1:*|H3N2:*|VIC:*) dataset_name="nextstrain/flu/\${official_subtype}/\${official_segment}" ;;
            H5*:HA) dataset_name='community/moncla-lab/iav-h5/ha/all-clades' ;;
        esac

        # Prefer a freshly downloaded official dataset. If the task has no
        # network access, retain support for the controlled local bundle.
        if [[ -n "\$dataset_name" ]]; then
            download_dir="${meta.id}_\${segment}_nextclade_dataset"
            echo "Downloading Nextclade dataset \$dataset_name" >&2
            if nextclade dataset get \
                --name "\$dataset_name" \
                --output-dir "\$download_dir"; then
                if [[ -f "\$download_dir/pathogen.json" ]]; then
                    dataset_dir="\$download_dir"
                else
                    echo "Downloaded dataset \$dataset_name contains no pathogen.json; trying local datasets" >&2
                fi
            else
                echo "Could not download \$dataset_name; trying local datasets" >&2
            fi
        fi

        if [[ -z "\$dataset_dir" ]]; then
            for candidate in "\${dataset_candidates[@]}"; do
                if [[ -f "\$candidate/pathogen.json" ]]; then
                    dataset_dir="\$candidate"
                    break
                fi

                # Resolve a single reference-specific child directory without
                # silently choosing between multiple installed reference datasets.
                nested_dataset=''
                for pathogen_file in "\$candidate"/*/pathogen.json; do
                    [[ -f "\$pathogen_file" ]] || continue
                    resolved_dataset=\${pathogen_file%/pathogen.json}
                    if [[ -n "\$nested_dataset" && "\$nested_dataset" != "\$resolved_dataset" ]]; then
                        echo "Multiple local Nextclade datasets found under \$candidate; use the flat <subtype>_<segment> layout to select one" >&2
                        exit 1
                    fi
                    nested_dataset="\$resolved_dataset"
                done
                if [[ -n "\$nested_dataset" ]]; then
                    dataset_dir="\$nested_dataset"
                    break
                fi
            done
        fi
        if [[ -z "\$dataset_dir" ]]; then
            echo "Unable to download or find a local Nextclade dataset for \${dataset_subtype}_\${dataset_segment}" >&2
            [[ -n "\$dataset_name" ]] && echo "Download attempted: \$dataset_name" >&2
            echo "Checked candidate roots:" >&2
            for candidate in "\${dataset_candidates[@]}"; do
                if [[ -d "\$candidate" ]]; then
                    echo "  exists but contains no resolvable pathogen.json: \$candidate" >&2
                else
                    echo "  missing: \$candidate" >&2
                fi
            done
            exit 1
        fi

        output_dir="${meta.id}_\${segment}_nextclade_output"
        nextclade run \
            --input-dataset "\$dataset_dir" \
            --output-all "\$output_dir" \
            "\$fasta_file"
        processed=1
        dataset_records="\${dataset_records}\${dataset_records:+,}\${segment_name}:\${dataset_dir}"

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
        datasets: "\$dataset_records"
    END_VERSIONS
    """
}

process NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'ignore'
    
    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory
    //container logic as needed

    input:
    tuple val(meta), path(fasta), path(subtype)
    path(dataset)
    

    output:
    tuple val(meta), path("*nextclade.csv"), emit: nextclade_csv, optional: true
    tuple val(meta), path("*translation*fasta"), path(subtype), emit: aminoacid_sequence
    tuple val(meta), path("*mutation.csv"), emit: nextclade_filtered, optional: true


    path("*summary.csv"), emit: nextclade_summary_rapport, optional: true
    path("*NC_mutation.csv"), emit: nextclade_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    for fasta_file in ${fasta}; do

        filename=\$(basename \$fasta_file)
        filename_no_ext=\${filename%.*}  
        segment_subtype=\${filename_no_ext#*-} 
        segment=\${segment_subtype%-*}  
        subtype_name=\${segment_subtype#*-} 
        
    
        dataset_sample=${dataset}/"\${subtype_name}_\${segment}"

    # Check for specific subtype and segment combinations
        if [[ "\$subtype_name" == *"H3"* ]]; then
            if [[ "\$segment" == *"HA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/ha/CY163680' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/na/EPI1857215' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            
            elif [[ "\$segment" == *"PB2"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/pb2' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"PB1"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/pb1' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"PA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/pa' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NP"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/np' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NS"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/ns' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            
            elif [[ "\$segment" == *"M"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h3n2/mp' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            else
                echo "Segment \$segment not recognized for subtype H1"
            fi

        elif [[ "\$subtype_name" == *"H1"* ]]; then
            if [[ "\$segment" == *"HA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/ha/california-7-2009' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/na/wisconsin-588-2019' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            
            elif [[ "\$segment" == *"PB2"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pb2' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"PB1"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pb1' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"PA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/pa' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NP"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/np' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NS"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/ns' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            
            elif [[ "\$segment" == *"M"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/h1n1pdm/mp' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            else
                echo "Segment \$segment not recognized for subtype H1"
            fi

        elif [[ "\$subtype_name" == *"VIC"* ]]; then
            if [[ "\$segment" == *"HA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/vic/ha/KX058884' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file

            elif [[ "\$segment" == *"NA"* ]]; then
                nextclade dataset get --name 'nextstrain/flu/vic/na/CY073894' --output-dir "${meta.id}_\${segment}_nextclade_dataset/"

                nextclade run \
                    --input-dataset "${meta.id}_\${segment}_nextclade_dataset/" \
                    --output-all=${meta.id}_\${segment}_nextclade_output/ \
                    \$fasta_file
            else
                echo "Segment \$segment not recognized for subtype VIC"
            fi
        else
            echo "Subtype \$subtype_name not recognized or not handled."
        fi

        # Save output files from Nextclade
        
        if compgen -G "${meta.id}_\${segment}_nextclade_output/*" > /dev/null; then
            for file in ${meta.id}_\${segment}_nextclade_output/*; do
                cat "\$file"
                basename=\$(basename \$file)
                if [[ "\$file" == *.csv ]]; then
                    mv "\$file" ./${meta.id}_\${segment}_\$basename
                else
                    mv "\$file" ./${meta.id}_\$basename
                fi
            done
        fi

        # Convert Nextclade output to mutations, frameshift and glyco files
        type=NC
        nextclade_csv="./${meta.id}_\${segment}_nextclade.csv"
        if [[ -f "\$nextclade_csv" ]]; then
            python3 /project-bin/nextclade_converter.py \
                "\$nextclade_csv" \
                ${meta.id} \
                \${segment} \
                \${type}
        fi    
        subtype=\$(cat ${subtype})

    done


    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS

    """
}

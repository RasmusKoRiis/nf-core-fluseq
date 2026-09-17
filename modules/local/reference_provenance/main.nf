process REFERENCE_PROVENANCE {
    tag 'runtime references'
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

    input:
    path references

    output:
    path 'reference_manifest.tsv', emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    printf 'path\tsha256\n' > reference_manifest.tsv
    find -L ${references} -type f -print \
        | sort \
        | while IFS= read -r reference_file; do
            checksum=\$(sha256sum "\$reference_file" | cut -d ' ' -f 1)
            printf '%s\t%s\n' "\$reference_file" "\$checksum" >> reference_manifest.tsv
        done

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        coreutils: 9.1
    END_VERSIONS
    """
}

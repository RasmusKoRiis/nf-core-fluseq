process REFERENCE_PROVENANCE {
    tag 'runtime references'
    label 'process_single'

    container 'docker.io/library/alpine@sha256:029a752048e32e843bd6defe3841186fb8d19a28dae8ec287f433bb9d6d1ad85'

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
        alpine: 3.20.3
    END_VERSIONS
    """
}

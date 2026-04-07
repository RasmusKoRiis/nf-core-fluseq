process SUBCLADE_NOMENCLATURE_RULES {
    tag 'github'
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'

    output:
    path("subclade_nomenclature"), emit: rules_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    '''
    set -euo pipefail

    python - <<'PY'
import os
import tarfile
import urllib.request

downloads = {
    "H3N2_HA": "https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/archive/refs/heads/main.tar.gz",
    "H1N1pdm_HA": "https://github.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA/archive/refs/heads/main.tar.gz",
    "B-Vic_HA": "https://github.com/influenza-clade-nomenclature/seasonal_B-Vic_HA/archive/refs/heads/main.tar.gz",
}

for profile, url in downloads.items():
    archive_path = f"{profile}.tar.gz"
    with urllib.request.urlopen(url, timeout=120) as response:
        payload = response.read()
    with open(archive_path, "wb") as handle:
        handle.write(payload)

    extracted = 0
    out_root = os.path.join("subclade_nomenclature", profile)
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            parts = member.name.split("/")
            if member.isfile() and len(parts) >= 3 and parts[1] in {"subclades", "clades"} and parts[-1].endswith(".yml"):
                out_dir = os.path.join(out_root, parts[1])
                os.makedirs(out_dir, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    continue
                with open(os.path.join(out_dir, parts[-1]), "wb") as handle:
                    handle.write(source.read())
                extracted += 1
    if extracted == 0:
        raise RuntimeError(f"No clade/subclade YAML files extracted from {url}")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: $(python --version 2>&1)
    END_VERSIONS
    '''
}

process SUBCLADE_NOMENCLATURE {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(fasta), path(subtype), path(coverage_csv)
    path rules_dir

    output:
    tuple val(meta), path("${meta.id}_subclade_nomenclature.csv"), emit: calls
    path("${meta.id}_subclade_nomenclature.csv"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python /project-bin/subclade_nomenclature.py \\
        --sample-id ${meta.id} \\
        --subtype-file ${subtype} \\
        --rules-dir ${rules_dir} \\
        --output ${meta.id}_subclade_nomenclature.csv \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

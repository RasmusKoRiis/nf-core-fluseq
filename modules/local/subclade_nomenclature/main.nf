process SUBCLADE_NOMENCLATURE_RULES {
    tag 'github'
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'

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
    "H3N2_HA": "https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/archive/740628ac709ac63cdee2a9c5450ae1172d6b5304.tar.gz",
    "H1N1pdm_HA": "https://github.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA/archive/ac496f507222e5f0c527a10277fd976191cc291c.tar.gz",
    "B-Vic_HA": "https://github.com/influenza-clade-nomenclature/seasonal_B-Vic_HA/archive/34c60b4f93067e72ad3e2b6b2f42e052cd162a8e.tar.gz",
}
references = {
    "H3N2_HA": "CY163680.1",
    "H1N1pdm_HA": "CY121680.1",
    "B-Vic_HA": "KX058884.1",
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

    ref_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={references[profile]}&rettype=fasta&retmode=text"
    )
    with urllib.request.urlopen(ref_url, timeout=120) as response:
        reference = response.read()
    if not reference.startswith(b">"):
        raise RuntimeError(f"Unexpected reference FASTA downloaded from {ref_url}")
    with open(os.path.join(out_root, "reference.fasta"), "wb") as handle:
        handle.write(reference)
PY

cat > versions.yml <<END_VERSIONS
"SUBCLADE_NOMENCLATURE_RULES":
    python: $(python --version 2>&1)
    h3n2_rules: 740628ac709ac63cdee2a9c5450ae1172d6b5304
    h1n1pdm_rules: ac496f507222e5f0c527a10277fd976191cc291c
    b_vic_rules: 34c60b4f93067e72ad3e2b6b2f42e052cd162a8e
END_VERSIONS
    '''
}

process SUBCLADE_NOMENCLATURE {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/rasmuskriis/blast_python_pandas@sha256:fd100d56162d663949f23a0c26bee52a6d4b0da66235ce0aa53407353185b66a'
    input:
    tuple val(meta), path(fasta), path(subtype), path(coverage_csv)
    path rules_dir
    path caller_script
    path characterisation_script
    path characterisation_guidelines

    output:
    tuple val(meta), path("${meta.id}_subclade_nomenclature.csv"), emit: calls
    path("${meta.id}_subclade_nomenclature.csv"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${caller_script} \\
        --sample-id ${meta.id} \\
        --subtype-file ${subtype} \\
        --rules-dir ${rules_dir} \\
        --output ${meta.id}_subclade_nomenclature_raw.csv \\
        ${fasta}

    python ${characterisation_script} \\
        --input ${meta.id}_subclade_nomenclature_raw.csv \\
        --subtype-file ${subtype} \\
        --guidelines-dir ${characterisation_guidelines} \\
        --output ${meta.id}_subclade_nomenclature.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}

"""Check Coverage -> Nextclade naming with synthetic tools and local datasets."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
@pytest.mark.parametrize(
    ("sample", "subtype", "segment", "dataset"),
    [
        ("2791130-INFA", "H3N2", "01-HA", "H3N2_HA"),
        ("sample", "H1N1", "02-NA", "H1N1_NA"),
        ("2786140-INFB", "VICVIC", "03-MP", "VIC_M"),
        ("sample_with_underscores-INFA", "H3N2", "05-PB2", "H3N2_PB2"),
    ],
)
def test_nextclade_uses_subtype_file_after_coverage(tmp_path, sample, subtype, segment, dataset):
    fasta_name = f"{sample}_{segment}-{subtype}.fa"
    (tmp_path / fasta_name).write_text(f">{sample}|{segment}-{subtype}\nACGT\n")
    (tmp_path / "subtype.txt").write_text(subtype + "\n")
    dataset_dir = tmp_path / "datasets" / dataset
    dataset_dir.mkdir(parents=True)
    (dataset_dir / "pathogen.json").write_text("{}\n")
    (tmp_path / "nextflow.config").write_text("process.errorStrategy = 'terminate'\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    scripts = {
        "coverage_finder.py": (
            '#!/bin/bash\nset -euo pipefail\n'
            'printf "Sample,Coverage\\n%s,100\\n" "$3" > "$2"\n'
            'printf "100\\n" > "${3}_${4}_coverage.txt"\n'
        ),
        "nextclade": (
            '#!/bin/bash\nset -euo pipefail\n'
            'if [[ "$1" == --version ]]; then echo "nextclade 3.0.0"; exit 0; fi\n'
            '[[ "$1" == run && "$2" == --input-dataset && "$4" == --output-all ]]\n'
            f'[[ "$3" == "datasets/{dataset}" ]]\n'
            '[[ -f "$6" ]]\n'
            'mkdir -p "$5"\n'
            'printf "%s\\n" "$3" > "$5/nextclade.csv"\n'
            'printf ">sample\\nACDE\\n" > "$5/nextclade.cds_translation.HA1.fasta"\n'
        ),
        "nextclade_converter.py": (
            '#!/bin/bash\nset -euo pipefail\n'
            f'[[ "$2" == "{sample}" && "$3" == "{segment}" && "$4" == NC ]]\n'
            '[[ -s "$1" ]]\n'
        ),
    }
    for name, script in scripts.items():
        path = bin_dir / name
        path.write_text(script)
        path.chmod(0o755)
    (tmp_path / "main.nf").write_text(
        f"include {{ COVERAGE }} from '{ROOT / 'modules/local/coverage/main'}'\n"
        f"include {{ NEXTCLADE }} from '{ROOT / 'modules/local/nextclade/main'}'\n"
        "workflow {\n"
        f"    COVERAGE(Channel.of(tuple([id: '{sample}'], file('{fasta_name}'), file('subtype.txt'))), 80)\n"
        "    inputs = COVERAGE.out.filtered_fasta.map { meta, fasta, subtype, report -> tuple(meta, fasta, subtype) }\n"
        "    NEXTCLADE(inputs, Channel.value(file('datasets')))\n"
        "    NEXTCLADE.out.nextclade_csv.toList().view { rows ->\n"
        "        assert rows.size() == 1\n"
        f"        assert rows[0][1].text.trim() == 'datasets/{dataset}'\n"
        f"        assert rows[0][1].name == '{sample}_{segment}_nextclade.csv'\n"
        "        'NEXTCLADE_COMPLETED'\n"
        "    }\n"
        "}\n"
    )
    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-ansi-log", "false"],
        cwd=tmp_path,
        env={**os.environ, "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_SYNTAX_PARSER": "v1"},
        capture_output=True,
        text=True,
        timeout=90,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "NEXTCLADE_COMPLETED" in output, output

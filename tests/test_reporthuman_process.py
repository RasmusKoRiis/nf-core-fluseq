"""Exercise REPORTHUMAN's embedded CSV-to-TSV conversion."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
def test_reporthuman_embedded_python_preserves_escaped_delimiters(tmp_path):
    csv_inputs = [
        "subtype.csv",
        "coverage.csv",
        "mutation_human.csv",
        "mutation_inhibition.csv",
        "lookup.csv",
        "nextclade_summary.csv",
        "nextclade_sample.csv",
        "mutation_vaccine.csv",
        "irma_depth.csv",
        "reassortment.csv",
        "subclade.csv",
    ]
    for name in csv_inputs:
        (tmp_path / name).write_text("Sample\nSAMPLE01\n", encoding="utf-8")
    (tmp_path / "filtered.fasta").write_text(">SAMPLE01\nACGT\n", encoding="utf-8")
    (tmp_path / "samplesheet.csv").write_text(
        "PCR-PlatePosition,SequenceID,Barcode,KonsCt\n" 'A1,SAMPLE01,barcode01,"24,56"\n',
        encoding="utf-8",
    )
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\n" "process.errorStrategy = 'terminate'\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    report = bin_dir / "report.py"
    report.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        "expected=$'A1\\tSAMPLE01\\tbarcode01\\t24,56'\n"
        "actual=$(sed -n '2p' \"$1\")\n"
        '[[ "$actual" == "$expected" ]]\n'
        "printf 'Sample,Subtype,Coverage-HA,Coverage-NA\\nSAMPLE01,H3N2,99,99\\n' > merged_report.csv\n",
        encoding="utf-8",
    )
    report.chmod(0o755)
    qc = bin_dir / "report_QC_calculation.py"
    qc.write_text(
        "#!/bin/bash\n" "set -euo pipefail\n" "input=$1\n" "[[ \"$2\" == '-o' ]]\n" 'cp "$input" "$3"\n',
        encoding="utf-8",
    )
    qc.chmod(0o755)

    args = ",\n        ".join(
        [
            "Channel.value(file('subtype.csv'))",
            "Channel.value(file('coverage.csv'))",
            "Channel.value(file('mutation_human.csv'))",
            "Channel.value(file('mutation_inhibition.csv'))",
            "Channel.value(file('lookup.csv'))",
            "Channel.value(file('nextclade_summary.csv'))",
            "Channel.value(file('nextclade_sample.csv'))",
            "Channel.value(file('mutation_vaccine.csv'))",
            "'TEST'",
            "'test-version'",
            "Channel.value(file('filtered.fasta'))",
            "Channel.value(file('irma_depth.csv'))",
            "'test-instrument'",
            "Channel.value(file('samplesheet.csv'))",
            "Channel.value(file('reassortment.csv'))",
            "Channel.value(file('subclade.csv'))",
        ]
    )
    (tmp_path / "main.nf").write_text(
        f"include {{ REPORTHUMAN }} from '{ROOT / 'modules/local/reporthuman/main'}'\n"
        "workflow {\n"
        "    REPORTHUMAN(\n"
        f"        {args}\n"
        "    )\n"
        "    REPORTHUMAN.out.report.view { report ->\n"
        "        assert report.name == 'TEST.csv'\n"
        "        'REPORTHUMAN_COMPLETED'\n"
        "    }\n"
        "}\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-ansi-log", "false"],
        cwd=tmp_path,
        env={
            **os.environ,
            "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
            "NXF_OFFLINE": "true",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_SYNTAX_PARSER": "v1",
        },
        capture_output=True,
        text=True,
        timeout=90,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "REPORTHUMAN_COMPLETED" in output

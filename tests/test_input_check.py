"""Exercise routine sample-sheet validation through the Nextflow harness."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
HARNESS = ROOT / "tests" / "input_check" / "main.nf"
CONFIG = ROOT / "tests" / "input_check" / "nextflow.config"
ROUTINE_SHEET = ROOT / "tests" / "data" / "input_check" / "routine.csv"
FASTQ_ROOT = ROOT / "tests" / "data" / "input_check" / "fastq"


def run_input_check(tmp_path, samples_dir):
    command = [
        NEXTFLOW,
        "-C",
        str(CONFIG),
        "run",
        str(HARNESS),
        "--input",
        str(ROUTINE_SHEET),
        "--samplesDir",
        str(samples_dir),
        "-ansi-log",
        "false",
    ]
    return subprocess.run(
        command,
        cwd=tmp_path,
        env={
            **os.environ,
            "NXF_OFFLINE": "true",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_SYNTAX_PARSER": "v1",
        },
        capture_output=True,
        text=True,
        timeout=90,
    )


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
def test_routine_sheet_accepts_legacy_samples_dir(tmp_path):
    result = run_input_check(tmp_path, FASTQ_ROOT)
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    assert "Validated input contract:" in output


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
def test_missing_fastq_root_reports_the_sample_and_path(tmp_path):
    missing_root = tmp_path / "missing-fastq-root"
    result = run_input_check(tmp_path, missing_root)
    output = result.stdout + result.stderr

    assert result.returncode != 0
    assert "No .fastq.gz or .fq.gz files were found for sample SAMPLE01" in output
    assert str(missing_root / "barcode01") in output
    assert "Unexpected error [InvocationTargetException]" not in output

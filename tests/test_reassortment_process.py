"""Exercise FASTA staging and N-content reporting through the Nextflow module."""

import csv
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
BLASTN = shutil.which("blastn")


@pytest.mark.skipif(NEXTFLOW is None or BLASTN is None, reason="Nextflow and BLASTN are required")
def test_reassortment_module_reports_full_query_n_content(tmp_path):
    rng = random.Random(42)
    queries = []
    references = []
    for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
        # Arbitrary synthetic bases, with Ns only outside the exact-match region.
        sequence = "".join(rng.choices("ACGT", k=200))
        queries.append(f">SAMPLE01_{segment}\n{sequence}{'N' * 50}\n")
        references.append(f">HUMAN-SEASONAL|H3N2|A/Synthetic/1/2026|{segment}|TEST_{segment}\n{sequence}\n")
    (tmp_path / "query.fasta").write_text("".join(queries), encoding="utf-8")
    (tmp_path / "references.fasta").write_text("".join(references), encoding="utf-8")
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\n"
        "process.errorStrategy = 'terminate'\n"
        "process.publishDir = [path: 'published', mode: 'copy']\n"
        "docker.enabled = false\n",
        encoding="utf-8",
    )
    (tmp_path / "main.nf").write_text(
        f"include {{ REASSORTMENT }} from '{ROOT / 'modules/local/reassortment/main'}'\n"
        "workflow {\n"
        "    REASSORTMENT(\n"
        "        Channel.of(tuple([id: 'SAMPLE01'], file('query.fasta'))),\n"
        "        Channel.value(file('references.fasta'))\n"
        "    )\n"
        "}\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-ansi-log", "false"],
        cwd=tmp_path,
        env={
            **os.environ,
            "PATH": os.pathsep.join([str(ROOT / "bin"), str(Path(sys.executable).parent), os.environ["PATH"]]),
            "NXF_OFFLINE": "true",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_SYNTAX_PARSER": "v1",
        },
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with (tmp_path / "published/SAMPLE01_reassortment_summary.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))

    assert row["Sample"] == "SAMPLE01"
    for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]:
        assert row[segment] == "HUMAN-SEASONAL|H3N2|A/Synthetic/1/2026(M:100.0%/N:20.0%)"
    assert row["Conclusion"] == "CONSISTENT - all eight segments match HUMAN-SEASONAL H3N2 references"

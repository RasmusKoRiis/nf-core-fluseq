"""Exercise FASTA header parsing, sample grouping, and file staging with synthetic data."""

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
@pytest.mark.parametrize("workflow", ["avian-fasta", "human-fasta"])
def test_full_fasta_headers_survive_grouping_and_staging(tmp_path, workflow):
    # Shared first words reproduce truncation; the third name checks legacy input.
    samples = ["A/Larus test/001/2026", "A/Larus other/002/2026", "A/Synthetic/003/2026"]
    records = {
        f"{sample}|{segment}": sequence
        for sample in samples
        for segment, sequence in [("HA", "ACGT"), ("NA", "TGCA")]
    }
    fasta = tmp_path / "input.fasta"
    fasta.write_text("".join(f">{header}\n{seq}\n" for header, seq in records.items()))
    expected = {
        sample: sorted(
            re.sub(r"[^A-Za-z0-9._-]+", "_", f"{sample}|{segment}") + ".fasta"
            for segment in ["HA", "NA"]
        )
        for sample in samples
    }
    (tmp_path / "expected.json").write_text(json.dumps(expected))

    # Execute the production channel setup, so this test catches use of rec.id.
    source = (ROOT / "workflows" / f"{workflow}.nf").read_text()
    preparation = source.split("/* 1) Split multi-FASTA", 1)[1].split("*/", 1)[1]
    preparation = preparation.split("/* 4)", 1)[0]
    (tmp_path / "main.nf").write_text(
        f"include {{ EMIT_FASTA_RECORD }} from '{ROOT / 'modules/local/emit_fasta_record/main'}'\n"
        f"include {{ WRITE_ID_MAP }} from '{ROOT / 'modules/local/write_id_map/main'}'\n"
        '''
        process STAGE_RECORDS {
            input:
            tuple val(meta), path(records)
            output:
            tuple val(meta), path('bundle.txt')
            script:
            """
            cat ${records} > bundle.txt
            """
        }

        workflow {
        '''
        + preparation
        + """
            def expected = new groovy.json.JsonSlurper().parse(file('expected.json').toFile())
            ch_sample_info.map { meta, files ->
                assert expected.containsKey(meta.orig): "Unexpected sample: ${meta.orig}"
                assert files.collect { it.name }.sort() == expected[meta.orig]
                files.each { f ->
                    assert f.readLines()[0] == ">${meta.id}"
                }
                tuple(meta, files)
            }.set { checked_records }
            STAGE_RECORDS(checked_records)
            STAGE_RECORDS.out.toList().view { rows ->
                assert rows.size() == expected.size()
                assert rows.collect { it[0].id }.unique().size() == expected.size()
                assert rows.collect { it[0].orig }.sort() == expected.keySet().sort()
                rows.each { meta, bundle ->
                    assert bundle.readLines().findAll { !it.startsWith('>') }.sort() == ['ACGT', 'TGCA']
                }
                'FASTA_HEADERS_OK'
            }
        }
        """
    )
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\n"
        "process.errorStrategy = 'terminate'\n"
        "docker.enabled = false\n"
        "conda.enabled = false\n"
    )
    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-lib", str(ROOT / "lib"),
         "--fasta", str(fasta), "-ansi-log", "false"],
        cwd=tmp_path,
        env={**os.environ, "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_SYNTAX_PARSER": "v1"},
        capture_output=True,
        text=True,
        timeout=90,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "FASTA_HEADERS_OK" in output, output

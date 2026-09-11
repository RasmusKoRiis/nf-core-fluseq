"""Check Nextclade translation names at the MUTATIONHUMAN boundary."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
def test_nextclade_h3_feature_names_map_to_mutation_reference_names(tmp_path):
    feature_names = ("HA", "HA2", "M", "NA2", "NP", "NS", "PA", "PB1", "PB2", "SIG")
    for feature in feature_names:
        (tmp_path / f"sample_nextclade.cds_translation.A_H3_{feature}.fasta").write_text(
            ">sample|translated\nACDE\n",
            encoding="utf-8",
        )
    (tmp_path / "subtype.txt").write_text("H3N2\n", encoding="utf-8")
    (tmp_path / "references").mkdir()
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\n" "process.errorStrategy = 'terminate'\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_mutation_finder = bin_dir / "mutation_finder.py"
    fake_mutation_finder.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        "segment=$3\n"
        "output=$5\n"
        "mutation_type=$6\n"
        "printf 'Sample,%s Differences %s\\nsample,No mutations found\\n' "
        '"$segment" "$mutation_type" > "$output"\n'
        'cp "$output" "${output%.csv}_report.csv"\n',
        encoding="utf-8",
    )
    fake_mutation_finder.chmod(0o755)

    fasta_inputs = ", ".join(
        f"file('sample_nextclade.cds_translation.A_H3_{feature}.fasta')" for feature in feature_names
    )
    expected_human = {
        "sample_HA1_human_mutation.csv",
        "sample_HA2_human_mutation.csv",
        "sample_M1_human_mutation.csv",
        "sample_NA_human_mutation.csv",
        "sample_NP_human_mutation.csv",
        "sample_NS1_human_mutation.csv",
        "sample_PA_human_mutation.csv",
        "sample_PB1_human_mutation.csv",
        "sample_PB2_human_mutation.csv",
        "sample_SigPep_human_mutation.csv",
    }
    expected_literal = ", ".join(f"'{name}'" for name in sorted(expected_human))
    (tmp_path / "main.nf").write_text(
        f"include {{ MUTATIONHUMAN }} from '{ROOT / 'modules/local/mutationhuman/main'}'\n"
        "workflow {\n"
        f"    translations = Channel.of(tuple([id: 'sample'], [{fasta_inputs}], file('subtype.txt')))\n"
        "    MUTATIONHUMAN(translations, Channel.value(file('references')))\n"
        "    MUTATIONHUMAN.out.human_mutation.view { meta, outputs, subtype ->\n"
        "        def names = outputs instanceof List ? outputs*.name : [outputs.name]\n"
        f"        assert names.toSet() == [{expected_literal}].toSet()\n"
        "        'MUTATIONHUMAN_COMPLETED'\n"
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
    assert "MUTATIONHUMAN_COMPLETED" in output

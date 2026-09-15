"""Regression tests for lookup input routing and filename handling."""

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.parametrize("script_name", ["table_lookup.py", "table_lookup_mammalian.py"])
def test_lookup_scripts_accept_header_only_mutation_tables(tmp_path, script_name):
    mutations = tmp_path / "header_only.csv"
    mutations.write_text("Sample,M2 inhibition full amino acid list\n", encoding="utf-8")
    workbook = tmp_path / "lookup.xlsx"
    pd.DataFrame(
        [{"segment": "M2", "subtype": "H5N1", "mutation": "S31N"}]
    ).to_excel(workbook, index=False)
    output = tmp_path / "lookup.csv"

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "bin" / script_name),
            str(mutations),
            str(output),
            str(workbook),
            "M2",
            "H5N1",
            "sample_with_underscores",
            "inhibition",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows == [
        {
            "M2 inhibition mutations": "No matching mutations found",
            "Sample": "sample_with_underscores",
        }
    ]


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
def test_lookup_modules_parse_segments_after_underscore_sample_ids(tmp_path):
    sample = "sample_with_underscores"
    inhibition_inputs = [
        f"{sample}_M2_inhibtion_mutation_full_mutation_list.csv",
        f"{sample}_PA_inhibtion_mutation.csv",
    ]
    mammalian_inputs = [
        f"{sample}_PB2_mamailian_mutation_full_mutation_list.csv",
        f"{sample}_NA_nextclade_lookup_mutations.csv",
    ]
    for input_name in inhibition_inputs + mammalian_inputs:
        (tmp_path / input_name).write_text("Sample,mutations\n", encoding="utf-8")
    (tmp_path / "subtype.txt").write_text("H5N1\n", encoding="utf-8")
    (tmp_path / "lookup.xlsx").write_text("synthetic lookup", encoding="utf-8")
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\nprocess.errorStrategy = 'terminate'\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_scripts = {
        "table_lookup.py": (
            "#!/bin/bash\n"
            "set -euo pipefail\n"
            'case "$1" in *_M2_*) expected=M2 ;; *_PA_*) expected=PA ;; *) exit 1 ;; esac\n'
            '[[ "$4" == "$expected" ]]\n'
            f'[[ "$6" == "{sample}" ]]\n'
            'printf "Sample,%s inhibtion mutations\\n%s,ok\\n" "$4" "$6" > "$2"\n'
        ),
        "table_lookup_mammalian.py": (
            "#!/bin/bash\n"
            "set -euo pipefail\n"
            'case "$1" in *_PB2_*) expected=PB2 ;; *_NA_*) expected=NA ;; *) exit 1 ;; esac\n'
            '[[ "$4" == "$expected" ]]\n'
            f'[[ "$6" == "{sample}" ]]\n'
            'printf "Sample,%s mammalian mutations\\n%s,ok\\n" "$4" "$6" > "$2"\n'
        ),
    }
    for name, contents in fake_scripts.items():
        script = bin_dir / name
        script.write_text(contents, encoding="utf-8")
        script.chmod(0o755)

    (tmp_path / "main.nf").write_text(
        f"include {{ TABLELOOKUP }} from '{ROOT / 'modules/local/tablelookup/main'}'\n"
        f"include {{ TABLELOOKUP_MAMMALIAN }} from '{ROOT / 'modules/local/tablelookup_mammalian/main'}'\n"
        "workflow {\n"
        f"    inhibition = Channel.of(tuple([id: '{sample}'], file('{inhibition_inputs[0]}'), file('subtype.txt')), tuple([id: '{sample}'], file('{inhibition_inputs[1]}'), file('subtype.txt')))\n"
        f"    mammalian = Channel.of(tuple([id: '{sample}'], file('{mammalian_inputs[0]}'), file('subtype.txt')), tuple([id: '{sample}'], file('{mammalian_inputs[1]}'), file('subtype.txt')))\n"
        "    lookup = Channel.value(file('lookup.xlsx'))\n"
        "    TABLELOOKUP(inhibition, lookup)\n"
        "    TABLELOOKUP_MAMMALIAN(mammalian, lookup)\n"
        "    TABLELOOKUP.out.inhibtion_mutations.map { meta, output -> output.name }.collect().view { outputs ->\n"
        f"        assert outputs.toSet() == ['{sample}_M2_inhibtion.csv', '{sample}_PA_inhibtion.csv'].toSet()\n"
        "        'INHIBITION_LOOKUP_COMPLETED'\n"
        "    }\n"
        "    TABLELOOKUP_MAMMALIAN.out.mammalian_mutations.map { meta, output -> output.name }.collect().view { outputs ->\n"
        f"        assert outputs.toSet() == ['{sample}_PB2_mammalian.csv', '{sample}_NA_mammalian.csv'].toSet()\n"
        "        'MAMMALIAN_LOOKUP_COMPLETED'\n"
        "    }\n"
        "}\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-ansi-log", "false"],
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
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "INHIBITION_LOOKUP_COMPLETED" in output
    assert "MAMMALIAN_LOOKUP_COMPLETED" in output


def test_avian_fasta_splits_mutation_lists_before_lookup_filters():
    source = (ROOT / "workflows" / "avian-fasta.nf").read_text(encoding="utf-8")
    lookup_handoff = source.split("def ch_full_mutation_files", 1)[1].split("TABLELOOKUP  (", 1)[0]

    assert ".full_mutation_list.flatMap" in source
    assert "tuple(meta, mutation_file, subtype)" in lookup_handoff
    assert lookup_handoff.count("ch_full_mutation_files.filter") == 2

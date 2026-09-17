"""Exercise file staging and output collection using synthetic tool outputs."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
@pytest.mark.parametrize("segments", [("FIRST",), ("FIRST", "SECOND"), ("MP",)])
def test_translation_handles_coverage_names_and_preserves_segment_outputs(tmp_path, segments):
    sample = "sample_with_underscores-INFA"
    inputs = []
    for index, segment in enumerate(segments, 1):
        name = f"{sample}_{index:02d}-{segment}-TEST.fa"
        (tmp_path / name).write_text(f">{sample}|{index:02d}-{segment}-TEST\nACGT\n")
        inputs.append(f"file('{name}')")
        dataset_segment = "M" if segment == "MP" else segment
        (tmp_path / "datasets" / f"TEST_{dataset_segment}").mkdir(parents=True)
    (tmp_path / "subtype.txt").write_bytes(b"TEST\r\n")
    (tmp_path / "nextflow.config").write_text(
        "process.executor = 'local'\n"
        "process { withName: AMINOACIDTRANSLATION { errorStrategy = 'terminate' } }\n"
    )
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    scripts = {
        "coverage_finder.py": (
            "#!/bin/bash\nset -euo pipefail\n"
            'printf "Sample,Coverage\\n%s,100\\n" "$3" > "$2"\n'
            'printf "100\\n" > "${3}_${4}_coverage.txt"\n'
        ),
        "nextclade": (
            "#!/bin/bash\nset -euo pipefail\n"
            'if [[ "$1" == --version ]]; then echo "nextclade 3.0.0"; exit 0; fi\n'
            '[[ "$1" == run && "$2" == --input-dataset && "$4" == --output-all ]]\n'
            '[[ -d "$3" && -f "$6" ]]\n'
            'dataset_segment="${3##*_}"\n'
            'segment="$dataset_segment"\n'
            '[[ "$segment" == M ]] && segment=MP\n'
            'mkdir -p "$5/nested"\n'
            'printf "%s\\n" "$segment" > "$5/nextclade.csv"\n'
            'printf "%s\\n" "$segment" > "$5/nextclade.json"\n'
            'printf "%s\\n" "$3" > "$5/dataset_path.txt"\n'
            # Identical basenames must survive both loop iterations.
            'printf ">synthetic\\nACDE\\n" > "$5/nextclade.cds_translation.shared.fasta"\n'
        ),
        "csv_conversion_nextclade.py": (
            "#!/bin/bash\nset -euo pipefail\n"
            'segment=$(cat "$1")\n'
            '[[ "$1" == "${2}_${segment}_nextclade.csv" ]]\n'
            'cp "$1" "${2}_${segment}_nextclade_mutations.csv"\n'
            'cp "$1" "${2}_${segment}_nextclade_lookup_mutations.csv"\n'
        ),
    }
    for name, script in scripts.items():
        path = bin_dir / name
        path.write_text(script)
        path.chmod(0o755)
    expected = ", ".join(f"'{segment}'" for segment in segments)
    (tmp_path / "main.nf").write_text(
        f"include {{ COVERAGE }} from '{ROOT / 'modules/local/coverage/main'}'\n"
        f"include {{ AMINOACIDTRANSLATION }} from '{ROOT / 'modules/local/aminoacidtranslation/main'}'\n"
        "workflow {\n"
        f"    COVERAGE(Channel.of(tuple([id: '{sample}'], [{', '.join(inputs)}], file('subtype.txt'))), 80)\n"
        "    inputs = COVERAGE.out.filtered_fasta.map { meta, fasta, subtype, report -> tuple(meta, fasta, subtype) }\n"
        "    AMINOACIDTRANSLATION(inputs, Channel.value(file('datasets')))\n"
        "    AMINOACIDTRANSLATION.out.nextclade_csv.view { outputs ->\n"
        "        def files = outputs instanceof List ? outputs : [outputs]\n"
        f"        assert files.collect {{ it.text.trim() }}.toSet() == [{expected}].toSet()\n"
        "        'CSV_COMPLETED'\n"
        "    }\n"
        "    AMINOACIDTRANSLATION.out.mutation_lookup_csv.view { meta, outputs, subtype ->\n"
        "        def files = outputs instanceof List ? outputs : [outputs]\n"
        f"        assert files.size() == {len(segments)}\n"
        "        assert subtype.text.trim() == 'TEST'\n"
        "        'LOOKUP_COMPLETED'\n"
        "    }\n"
        "    AMINOACIDTRANSLATION.out.aminoacid_sequence.view { meta, outputs, subtype ->\n"
        "        def files = outputs instanceof List ? outputs : [outputs]\n"
        f"        assert files.size() == {len(segments)}\n"
        "        assert files.every { it.text.startsWith('>synthetic') }\n"
        "        assert subtype.text.trim() == 'TEST'\n"
        "        'TRANSLATION_COMPLETED'\n"
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
    for marker in ("CSV_COMPLETED", "LOOKUP_COMPLETED", "TRANSLATION_COMPLETED"):
        assert marker in output, output
    assert "Input tuple does not match" not in output, output
    task_scripts = list((tmp_path / "work").glob("*/*/.command.sh"))
    task_dir = next(path.parent for path in task_scripts if "nextclade run" in path.read_text())
    for segment in segments:
        assert (task_dir / f"{sample}_{segment}_nextclade.json").read_text().strip() == segment
        dataset_segment = "M" if segment == "MP" else segment
        assert (task_dir / f"{sample}_{segment}_dataset_path.txt").read_text().strip().endswith(
            f"TEST_{dataset_segment}"
        )

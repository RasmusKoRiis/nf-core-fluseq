"""Exercise the human workflow's MultiQC inputs without containers or biological data."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow is not installed")
@pytest.mark.parametrize(
    "available_versions", [(), ("SUBTYPEFINDER",), ("NEXTCLADE", "SUBTYPEFINDER", "MUTATIONHUMAN")]
)
def test_multiqc_with_missing_upstream_versions(tmp_path, available_versions):
    source = (ROOT / "workflows/human.nf").read_text()
    # Execute the actual input assembly so this test detects regressions in the workflow.
    assembly = source.split("    ch_multiqc_files = Channel.empty()", 1)[1].split("    MULTIQC (", 1)[0]
    assembly = "    ch_multiqc_files = Channel.empty()" + assembly
    declarations = []
    for name in ("NEXTCLADE", "SUBTYPEFINDER", "MUTATIONHUMAN"):
        channel = "Channel.empty()"
        if name in available_versions:
            filename = f"{name}_versions.yml"
            (tmp_path / filename).write_text("test: version\n")
            channel = f"Channel.value(file('{filename}'))"
        declarations.append(f"    def {name} = [out: [versions: {channel}]]")

    (tmp_path / "fastqc.zip").write_text("synthetic FastQC input\n")
    (tmp_path / "software_versions_mqc.yml").write_text("test: version\n")
    (tmp_path / "nextflow.config").write_text("process.errorStrategy = 'ignore'\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_multiqc = bin_dir / "multiqc"
    fake_multiqc.write_text('#!/bin/sh\n[ "$1" = "--version" ] || exit 1\necho "multiqc, version 1.14"\n')
    fake_multiqc.chmod(0o755)
    (tmp_path / "main.nf").write_text(
        f"include {{ MULTIQC }} from '{ROOT / 'modules/nf-core/multiqc/main'}'\n"
        "workflow {\n"
        + "\n".join(declarations)
        + "\n    def CUSTOM_DUMPSOFTWAREVERSIONS = [out: [mqc_yml: Channel.value(file('software_versions_mqc.yml'))]]\n"
        "    def FASTQC = [out: [zip: Channel.of(tuple([id: 'sample'], file('fastqc.zip')))]]\n"
        "    ch_workflow_summary = Channel.value('workflow: test')\n"
        "    ch_methods_description = Channel.value('methods: test')\n"
        + assembly
        + "    inputs = ch_multiqc_files.collect().map { files ->\n"
        "        assert files.every { it instanceof java.nio.file.Path }\n"
        f"        assert files.size() == {4 + len(available_versions)}\n"
        "        files\n"
        "    }\n"
        "    MULTIQC(inputs, [], [], [])\n"
        "    MULTIQC.out.report.toList().view { reports ->\n"
        "        assert reports.size() == 1\n"
        "        assert reports[0].name == 'multiqc_report.html'\n"
        "        'MULTIQC_COMPLETED'\n"
        "    }\n"
        "}\n"
    )
    result = subprocess.run(
        [NEXTFLOW, "-C", "nextflow.config", "run", "main.nf", "-stub-run", "-ansi-log", "false"],
        cwd=tmp_path,
        env={**os.environ, "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_SYNTAX_PARSER": "v1"},
        capture_output=True,
        text=True,
        timeout=90,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "MULTIQC_COMPLETED" in output, output

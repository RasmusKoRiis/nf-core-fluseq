"""Final FASTA exports survive staging and are published beside the reports."""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow is not installed')


@pytest.mark.parametrize('process,module,port,destination', [
    ('REPORTHUMANFASTA', 'reporthumanfasta', 'filtered_fasta', 'reporthuman'),
    ('DRUG_RESISTANCE_REPORT', 'drug_resistance_report', 'final_fasta', 'report'),
])
@pytest.mark.parametrize('empty', [False, True])
def test_final_fasta_export(tmp_path, process, module, port, destination, empty):
    csv_names = ['subtype', 'coverage', 'mutation', 'resistance', 'lookup',
                 'summary', 'sample', 'vaccine', 'reassortment', 'subclade']
    for name in csv_names:
        (tmp_path / f'{name}.csv').write_text('Sample,Subtype\nUID1,H1N1\n')
    (tmp_path / 'id_map.tsv').write_text('SampleID\tOriginalName\nUID1\tSynthetic\n')
    fastas = []
    expected = ''
    if not empty:
        for index, segment in enumerate(['HA', 'NA'], 1):
            directory = tmp_path / str(index)
            directory.mkdir()
            sequence = f'>UID1|0{index}-{segment}-H1N1\nACGT\n'
            (directory / 'same.fasta').write_text(sequence)
            fastas.append(f"file('{index}/same.fasta')")
            expected += sequence
    fasta_input = 'Channel.value([' + ', '.join(fastas) + '])'
    inputs = ["file('subtype.csv')", "file('resistance.csv')", "file('id_map.tsv')", "'TEST'", fasta_input]
    if process == 'REPORTHUMANFASTA':
        inputs = [f"file('{name}.csv')" for name in csv_names[:8]] + [
            "file('id_map.tsv')", "'TEST'", "'test-version'", fasta_input,
            "'test-instrument'", "'unused-samplesheet'", "file('reassortment.csv')", "file('subclade.csv')",
        ]
    (tmp_path / 'main.nf').write_text(
        f"include {{ {process} }} from '{ROOT}/modules/local/{module}/main'\n"
        f"workflow {{\n    {process}({', '.join(inputs)})\n"
        f"    {process}.out.{port}.view {{ 'FINAL_FASTA_READY' }}\n}}\n"
    )
    (tmp_path / 'nextflow.config').write_text(
        f"params.outdir = '{tmp_path}/published'\nparams.publish_dir_mode = 'copy'\n"
        "params.multiqc_title = null\n"
        f"includeConfig '{ROOT}/conf/modules.config'\nprocess.executor = 'local'\n"
        "docker.enabled = false\nconda.enabled = false\n"
    )
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    helpers = {
        'reportfasta.py': "#!/bin/sh\nprintf 'Sample,Subtype\\nUID1,H1N1\\n' > merged_report.csv\n",
        'report_QC_calculation.py': '#!/bin/sh\ncp "$1" "$3"\n',
    }
    for name, source in helpers.items():
        (bin_dir / name).write_text(source)
        (bin_dir / name).chmod(0o755)
    shutil.copy(ROOT / 'bin/drug_resistance_report.py', bin_dir)
    (bin_dir / 'drug_resistance_report.py').chmod(0o755)
    result = subprocess.run(
        ['nextflow', '-C', 'nextflow.config', 'run', 'main.nf', '-ansi-log', 'false'],
        cwd=tmp_path, env={**os.environ, 'NXF_OFFLINE': 'true'},
        capture_output=True, text=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'FINAL_FASTA_READY' in result.stdout
    assert (tmp_path / 'published' / destination / 'TEST.fasta').read_text() == expected

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

BIOPYTHON_AVAILABLE = importlib.util.find_spec("Bio") is not None


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MUTATION_FINDER = REPOSITORY_ROOT / "bin" / "mutation_finder.py"
REPORT_SCRIPTS = [
    REPOSITORY_ROOT / "bin" / "report.py",
    REPOSITORY_ROOT / "bin" / "reportfasta.py",
    REPOSITORY_ROOT / "bin" / "reportavian.py",
]
QC_REPORT_SCRIPT = REPOSITORY_ROOT / "bin" / "report_QC_calculation.py"


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def run_mutation_finder(tmp_path, mutation_type, segment, reference_header, sequence="ACDE"):
    subtype = "H3N2"
    reference_dir = tmp_path / "references" / mutation_type / subtype
    reference_dir.mkdir(parents=True)
    (reference_dir / f"{segment}.fasta").write_text(
        f">{reference_header}\nACDE\n",
        encoding="utf-8",
    )

    sequence_file = tmp_path / "sample.fasta"
    if sequence is not None:
        sequence_file.write_text(f">sample|translated\n{sequence}\n", encoding="utf-8")
    else:
        sequence_file.write_text("", encoding="utf-8")

    output_file = tmp_path / "mutations.csv"
    subprocess.run(
        [
            sys.executable,
            str(MUTATION_FINDER),
            str(sequence_file),
            str(tmp_path / "references"),
            segment,
            subtype,
            str(output_file),
            mutation_type,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output_file


@pytest.mark.skipif(not BIOPYTHON_AVAILABLE, reason="Biopython is not installed")
@pytest.mark.parametrize(
    ("mutation_type", "segment", "reference_header", "column", "expected"),
    [
        ("human", "HA1", "A/Cambodia/e0826360/2020_HA1", "Mutation reference", "A/Cambodia/e0826360/2020"),
        ("mamailian", "NS", "A/Goose/Guangdong/1/96_NS", "Mutation reference", "A/Goose/Guangdong/1/96"),
        (
            "human_vaccine",
            "NA",
            "A/Croatia/1XXXV/2023_NA",
            "Vaccine mutation reference",
            "A/Croatia/1XXXV/2023",
        ),
    ],
)
def test_mutation_report_uses_actual_reference_header(
    tmp_path, mutation_type, segment, reference_header, column, expected
):
    output_file = run_mutation_finder(tmp_path, mutation_type, segment, reference_header)

    report_fields, report_rows = read_csv(output_file.with_name("mutations_report.csv"))
    assert column in report_fields
    assert report_rows[0][column] == expected

    full_report_fields, full_report_rows = read_csv(output_file.with_name("mutations_full_mutation_list_report.csv"))
    assert column in full_report_fields
    assert full_report_rows[0][column] == expected


@pytest.mark.skipif(not BIOPYTHON_AVAILABLE, reason="Biopython is not installed")
@pytest.mark.parametrize(
    ("mutation_type", "column"),
    [
        ("human", "Mutation reference"),
        ("human_vaccine", "Vaccine mutation reference"),
    ],
)
def test_empty_mutation_reports_include_reference_column(tmp_path, mutation_type, column):
    output_file = run_mutation_finder(
        tmp_path,
        mutation_type,
        "HA1",
        "A/Example/1/2020_HA1",
        sequence=None,
    )

    report_fields, report_rows = read_csv(output_file.with_name("mutations_report.csv"))
    assert column in report_fields
    assert report_rows == []

    full_report_fields, full_report_rows = read_csv(output_file.with_name("mutations_full_mutation_list_report.csv"))
    assert column in full_report_fields
    assert full_report_rows == []


@pytest.mark.parametrize("report_script", REPORT_SCRIPTS, ids=lambda path: path.name)
def test_report_mergers_preserve_reference_columns(tmp_path, report_script):
    (tmp_path / "human_report.csv").write_text(
        "Sample,Mutation reference,Characterisation_Result\n"
        "sample,A/Cambodia/e0826360/2020,A/Victoria/4897/2022-like + R45K\n",
        encoding="utf-8",
    )
    (tmp_path / "vaccine_report.csv").write_text(
        "Sample,Vaccine mutation reference\nsample,A/Croatia/1XXXV/2023\n",
        encoding="utf-8",
    )

    command = [sys.executable, str(report_script)]
    if report_script.name == "report.py":
        samplesheet = tmp_path / "samplesheet.tsv"
        samplesheet.write_text("SequenceID\nsample\n", encoding="utf-8")
        command.append(str(samplesheet))

    subprocess.run(command, cwd=tmp_path, check=True, capture_output=True, text=True)

    fields, rows = read_csv(tmp_path / "merged_report.csv")
    assert "Mutation reference" in fields
    assert "Vaccine mutation reference" in fields
    assert "Characterisation_Result" in fields
    assert rows[0]["Mutation reference"] == "A/Cambodia/e0826360/2020"
    assert rows[0]["Vaccine mutation reference"] == "A/Croatia/1XXXV/2023"
    assert rows[0]["Characterisation_Result"] == "A/Victoria/4897/2022-like + R45K"


@pytest.mark.parametrize("report_script", REPORT_SCRIPTS, ids=lambda path: path.name)
def test_empty_report_mergers_include_reference_columns(tmp_path, report_script):
    command = [sys.executable, str(report_script)]
    if report_script.name == "report.py":
        samplesheet = tmp_path / "samplesheet.tsv"
        samplesheet.write_text("SequenceID\nsample\n", encoding="utf-8")
        command.append(str(samplesheet))

    subprocess.run(command, cwd=tmp_path, check=True, capture_output=True, text=True)

    fields, _ = read_csv(tmp_path / "merged_report.csv")
    assert "Mutation reference" in fields
    assert "Vaccine mutation reference" in fields


def test_final_qc_preserves_characterisation_columns(tmp_path):
    input_path = tmp_path / "merged_report.csv"
    output_path = tmp_path / "final_report.csv"
    input_path.write_text(
        "Sample,Subtype,Coverage-HA,Coverage-NA,Characterisation_Result\n"
        "sample,H1N1,100,100,A/Victoria/4897/2022-like + R45K\n",
        encoding="utf-8",
    )

    subprocess.run(
        [sys.executable, str(QC_REPORT_SCRIPT), str(input_path), "-o", str(output_path)],
        check=True,
        capture_output=True,
        text=True,
    )

    fields, rows = read_csv(output_path)
    assert "Characterisation_Result" in fields
    assert rows[0]["Characterisation_Result"] == "A/Victoria/4897/2022-like + R45K"

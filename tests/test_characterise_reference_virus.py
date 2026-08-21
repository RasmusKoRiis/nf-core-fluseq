import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPOSITORY_ROOT / "bin" / "characterise_reference_virus.py"
GUIDELINES = REPOSITORY_ROOT / "assets" / "characterisation_guidelines"

SPEC = importlib.util.spec_from_file_location("characterise_reference_virus", SCRIPT_PATH)
CHARACTERISE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHARACTERISE)


def nomenclature_row(profile, clade, subclade, unique="NA", missing="NA", fraction="1.000"):
    return {
        "Sample": "sample",
        "Subclade_Nomenclature_Profile": profile,
        "Subclade_Nomenclature_Clade": clade,
        "Subclade_Nomenclature_Subclade": subclade,
        "Subclade_Nomenclature_Closest_Subclade_Missing_Mutations": missing,
        "Subclade_Nomenclature_Unique_Mutations": unique,
        "Subclade_Nomenclature_Subclade_Match_Fraction": fraction,
    }


@pytest.mark.parametrize(
    ("subtype", "row", "reference", "status", "guideline_extras"),
    [
        (
            "H1N1",
            nomenclature_row("H1N1pdm_HA", "5a.2a.1", "D.3.1"),
            "A/Missouri/11/2025",
            "Exact reporting category",
            "NA",
        ),
        (
            "H1N1",
            nomenclature_row("H1N1pdm_HA", "5a.2a.1", "D.1"),
            "A/Victoria/4897/2022",
            "Derived reporting category",
            "R45K",
        ),
        (
            "H1N1",
            nomenclature_row("H1N1pdm_HA", "5a.2a", "C.1.9.1"),
            "A/Lisboa/188/2023",
            "Derived reporting category",
            "P137S",
        ),
        (
            "H3N2",
            nomenclature_row("H3N2_HA", "2a.3a.1", "J.1.1"),
            "A/Thailand/8/2022",
            "Derived reporting category",
            "I25V; V347M; S145N",
        ),
        (
            "H3N2",
            nomenclature_row("H3N2_HA", "2a.3a.1", "K"),
            "A/Norway/8765/2025",
            "Exact reporting category",
            "NA",
        ),
        (
            "VIC",
            nomenclature_row("B-Vic_HA", "V1A.3a.2", "C.5.4"),
            "B/Stockholm/3/2022",
            "Derived reporting category",
            "V117I; E128K; A154T; K326R",
        ),
    ],
)
def test_reporting_reference_classification(subtype, row, reference, status, guideline_extras):
    result = CHARACTERISE.characterise_row(row, subtype, str(GUIDELINES))

    assert result["Characterisation_Reference_Virus"] == reference
    assert result["Characterisation_Status"] == status
    assert result["Characterisation_Guideline_Extra_Mutations"] == guideline_extras
    assert result["Characterisation_Result"].startswith(f"{reference}-like")


def test_nonreporting_reference_is_kept_as_closest_context():
    row = nomenclature_row("H1N1pdm_HA", "5a.2a.1", "D.1")
    result = CHARACTERISE.characterise_row(row, "H1N1", str(GUIDELINES))

    assert result["Characterisation_Reference_Virus"] == "A/Victoria/4897/2022"
    assert result["Characterisation_Closest_Guideline_Reference"] == "A/Netherlands/10481/2024"


def test_sample_specific_amino_acid_mutations_are_appended():
    row = nomenclature_row(
        "H1N1pdm_HA",
        "5a.2a.1",
        "D.1",
        unique="nuc:10A;HA1:R45K;HA1:A99T;HA2:S12N",
    )
    result = CHARACTERISE.characterise_row(row, "H1N1", str(GUIDELINES))

    assert result["Characterisation_Sample_Extra_Mutations"] == "HA1:R45K; HA1:A99T; HA2:S12N"
    assert result["Characterisation_All_Extra_Mutations"] == "R45K; HA1:A99T; HA2:S12N"
    assert "nuc:" not in result["Characterisation_Result"]


def test_incomplete_subclade_call_is_provisional():
    row = nomenclature_row(
        "H1N1pdm_HA",
        "5a.2a.1",
        "D.1",
        missing="HA1:45K",
        fraction="0.500",
    )
    result = CHARACTERISE.characterise_row(row, "H1N1", str(GUIDELINES))

    assert result["Characterisation_Status"] == "Review - incomplete subclade call"
    assert result["Characterisation_Result"].startswith("Provisional: ")


def test_yamagata_table_does_not_invent_a_reporting_category():
    row = nomenclature_row("Unsupported subtype", "Y3", "NA")
    result = CHARACTERISE.characterise_row(row, "YAM", str(GUIDELINES))

    assert result["Characterisation_Profile"] == "B/Yamagata"
    assert result["Characterisation_Reference_Virus"] == "NA"
    assert result["Characterisation_Closest_Guideline_Reference"] == "B/Phuket/3073/2013"
    assert result["Characterisation_Status"] == "No reporting categories in guideline"
    assert result["Characterisation_Result"] == (
        "No reporting category; closest guideline reference B/Phuket/3073/2013-like"
    )


def test_cli_preserves_subclade_columns_and_appends_characterisation(tmp_path):
    input_path = tmp_path / "subclade.csv"
    subtype_path = tmp_path / "subtype.txt"
    output_path = tmp_path / "characterisation.csv"
    row = nomenclature_row("B-Vic_HA", "V1A.3a.2", "C.5.6.1")

    with input_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    subtype_path.write_text("VIC\n", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--input",
            str(input_path),
            "--subtype-file",
            str(subtype_path),
            "--guidelines-dir",
            str(GUIDELINES),
            "--output",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with output_path.open(encoding="utf-8", newline="") as handle:
        output = next(csv.DictReader(handle))
    assert output["Subclade_Nomenclature_Subclade"] == "C.5.6.1"
    assert output["Characterisation_Reference_Virus"] == "B/ENG/120/2025"
    assert output["Characterisation_Status"] == "Exact reporting category"

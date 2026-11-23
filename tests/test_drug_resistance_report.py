import csv
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
REPORT_SCRIPT = REPOSITORY_ROOT / "bin" / "drug_resistance_report.py"


def write(path, content):
    path.write_text(content, encoding="utf-8")
    return path


def test_drug_resistance_report_contains_only_focused_results(tmp_path):
    id_map = write(
        tmp_path / "id_map.tsv",
        "SampleID\tOriginalName\nUID1\tA/Test/1/2026\nUID2\tB/Test/2/2026\n",
    )
    subtype = write(
        tmp_path / "subtype.csv",
        "Sample,Subtype\nUID1,H1N1\nUID2,VIC\n",
    )
    m2 = write(
        tmp_path / "m2.csv",
        "Sample,M2 inhibtion mutations\nUID1,No matching mutations found\n",
    )
    na = write(
        tmp_path / "na.csv",
        "Sample,NA1 inhibtion mutations\nUID1,H275Y\n",
    )
    pa = write(
        tmp_path / "pa.csv",
        "Sample,PA inhibtion mutations\nUID1,No matching mutations found\n",
    )
    output = tmp_path / "run42_drug_resistance_report.csv"

    subprocess.run(
        [
            sys.executable,
            str(REPORT_SCRIPT),
            "--id-map",
            str(id_map),
            "--subtype",
            str(subtype),
            "--resistance",
            str(m2),
            str(na),
            str(pa),
            "--output",
            str(output),
        ],
        check=True,
    )

    with output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert reader.fieldnames == [
        "SampleID",
        "OriginalName",
        "Subtype",
        "M2 inhibtion mutations",
        "NA inhibtion mutations",
        "PA inhibtion mutations",
        "DR_Res_Adamantine",
        "DR_Res_Oseltamivir",
        "DR_Res_Zanamivir",
        "DR_Res_Peramivir",
        "DR_Res_Laninamivir",
        "DR_Res_Baloxavir",
        "DR_M2_Mut",
        "DR_NA_Mut",
        "DR_PA_Mut",
    ]

    assert rows[0]["OriginalName"] == "A/Test/1/2026"
    assert rows[0]["Subtype"] == "H1N1"
    assert rows[0]["NA inhibtion mutations"] == "H275Y"
    assert rows[0]["DR_Res_Adamantine"] == "AANI"
    assert rows[0]["DR_Res_Oseltamivir"] == "AARI"
    assert rows[0]["DR_Res_Baloxavir"] == "AANS"
    assert rows[0]["DR_NA_Mut"] == "H275Y"
    assert rows[1]["OriginalName"] == "B/Test/2/2026"
    assert rows[1]["Subtype"] == "VIC"
    assert rows[1]["DR_Res_Oseltamivir"] == "NA"


@pytest.mark.parametrize("script", ["report.py", "reportfasta.py", "drug_resistance_report.py"])
def test_reporting_codes_preserve_mutation_details(tmp_path, script):
    source = write(
        tmp_path / "results.csv",
        "Sample,M2 inhibtion mutations,NA inhibtion mutations,PA inhibtion mutations\n"
        "matched,M2 test match,NA test match,PA test match\n"
        "unmatched,No matching mutations found,No matching mutations found,No matching mutations found\n"
        "missing,NA,NA,NA\n",
    )
    command = [sys.executable, str(REPOSITORY_ROOT / "bin" / script)]
    output = tmp_path / "merged_report.csv"
    sample_key = "Sample"
    if script == "report.py":
        samplesheet = write(tmp_path / "samples.tsv", "SequenceID\nmatched\nunmatched\nmissing\n")
        command.append(str(samplesheet))
    elif script == "drug_resistance_report.py":
        id_map = write(
            tmp_path / "id_map.tsv",
            "SampleID\tOriginalName\nmatched\tmatched\nunmatched\tunmatched\nmissing\tmissing\n",
        )
        command.extend(["--id-map", str(id_map), "--resistance", str(source), "--output", str(output)])
        sample_key = "SampleID"

    subprocess.run(command, cwd=tmp_path, check=True, capture_output=True, text=True)
    with output.open(newline="", encoding="utf-8") as handle:
        rows = {row[sample_key]: row for row in csv.DictReader(handle)}

    for drug, match_code, no_match_code in [
        ("Adamantine", "AAHRI", "AANI"),
        ("Oseltamivir", "AARI", "AANI"),
        ("Zanamivir", "AARI", "AANI"),
        ("Peramivir", "AARI", "AANI"),
        ("Laninamivir", "AARI", "AANI"),
        ("Baloxavir", "AARS", "AANS"),
    ]:
        column = f"DR_Res_{drug}"
        assert rows["matched"][column] == match_code
        assert rows["unmatched"][column] == no_match_code
        assert rows["missing"][column] == "NA"

    for source_name in ["M2", "NA", "PA"]:
        column = f"DR_{source_name}_Mut"
        assert rows["matched"][column] == f"{source_name} test match"
        assert rows["unmatched"][column] == "No Mutations"
        assert rows["missing"][column] == "NA"

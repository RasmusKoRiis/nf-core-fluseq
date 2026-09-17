import csv
import subprocess
import sys
from pathlib import Path


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
    assert rows[0]["DR_Res_Oseltamivir"] == "Review"
    assert rows[0]["DR_Res_Baloxavir"] == "AANS"
    assert rows[0]["DR_NA_Mut"] == "H275Y"
    assert rows[1]["OriginalName"] == "B/Test/2/2026"
    assert rows[1]["Subtype"] == "VIC"
    assert rows[1]["DR_Res_Oseltamivir"] == "NA"

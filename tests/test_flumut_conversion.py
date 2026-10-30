"""Check formatting of grouped source entries in the FluMut CSV export."""

import csv
from pathlib import Path
import subprocess
import sys


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "flumut_conversion.py"


def test_grouped_entries_keep_effect_alignment_and_csv_columns(tmp_path):
    source = tmp_path / "markers.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Marker", "Mutations in your sample", "Effect"])
        writer.writerow(["sample", "HA1", "HA1:first", "first effect"])
        writer.writerow(["sample", "HA1", "HA1:second,HA1:third", "second effect, detail"])
        writer.writerow(["sample", "HA1", "HA1:second,HA1:third", "duplicate"])
        writer.writerow(["sample", "PB1", "PB1:other,PB1:value", "other effect"])

    subprocess.run([sys.executable, str(SCRIPT), str(source), "sample"], cwd=tmp_path, check=True)
    output = (tmp_path / "sample_flumut_report.csv").read_text(encoding="utf-8")
    header, record = output.splitlines()
    assert len(header.split(",")) == len(record.split(","))
    rows = list(csv.DictReader(output.splitlines()))
    assert len(rows) == 1
    assert rows[0]["Flumut_HA1"] == "first;second+HA1:third"
    assert rows[0]["Effect_HA1"] == "first effect;second effect / detail"
    assert rows[0]["Flumut_PB1"] == "other+PB1:value"
    assert rows[0]["Effect_PB1"] == "other effect"

import csv
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY_ROOT / "bin" / "detect_reassortment.py"
METADATA = REPOSITORY_ROOT / "bin" / "reassortment_reference_metadata.csv"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def run_screen(tmp_path, hits, use_metadata=True):
    blast_path = tmp_path / "hits.tsv"
    output_path = tmp_path / "summary.csv"
    lines = []
    for segment, (subject, identity) in hits.items():
        fields = [
            f"sample_{segment}",
            subject,
            identity,
            1000,
            0,
            0,
            1,
            1000,
            1,
            1000,
            "1e-50",
            500,
        ]
        lines.append("\t".join(str(value) for value in fields))
    blast_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")

    command = [
        sys.executable,
        str(SCRIPT),
        "--blast",
        str(blast_path),
        "--output",
        str(output_path),
        "--sample",
        "sample",
    ]
    if use_metadata:
        command.extend(["--metadata", str(METADATA)])
    subprocess.run(command, check=True, capture_output=True, text=True)

    with output_path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle))


def enriched_hits(origin="HUMAN", subtype="H1N1", strain="A/Test/1/2026"):
    return {
        segment: (f"{origin}|{subtype}|{strain}|EPI_TEST", 99.0)
        for segment in SEGMENTS
    }


@pytest.mark.parametrize(
    ("subject", "expected"),
    [
        (
            "A/Victoria/2570/2019|pdm09|{segment}|EPI_ISL_528951",
            "HUMAN|H1N1|A/Victoria/2570/2019(99.0%)",
        ),
        (
            "B/Brisbane/60/2008|{segment}|EPI_ISL_246494",
            "HUMAN|B/Victoria|B/Brisbane/60/2008(99.0%)",
        ),
    ],
)
def test_legacy_headers_are_enriched_from_accession_metadata(tmp_path, subject, expected):
    hits = {
        segment: (subject.format(segment=segment), 99.0)
        for segment in SEGMENTS
    }
    row = run_screen(tmp_path, hits)

    assert row["PB2"] == expected
    assert row["Origins"] == "HUMAN"
    assert row["Reassortment"] == "No"


@pytest.mark.parametrize(
    ("change", "status", "conclusion"),
    [
        (
            None,
            "No",
            "CONSISTENT - all segments match one HUMAN H1N1 reference strain",
        ),
        (
            "strain",
            "Yes",
            "REVIEW - possible within-subtype reassortment; multiple HUMAN H1N1 reference strains",
        ),
        (
            "subtype",
            "Yes",
            "ALERT - human subtype discordance (H1N1,H3N2)",
        ),
        (
            "origin",
            "Yes",
            "ALERT - mixed origins (AVIAN,HUMAN) and subtypes (H1N1,H5N1)",
        ),
    ],
)
def test_actionable_conclusion_categories(tmp_path, change, status, conclusion):
    hits = enriched_hits()
    if change == "strain":
        hits["NS"] = ("HUMAN|H1N1|A/Test/2/2026|EPI_OTHER", 99.0)
    elif change == "subtype":
        hits["NS"] = ("HUMAN|H3N2|A/Test/3/2026|EPI_OTHER", 99.0)
    elif change == "origin":
        hits["NS"] = ("AVIAN|H5N1|A/duck/Test/1/2026|EPI_OTHER", 99.0)

    row = run_screen(tmp_path, hits)

    assert row["Reassortment"] == status
    assert row["Conclusion"] == conclusion


def test_all_non_human_segments_are_flagged(tmp_path):
    row = run_screen(
        tmp_path,
        enriched_hits(origin="AVIAN", subtype="H5N1", strain="A/duck/Test/1/2026"),
    )

    assert row["Reassortment"] == "No"
    assert row["Conclusion"] == "FLAG - all segments match one non-human AVIAN H5N1 reference strain"


def test_missing_and_low_identity_segments_are_inconclusive(tmp_path):
    hits = enriched_hits()
    del hits["NS"]
    hits["NA"] = (hits["NA"][0], 79.0)

    row = run_screen(tmp_path, hits)

    assert row["NS"] == "Missing"
    assert row["NA"].startswith("TooLow(79.0%):HUMAN|H1N1|")
    assert row["Reassortment"] == "Unknown"
    assert row["Conclusion"] == "INCONCLUSIVE - missing segments: NS; low-identity segments: NA"


def test_empty_blast_result_reports_all_segments_missing(tmp_path):
    row = run_screen(tmp_path, {})

    assert all(row[segment] == "Missing" for segment in SEGMENTS)
    assert row["Reassortment"] == "Unknown"
    assert row["Conclusion"].startswith("INCONCLUSIVE - missing segments:")


def test_reference_metadata_has_unique_accessions():
    with METADATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    accessions = [row["accession"] for row in rows]
    assert len(rows) == 37
    assert len(accessions) == len(set(accessions))
    assert {row["origin"] for row in rows} == {"HUMAN"}
    assert {row["subtype"] for row in rows} == {"H1N1", "H3N2", "B/Victoria"}

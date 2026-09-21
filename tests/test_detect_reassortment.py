import csv
import runpy
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY_ROOT / "bin" / "detect_reassortment.py"
METADATA = REPOSITORY_ROOT / "bin" / "reassortment_reference_metadata.csv"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def run_screen(tmp_path, hits, use_metadata=True, query_sequences=None):
    blast_path = tmp_path / "hits.tsv"
    fasta_path = tmp_path / "query.fasta"
    output_path = tmp_path / "summary.csv"
    if query_sequences is None:
        query_sequences = {segment: "ACGT" * 250 for segment in hits}
    fasta_path.write_text(
        "".join(f">sample_{segment} synthetic record\n{sequence}\n" for segment, sequence in query_sequences.items()),
        encoding="utf-8",
    )
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
        "--fasta",
        str(fasta_path),
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


def enriched_hits(origin="HUMAN-SEASONAL", subtype="H1N1", strain="A/Test/1/2026"):
    return {segment: (f"{origin}|{subtype}|{strain}|EPI_TEST", 99.0) for segment in SEGMENTS}


@pytest.mark.parametrize(
    ("subject", "expected"),
    [
        (
            "A/Victoria/2570/2019|pdm09|{segment}|EPI_ISL_528951",
            "HUMAN-SEASONAL|H1N1|A/Victoria/2570/2019(M:99.0%/N:0.0%)",
        ),
        (
            "B/Brisbane/60/2008|{segment}|EPI_ISL_246494",
            "HUMAN-SEASONAL|B/Victoria|B/Brisbane/60/2008(M:99.0%/N:0.0%)",
        ),
    ],
)
def test_legacy_headers_are_enriched_from_accession_metadata(tmp_path, subject, expected):
    hits = {segment: (subject.format(segment=segment), 99.0) for segment in SEGMENTS}
    row = run_screen(tmp_path, hits)

    assert row["PB2"] == expected
    assert row["Origins"] == "HUMAN-SEASONAL"
    assert row["Reassortment"] == "No"


@pytest.mark.parametrize(
    ("change", "status", "conclusion"),
    [
        (
            None,
            "No",
            "CONSISTENT - all eight segments match HUMAN-SEASONAL H1N1 references",
        ),
        (
            "strain",
            "Yes",
            "CONSISTENT - all eight segments match HUMAN-SEASONAL H1N1 references",
        ),
        (
            "subtype",
            "Yes",
            "ALERT - subtype discordance: H1N1,H3N2",
        ),
        (
            "origin",
            "Yes",
            "ALERT - origin(s) not confirmed seasonal human: AVIAN; "
            "mixed origins: AVIAN,HUMAN-SEASONAL; subtype discordance: H1N1,H5N1",
        ),
    ],
)
def test_actionable_conclusion_categories(tmp_path, change, status, conclusion):
    hits = enriched_hits()
    if change == "strain":
        hits["NS"] = ("HUMAN-SEASONAL|H1N1|A/Test/2/2026|EPI_OTHER", 99.0)
    elif change == "subtype":
        hits["NS"] = ("HUMAN-SEASONAL|H3N2|A/Test/3/2026|EPI_OTHER", 99.0)
    elif change == "origin":
        hits["NS"] = ("AVIAN|H5N1|A/duck/Test/1/2026|EPI_OTHER", 99.0)

    row = run_screen(tmp_path, hits)

    assert row["Reassortment"] == status
    assert row["Conclusion"] == conclusion


def test_all_non_human_segments_alert(tmp_path):
    row = run_screen(
        tmp_path,
        enriched_hits(origin="AVIAN", subtype="H5N1", strain="A/duck/Test/1/2026"),
    )

    assert row["Reassortment"] == "No"
    assert row["Conclusion"] == "ALERT - origin(s) not confirmed seasonal human: AVIAN"


def test_missing_and_low_identity_segments_alert(tmp_path):
    hits = enriched_hits()
    del hits["NS"]
    hits["NA"] = (hits["NA"][0], 79.0)

    row = run_screen(tmp_path, hits)

    assert row["NS"] == "Missing"
    assert row["NA"].startswith("TooLow(M:79.0%/N:0.0%):HUMAN-SEASONAL|H1N1|")
    assert row["Reassortment"] == "Unknown"
    assert row["Conclusion"] == "ALERT - missing segments: NS; low-identity segments: NA"


def test_empty_blast_result_reports_all_segments_missing(tmp_path):
    row = run_screen(tmp_path, {})

    assert all(row[segment] == "Missing" for segment in SEGMENTS)
    assert row["Reassortment"] == "Unknown"
    assert row["Conclusion"].startswith("ALERT - missing segments:")


@pytest.mark.parametrize("subtype", ["H1N1", "H3N2", "B/Victoria", "B/Yamagata"])
def test_complete_seasonal_profile_is_consistent(tmp_path, subtype):
    row = run_screen(tmp_path, enriched_hits(subtype=subtype))

    assert row["Conclusion"] == f"CONSISTENT - all eight segments match HUMAN-SEASONAL {subtype} references"


@pytest.mark.parametrize("origin", ["HUMAN", "HUMAN-NONSEASONAL", "SWINE", "AVIAN"])
def test_generic_human_or_other_origin_is_not_assumed_seasonal(tmp_path, origin):
    row = run_screen(tmp_path, enriched_hits(origin=origin))

    assert row["Conclusion"] == f"ALERT - origin(s) not confirmed seasonal human: {origin}"


def test_known_accession_refines_generic_human_origin(tmp_path):
    hits = {
        segment: ("HUMAN|H3N2|A/Singapore/GP20238/2024|EPI_ISL_20307689", 99.9)
        for segment in SEGMENTS
    }
    sequences = {segment: "ACGT" * 250 + "N" * 125 + "n" * 125 for segment in SEGMENTS}
    row = run_screen(tmp_path, hits, query_sequences=sequences)

    assert row["PB2"] == "HUMAN-SEASONAL|H3N2|A/Singapore/GP20238/2024(M:99.9%/N:20.0%)"
    assert row["Conclusion"] == "CONSISTENT - all eight segments match HUMAN-SEASONAL H3N2 references"


@pytest.mark.parametrize(("origin", "subtype"), [("AVIAN", "H3N2"), ("HUMAN", "H1N1")])
def test_accession_metadata_does_not_override_conflicting_explicit_annotation(tmp_path, origin, subtype):
    hits = {
        segment: (f"{origin}|{subtype}|A/Test/1/2026|EPI_ISL_20307689", 99.0)
        for segment in SEGMENTS
    }
    row = run_screen(tmp_path, hits)

    assert row["Origins"] == origin
    assert row["Conclusion"].startswith("ALERT - ")


@pytest.mark.parametrize("field", ["origin", "subtype", "strain"])
def test_incomplete_reference_metadata_alerts(tmp_path, field):
    annotation = {"origin": "HUMAN-SEASONAL", "subtype": "H1N1", "strain": "A/Test/1/2026"}
    annotation[field] = "UNKNOWN"
    row = run_screen(tmp_path, enriched_hits(**annotation), use_metadata=False)

    assert row["Reassortment"] == "Unknown"
    assert row["Conclusion"].startswith("ALERT - reference metadata missing for:")


@pytest.mark.parametrize(("identity", "accepted"), [(79.94, False), (79.99, False), (80.0, True), (80.01, True)])
def test_identity_threshold_uses_unrounded_blast_value(tmp_path, identity, accepted):
    hits = enriched_hits()
    hits["NS"] = (hits["NS"][0], identity)
    row = run_screen(tmp_path, hits)

    assert row["NS"].startswith("TooLow") is not accepted
    assert row["Reassortment"] == ("No" if accepted else "Unknown")
    assert row["Conclusion"].startswith("CONSISTENT" if accepted else "ALERT")


def test_n_content_is_per_full_query_record_not_alignment_length(tmp_path):
    sequences = {segment: "ACGT" * 250 for segment in SEGMENTS}
    # The mocked BLAST alignment covers positions 1..1000; Ns are outside it.
    sequences["PB2"] += "\nnN" * 125
    sequences["NA"] = "acgtNN\nacgt nn\n"  # 4 Ns / 12 sequence characters.
    row = run_screen(tmp_path, enriched_hits(), query_sequences=sequences)

    assert row["PB2"].endswith("(M:99.0%/N:20.0%)")
    assert row["NA"].endswith("(M:99.0%/N:33.3%)")
    assert row["NS"].endswith("(M:99.0%/N:0.0%)")


def test_n_content_uses_the_selected_blast_query(tmp_path):
    hits = enriched_hits()
    hits["PB2_alternative"] = (hits["PB2"][0], 99.9)
    sequences = {segment: "ACGT" * 250 for segment in hits}
    sequences["PB2_alternative"] += "N" * 250
    row = run_screen(tmp_path, hits, query_sequences=sequences)

    assert row["PB2"].endswith("(M:99.9%/N:20.0%)")


def test_low_identity_hit_still_includes_n_content(tmp_path):
    hits = enriched_hits()
    hits["MP"] = (hits["MP"][0], 75.0)
    sequences = {segment: "ACGT" * 250 for segment in hits}
    sequences["MP"] += "N" * 250
    row = run_screen(tmp_path, hits, query_sequences=sequences)

    assert row["MP"].startswith("TooLow(M:75.0%/N:20.0%):")


@pytest.mark.parametrize("sequence", [None, ""])
def test_missing_or_empty_query_sequence_is_not_reported_as_zero_ns(tmp_path, sequence):
    sequences = {segment: "ACGT" * 250 for segment in SEGMENTS}
    if sequence is None:
        del sequences["NS"]
    else:
        sequences["NS"] = sequence

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_screen(tmp_path, enriched_hits(), query_sequences=sequences)
    assert "BLAST query 'sample_NS' is absent or empty" in error.value.stderr


def test_duplicate_query_identifiers_are_rejected(tmp_path):
    sequences = {segment: "ACGT" * 250 for segment in SEGMENTS}
    sequences["NS"] += "\n>sample_NS\nNNNN"

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_screen(tmp_path, enriched_hits(), query_sequences=sequences)
    assert "Duplicate query FASTA identifier: sample_NS" in error.value.stderr


def test_report_merge_and_surveillance_preserve_match_and_n_percentages(tmp_path):
    sequences = {segment: "ACGT" * 250 for segment in SEGMENTS}
    sequences["PB2"] += "N" * 250
    row = run_screen(tmp_path, enriched_hits(), query_sequences=sequences)
    subprocess.run(
        [sys.executable, str(REPOSITORY_ROOT / "bin/reportfasta.py")],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    with (tmp_path / "merged_report.csv").open(newline="", encoding="utf-8") as handle:
        merged = next(csv.DictReader(handle))
    for column in [*SEGMENTS, "Conclusion", "Reassortment", "Origins"]:
        assert merged[column] == row[column]

    surveillance = runpy.run_path(str(REPOSITORY_ROOT / "bin/surveillance_summary.py"))
    parsed = surveillance["parse_reassortment"]([str(tmp_path / "summary.csv")])["sample"]
    assert parsed["segments"]["PB2"] == row["PB2"]
    assert parsed["lineages"] == ["HUMAN-SEASONAL|H1N1|A/Test/1/2026"]


def test_reference_metadata_has_unique_accessions():
    with METADATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    accessions = [row["accession"] for row in rows]
    assert len(rows) == 37
    assert len(accessions) == len(set(accessions))
    assert {row["origin"] for row in rows} == {"HUMAN-SEASONAL"}
    assert {row["subtype"] for row in rows} == {"H1N1", "H3N2", "B/Victoria"}

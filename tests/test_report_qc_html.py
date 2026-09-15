"""QC rendering and real Nextflow module execution using invented report records."""

import csv
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin/report_qc_html.py"
TEMPLATE = ROOT / "bin/templates/report_qc.html"
spec = importlib.util.spec_from_file_location("report_qc_html", SCRIPT)
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


def sample(identifier="SYNTHETIC-01", **overrides):
    row = {
        "Sample": identifier,
        "RunID": "SYNTHETIC_RUN",
        "Instrument ID": "TEST",
        "Date": "2026-01-01",
        "Release Version": "test",
        "NGS_QC_Sum": "",
        "GISAID_Comment": "",
        "Subtype": "H3N2",
        "Subclade_Nomenclature_Subclade": "TEST_LABEL",
        "Subclade_Nomenclature_Subclade_Match_Fraction": "1",
        "Characterisation_Status": "Exact reporting category",
    }
    row.update({f"Coverage-{seg}": "100" for seg in qc.SEGMENTS})
    row.update({f"Nextclade QC {protein}": "good" for protein in qc.PROTEINS})
    row.update(overrides)
    return row


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def assess(**overrides):
    row = sample(**overrides)
    return qc.summarise_sample(row, list(row))


def test_existing_flags_and_additional_review_are_distinct(tmp_path):
    rows = [sample(), sample("SYNTHETIC-02", NGS_QC_Sum="PB2:FS"), sample("SYNTHETIC-03", **{"Nextclade QC NS": "NA"})]
    report = qc.build_report(write_csv(tmp_path / "run.csv", rows))
    assert report["counts"]["reported_qc"] == {"No flags": 2, "Review": 1}
    assert report["counts"]["assessment"] == {"Review": 2, "No flags": 1}
    assert report["samples"][1]["reported_qc"] == "PB2:FS"
    assert report["assessment"] == "Review required"
    assert report["samples"][2]["missing_nextclade_qc"] == ["NS"]


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("80", "No flags"),
        ("79.99", "Review"),
        ("0", "Review"),
        ("NA", "Review"),
        ("", "Review"),
        ("nan", "Review"),
        ("inf", "Review"),
        ("-1", "Review"),
        ("101", "Review"),
        ("bad", "Review"),
    ],
)
def test_coverage_boundary_and_invalid_values(raw, expected):
    result = assess(**{"Coverage-HA": raw})
    assert result["qc_status"] == expected
    assert result["reported_qc_status"] == "No flags"


def test_matrix_alias_precedence_and_fallback():
    assert assess(**{"Coverage-MP": "NA", "Coverage-M": "90"})["coverage"]["MP"] == 90
    assert assess(**{"Coverage-MP": "0", "Coverage-M": "100"})["coverage"]["MP"] == 0
    assert assess(**{"Coverage-MP": "bad", "Coverage-M": "100"})["coverage"]["MP"] is None


def test_display_preserves_coverage_boundary_and_noise_precision():
    assert qc.fmt(79.99, "%") == "79.99%"
    assert qc.fmt(0.00016, precision=5) == "0.00016"
    assert qc.fmt(1000) == "1,000"
    assert qc.fmt(0) == "0"
    assert qc.fmt(None) == "—"


def test_missing_qc_summary_is_not_a_pass():
    row = sample()
    del row["NGS_QC_Sum"]
    assert qc.summarise_sample(row, list(row))["reported_qc_status"] == "Unavailable"
    assert assess(NGS_QC_Sum="NA")["reported_qc_status"] == "Unavailable"
    row = {"Sample": "SYNTHETIC-EMPTY"}
    assert qc.summarise_sample(row, list(row))["qc_status"] == "Unavailable"


def test_source_review_and_nextclade_warnings():
    assert assess(GISAID_Comment="Review")["qc_status"] == "Review"
    assert assess(**{"Nextclade QC HA1": "mediocre"})["qc_status"] == "Review"
    assert assess(**{"Nextclade QC HA1": "unknown-status"})["qc_status"] == "Review"


def test_provisional_and_unknown_subclades_remain_labelled():
    assert assess(Subclade_Nomenclature_Subclade_Match_Fraction="0.95")["subclade_status"] == "Provisional"
    assert assess(Characterisation_Status="Review - incomplete subclade call")["subclade_status"] == "Provisional"
    assert assess(Subclade_Nomenclature_Subclade_Match_Fraction="NA")["subclade_status"] == "Match not assessed"
    assert assess(Subclade_Nomenclature_Subclade="NA")["subclade_status"] == "Unavailable"


def test_aliases_denominators_zero_and_missing_depth(tmp_path):
    rows = [
        sample("SYNTHETIC-01", Subtype="VIC", DEPTH_HA="0"),
        sample("SYNTHETIC-02", Subtype="VICVIC", DEPTH_HA="100"),
        sample("SYNTHETIC-03", Subtype="NA", DEPTH_HA="NA", **{"Coverage-HA": "NA"}),
    ]
    report = qc.build_report(write_csv(tmp_path / "run.csv", rows))
    assert report["counts"]["subtypes"] == {"B/Victoria": 2, "Unavailable": 1}
    assert report["segments"]["HA"] == {
        "at_least_80": 2,
        "below_80": 0,
        "missing_or_invalid": 1,
        "median_coverage": 100,
        "median_segment_reads": 50,
        "segment_reads_available": 2,
    }


def test_html_escapes_all_csv_text_and_preserves_quoted_cells(tmp_path):
    payload = '<script>alert("test")</script><img src=x onerror=alert(1)>@@TITLE@@'
    row = sample(
        payload,
        RunID=payload,
        Conclusion='REVIEW, quoted "text"\nand newline',
        Subtype=payload,
        **{"Characterisation_Reference_Virus": payload},
    )
    path = write_csv(tmp_path / "run.csv", [row])
    report = qc.build_report(path)
    rendered = qc.render(report, TEMPLATE.read_text())
    assert payload not in rendered
    assert qc.esc(payload) in rendered
    assert rendered.count("<script>") == 1
    assert report["samples"][0]["reassortment_conclusion"] == row["Conclusion"]
    assert report["source"]["sha256"] == qc.hashlib.sha256(path.read_bytes()).hexdigest()
    assert "<script src=" not in rendered and "<link " not in rendered


def test_qc_filter_exposes_every_assessment(tmp_path):
    class Options(HTMLParser):
        def __init__(self):
            super().__init__()
            self.in_qc = False
            self.in_option = False
            self.labels = []

        def handle_starttag(self, tag, attrs):
            if tag == "select":
                self.in_qc = dict(attrs).get("id") == "qc"
            if tag == "option" and self.in_qc:
                self.in_option = True

        def handle_data(self, text):
            if self.in_qc and self.in_option:
                self.labels.append(text)

        def handle_endtag(self, tag):
            if tag == "option":
                self.in_option = False
            if tag == "select":
                self.in_qc = False

    parser = Options()
    report = qc.build_report(write_csv(tmp_path / "run.csv", [sample()]))
    parser.feed(qc.render(report, TEMPLATE.read_text()))
    assert parser.labels == ["All samples", "Review", "No flags", "Unavailable"]


@pytest.mark.parametrize(
    "text", ["X\n1\n", "Sample,Sample\na,a\n", "Sample,X\na\n", "Sample\na,b\n", "Sample\na\na\n", "Sample,X\n,1\n"]
)
def test_malformed_or_ambiguous_input_rejected(tmp_path, text):
    path = tmp_path / "bad.csv"
    path.write_text(text)
    with pytest.raises(ValueError):
        qc.build_report(path)


def test_empty_report_and_mixed_metadata(tmp_path):
    path = write_csv(tmp_path / "empty.csv", [], ["Sample"])
    report = qc.build_report(path)
    assert report["assessment"] == "Not assessable"
    assert "0 samples" in qc.render(report, TEMPLATE.read_text())
    path = write_csv(tmp_path / "mixed.csv", [sample(), sample("SYNTHETIC-02", RunID="OTHER")])
    assert any("Multiple values in RunID" in message for message in qc.build_report(path)["warnings"])


def test_cli_outputs_are_reproducible_and_do_not_change_input(tmp_path):
    path = write_csv(tmp_path / "run with spaces.csv", [sample()])
    before = path.read_bytes()
    command = [sys.executable, str(SCRIPT), str(path), "--outdir", str(tmp_path / "output")]
    subprocess.run(command, check=True, capture_output=True, text=True)
    html_path = tmp_path / "output/run with spaces_qc.html"
    original_html = html_path.read_bytes()
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert original_html == html_path.read_bytes()
    assert before == path.read_bytes()
    report = json.loads((tmp_path / "output/run with spaces_qc.json").read_text())
    assert report["counts"]["samples"] == 1
    assert report["assessment"] == "No QC flags detected"


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="Nextflow is not installed")
def test_nextflow_module_execution_and_publish(tmp_path):
    path = write_csv(tmp_path / "synthetic.csv", [sample(), sample("SYNTHETIC-02", NGS_QC_Sum="PB1:LC")])
    config = tmp_path / "nextflow.config"
    config.write_text(
        f"params.outdir = '{tmp_path / 'published'}'\n"
        "params.publish_dir_mode = 'copy'\n"
        "params.multiqc_title = null\n"
        f"includeConfig '{ROOT / 'conf/modules.config'}'\n"
        "process.executor = 'local'\nprocess.errorStrategy = 'terminate'\n"
        "docker.enabled = false\nsingularity.enabled = false\n"
    )
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"include {{ REPORT_QC_HTML }} from '{ROOT / 'modules/local/report_qc_html/main'}'\n"
        "workflow {\n"
        f"    REPORT_QC_HTML(Channel.value(file('{path}')), Channel.value(file('{SCRIPT}')), Channel.value(file('{TEMPLATE}')))\n"
        "    REPORT_QC_HTML.out.html.view { report ->\n"
        "        assert report.name == 'synthetic_qc.html'\n"
        "        assert report.text.contains('SYNTHETIC-02')\n"
        "        'HTML_QC_COMPLETED'\n"
        "    }\n"
        "    REPORT_QC_HTML.out.summary.view { summary ->\n"
        "        assert summary.name == 'synthetic_qc.json'\n"
        "        'JSON_QC_COMPLETED'\n"
        "    }\n"
        "    REPORT_QC_HTML.out.versions.view { versions ->\n"
        "        assert versions.text.contains('report_qc_html: 1.0.0')\n"
        "        'VERSIONS_COMPLETED'\n"
        "    }\n"
        "}\n"
    )
    result = subprocess.run(
        [shutil.which("nextflow"), "-C", str(config), "run", str(workflow), "-ansi-log", "false"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=90,
        env={**os.environ, "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_SYNTAX_PARSER": "v1"},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for marker in ("HTML_QC_COMPLETED", "JSON_QC_COMPLETED", "VERSIONS_COMPLETED"):
        assert marker in result.stdout
    published = tmp_path / "published/reporthuman"
    assert sorted(p.name for p in published.iterdir()) == ["synthetic_qc.html", "synthetic_qc.json"]
    assert json.loads((published / "synthetic_qc.json").read_text())["counts"]["samples"] == 2

"""Master CSV import and transaction tests use synthetic administrative data only."""

import csv
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import merge_master_csv as merge
import update_master_csv as update


def table(headers, *rows):
    return headers.split(";"), [dict(zip(headers.split(";"), row)) for row in rows]


def save(path, data, delimiter=";", encoding="utf-8-sig"):
    with path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=data[0], delimiter=delimiter)
        writer.writeheader()
        writer.writerows(data[1])
    return path


def inputs(tmp_path):
    master = save(tmp_path / "master.csv", table(
        "Key;NGS_LW_ID;Comment_GC_Tessy;NGS_Reads_HA;NC_HA_Subclade;Custom",
        ["001", "SEQ-001", "keep; this", "100", "old annotation", "old"],
        ["002", "SEQ-002", "untouched", "88", "existing", "retain"],
    ))
    metadata = save(tmp_path / "metadata.csv", table(
        "Key;RekvNr;NGS_LW_ID;Pasient_Kjonn;Prove_Maanad;INF_Res;Added_metadata;NGS_Sekvens_Resultat",
        ["source-001", "001", "SEQ-001", "value", "01", "fresh", "æøå;\nquoted \"text\"", "stale"],
        ["source-003", "003", "SEQ-003", "value", "02", "fresh", "new", "stale"],
        ["source-004", "004", "SEQ-004", "value", "02", "not analysed", "new", "stale"],
    ))
    report = save(tmp_path / "report.csv", table(
        "Sample;RunID;Coverage-MP;subclade;Sekvens_Resultat;New_report_column",
        ["SEQ-001", "RUN02", "0", "", "current", "0010"],
        ["SEQ-003", "RUN02", "99.5", "annotation", "current", "NA"],
    ), delimiter=",")
    return master, metadata, report


def test_fastq_merge_update_insert_preserve_and_idempotence(tmp_path):
    master, metadata, report = inputs(tmp_path)
    original = master.read_bytes()
    stats = merge.merge_files(master, metadata, report, "fastq")
    headers, rows = merge.read_table(master)
    assert headers[:len(merge.schema()[0])] == merge.schema()[0]
    assert len(headers) == len(set(headers))
    assert [row["Key"] for row in rows] == ["001", "002", "003"]
    assert stats["added"] == 1 and stats["updated"] == 1 and stats["total"] == 3
    first, untouched, added = rows
    assert first["NGS_Run_ID"] == "RUN02"
    assert first["NGS_Coverage_M"] == "0"
    assert first["NGS_Reads_HA"] == "100"  # Column absent in report.
    assert first["NC_HA_Subclade"] == ""  # Column present but blank in report.
    assert first["Comment_GC_Tessy"] == "keep; this"
    assert first["Pasient_Kjnn"] == "value"
    assert first["Prove_M nad"] == "01"
    assert first["Metadata_Key"] == "source-001"
    assert first["Added_metadata"] == 'æøå;\nquoted "text"'
    assert first["NGS_Sekvens_Resultat"] == "current"
    assert first["New_report_column"] == "0010"
    assert untouched["Custom"] == "retain" and untouched["NGS_Run_ID"] == ""
    assert added["New_report_column"] == "NA"
    assert Path(stats["backup"]).read_bytes() == original
    assert master.read_bytes().startswith(b"\xef\xbb\xbf")
    assert not Path(str(master) + ".lock").exists()
    after = master.read_bytes()
    again = merge.merge_files(master, metadata, report, "fastq")
    assert again["unchanged"] == 2 and not again["written"]
    assert master.read_bytes() == after
    assert len(list(tmp_path.glob("master.csv.bak.*"))) == 1


def test_fasta_uses_original_name_keeps_depth_and_never_joins_uid(tmp_path):
    master, metadata, _ = inputs(tmp_path)
    report = save(tmp_path / "fasta.csv", table(
        "SampleID;OriginalName;Sample;RunID;Coverage-M;NC_HA_Subclade",
        ["UID123", "001", "UID123", "FASTA02", "98", "annotation"],
        ["UID456", "unmatched-name", "UID456", "FASTA02", "99", "annotation"],
    ), delimiter=",")
    stats = merge.merge_files(master, metadata, report, "fasta")
    _, rows = merge.read_table(master)
    assert [row["Key"] for row in rows] == ["001", "002", "unmatched-name"]
    assert rows[0]["NGS_LW_ID"] == "SEQ-001"
    assert rows[0]["SampleID"] == "UID123" and rows[0]["Sample"] == "UID123"
    assert rows[0]["NGS_Reads_HA"] == "100"
    assert rows[0]["NGS_Coverage_M"] == "98"
    assert all(rows[2][column] == "" for column in ("Metadata_Key", "INF_Res", "Pasient_Kjnn", "Added_metadata"))
    assert stats["without_metadata"] == 1


def test_blank_metadata_preserved_and_optional_nonempty_report(tmp_path):
    master = save(tmp_path / "master.csv", table("Key;INF_Res;NC_HA_Subclade", ["001", "retained", "old"]))
    metadata = save(tmp_path / "metadata.csv", table("RekvNr;INF_Res", ["001", ""]))
    report = save(tmp_path / "report.csv", table("OriginalName;subclade", ["001", ""]))
    merge.merge_files(master, metadata, report, "fasta", nonempty_report=True)
    row = merge.read_table(master)[1][0]
    assert row["INF_Res"] == "retained" and row["NC_HA_Subclade"] == "old"


@pytest.mark.parametrize("encoding", ["utf-8-sig", "cp1252", "utf-16"])
@pytest.mark.parametrize("delimiter", [";", ",", "\t"])
def test_csv_roundtrip(encoding, delimiter, tmp_path):
    expected = table("Key;Field", ["001", 'Blå; comma, tab\t\n"quoted"'])
    path = save(tmp_path / "input.csv", expected, delimiter, encoding)
    assert merge.read_table(path) == expected


@pytest.mark.parametrize("bad", [
    "duplicate_master", "duplicate_metadata", "duplicate_report", "ambiguous_sequence", "missing_sample",
    "missing_key", "missing_source", "empty_report", "column_collision", "duplicate_header", "ragged_row",
])
def test_invalid_inputs_never_modify_master(tmp_path, bad):
    master, metadata, report = inputs(tmp_path)
    path = report
    if bad == "duplicate_master":
        path = master
    elif bad == "duplicate_metadata":
        path = metadata
    header, rows = merge.read_table(path)
    if bad.startswith("duplicate_") and bad != "duplicate_header":
        rows.append(rows[0].copy())
    elif bad == "ambiguous_sequence":
        path = metadata
        header, rows = merge.read_table(path)
        rows[1]["NGS_LW_ID"] = rows[0]["NGS_LW_ID"]
    elif bad == "missing_sample":
        rows[0]["Sample"] = ""
    elif bad == "missing_key":
        path = metadata
        header, rows = merge.read_table(path)
        rows[0]["RekvNr"] = ""
    elif bad == "missing_source":
        header.remove("Sample")
        for row in rows:
            row.pop("Sample")
    elif bad == "empty_report":
        rows = []
    elif bad == "column_collision":
        header.append("NGS_Coverage_M")
        for row in rows:
            row["NGS_Coverage_M"] = "conflicting"
    save(path, (header, rows))
    if bad == "duplicate_header":
        path.write_text("Sample,Sample\nSEQ-001,SEQ-001\n")
    elif bad == "ragged_row":
        path.write_text("Sample,RunID\nSEQ-001,RUN,extra\n")
    original = master.read_bytes()
    with pytest.raises(merge.MergeError):
        merge.merge_files(master, metadata, report, "fastq")
    assert master.read_bytes() == original
    assert not list(tmp_path.glob("master.csv.bak.*"))
    assert not Path(str(master) + ".lock").exists()


def test_lock_blocks_another_writer(tmp_path):
    master, metadata, report = inputs(tmp_path)
    with merge.master_lock(master):
        with pytest.raises(merge.MergeError, match="locked"):
            merge.merge_files(master, metadata, report, "fastq")
        assert Path(str(master) + ".lock").is_dir()


def test_preview_cli_and_failure_exit(tmp_path):
    master, metadata, report = inputs(tmp_path)
    before = master.read_bytes()
    output = tmp_path / "preview.csv"
    command = [sys.executable, str(ROOT / "bin/merge_master_csv.py"), "--master", str(master),
               "--metadata", str(metadata), "--report", str(report), "--workflow", "fastq", "--output", str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert len(merge.read_table(output)[1]) == 3 and master.read_bytes() == before
    assert not list(tmp_path.glob("master.csv.bak.*"))
    assert subprocess.run(command, capture_output=True).returncode != 0


def test_xml_bindings_and_exact_base_schema():
    columns, canonical = merge.schema()
    assert len(columns) == 193
    assert len(set(columns)) == len(columns)
    for name, expected in (("FLUSEQ", ("Sample", "NGS_LW_ID")),
                           ("FLUSEQFASTA", ("OriginalName", "Key")),
                           ("LABWAREDIRECT", ("RekvNr", "Key"))):
        mapping = merge.read_mapping(merge.ASSETS / (name + ".xml"), canonical)
        assert mapping[expected[0]] == expected[1]
        # Every mapped XML field is either in the required schema or retained by name.
        assert all(dest in columns or source == dest for source, dest in mapping.items())


def test_unsupported_xml_transform_fails(tmp_path):
    path = tmp_path / "mapping.xml"
    path.write_text((merge.ASSETS / "FLUSEQ.xml").read_text().replace('IsRegex="0"', 'IsRegex="1"', 1))
    with pytest.raises(merge.MergeError, match="Unsupported transformation"):
        merge.read_mapping(path, merge.schema()[1])


def test_external_edit_detected_before_replace(tmp_path, monkeypatch):
    master, metadata, report = inputs(tmp_path)
    original_write = merge.write_table
    def write_and_edit(*args):
        original_write(*args)
        master.write_text("Key;manual\nchanged;yes\n")
    monkeypatch.setattr(merge, "write_table", write_and_edit)
    with pytest.raises(merge.MergeError, match="changed during merge"):
        merge.merge_files(master, metadata, report, "fastq")
    assert master.read_text() == "Key;manual\nchanged;yes\n"


def test_failed_replace_preserves_master_and_backup(tmp_path, monkeypatch):
    master, metadata, report = inputs(tmp_path)
    original = master.read_bytes()
    def fail(*args):
        raise OSError("simulated replace failure")
    monkeypatch.setattr(merge.os, "replace", fail)
    with pytest.raises(OSError):
        merge.merge_files(master, metadata, report, "fastq")
    assert master.read_bytes() == original
    assert next(tmp_path.glob("master.csv.bak.*")).read_bytes() == original
    assert not list(tmp_path.glob("*.tmp"))


class FakeSMB(update.SMB):
    def __init__(self, contents, failure=None):
        self.files = {"master.csv": contents}
        self.locked = False
        self.failure = failure
        self.commands = []

    def command(self, verb, *args):
        self.commands.append((verb, *args))
        if self.failure == verb:
            raise merge.MergeError("simulated SMB failure")
        if verb == "mkdir":
            assert not self.locked
            self.locked = True
        elif verb == "rmdir":
            self.locked = False
        elif verb == "get":
            Path(args[1]).write_bytes(self.files[args[0]])
        elif verb == "put":
            self.files[args[1]] = Path(args[0]).read_bytes()
            if self.failure == "corrupt_upload":
                self.files[args[1]] += b"corrupt"
        elif verb == "rename":
            assert self.locked and args[2] == "-f"
            self.files[args[1]] = self.files.pop(args[0])
        else:
            raise AssertionError(verb)


@pytest.mark.parametrize("failure", [None, "put", "corrupt_upload", "rename"])
def test_smb_transaction(tmp_path, failure):
    master, metadata, report = inputs(tmp_path)
    original = master.read_bytes()
    remote = FakeSMB(original, failure)
    config = {"master": None, "metadata": str(metadata), "master_smb": remote,
              "master_name": "master.csv", "nonempty_report": False}
    if failure:
        with pytest.raises(merge.MergeError):
            update.update(config, report, "fastq")
        assert remote.files["master.csv"] == original
    else:
        stats = update.update(config, report, "fastq")
        assert stats["added"] == 1 and stats["updated"] == 1
        assert remote.files[stats["backup"]] == original
        assert remote.files["master.csv"] != original
    assert not remote.locked


def test_smb_zero_exit_status_with_error_is_failure(monkeypatch):
    monkeypatch.setattr(update.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        [], 0, "NT_STATUS_ACCESS_DENIED opening remote file", ""))
    with pytest.raises(merge.MergeError, match="NT_STATUS_ACCESS_DENIED"):
        update.SMB("//server/share", "directory with spaces", "/unused").command("get", "a", "b")


@pytest.mark.parametrize("path", ['file; rm master.csv', '"quoted"', 'line\nbreak', '*', '?'])
def test_reject_smb_command_injection(path):
    with pytest.raises(merge.MergeError):
        update.smb_quote(path)


def test_smb_quoting_and_subprocess_arguments(monkeypatch):
    calls = []
    def run(args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(args, 0, "", "")
    monkeypatch.setattr(update.subprocess, "run", run)
    update.SMB("//server/share", "directory with spaces", "/auth").command("rename", "stage.csv", "master.csv", "-f")
    args, options = calls[0]
    assert args[-1] == 'rename "stage.csv" "master.csv" "-f"'
    assert args[args.index("-D") + 1] == "directory with spaces"
    assert not options.get("shell")


def test_wrapper_configuration_and_validation_skip(tmp_path):
    assert update.configuration({})[0] is None
    assert update.configuration({"FLUSEQ_MASTER_UPDATE": "true"}, "VER")[0] is None
    with pytest.raises(merge.MergeError, match="metadata share is not yet known"):
        update.configuration({"FLUSEQ_MASTER_UPDATE": "true"})
    master, metadata, _ = inputs(tmp_path)
    config, reason = update.configuration({"FLUSEQ_MASTER_CSV": str(master), "FLUSEQ_METADATA_CSV": str(metadata)})
    assert reason is None and config["master"] == str(master)


def test_wrapper_entrypoint_runs_merge_and_skips_validation(tmp_path):
    master, metadata, report = inputs(tmp_path)
    env = dict(os.environ, FLUSEQ_MASTER_UPDATE="true", FLUSEQ_MASTER_CSV=str(master), FLUSEQ_METADATA_CSV=str(metadata))
    command = [sys.executable, str(ROOT / "bin/update_master_csv.py"), "--workflow", "fastq", "--report", str(report)]
    original = master.read_bytes()
    assert subprocess.run(command + ["--check-config"], env=env, capture_output=True).returncode == 0
    assert master.read_bytes() == original
    result = subprocess.run(command + ["--validation", "VER"], env=env, capture_output=True, text=True)
    assert result.returncode == 0 and "validation run" in result.stdout
    assert master.read_bytes() == original
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert '"added": 1' in result.stdout
    assert master.read_bytes() != original


def test_empty_master_and_fastq_match_from_master_without_metadata(tmp_path):
    master, metadata, report = inputs(tmp_path)
    empty = save(tmp_path / "empty.csv", table("Key"))
    stats = merge.merge_files(empty, metadata, report, "fastq")
    assert stats["added"] == 2
    # An existing sequencing-ID association can supply the key when a metadata
    # export lacks that sample. No fabricated key is introduced.
    empty_metadata = save(tmp_path / "empty_metadata.csv", table("RekvNr;NGS_LW_ID"))
    report = save(report, table("Sample;RunID", ["SEQ-001", "RUN03"]))
    stats = merge.merge_files(master, empty_metadata, report, "fastq")
    assert stats["without_metadata"] == 1 and stats["updated"] == 1


def test_master_aliases_and_all_xml_columns_keep_values():
    _, canonical = merge.schema()
    for name in ("FLUSEQ", "FLUSEQFASTA", "LABWAREDIRECT"):
        mapping = merge.read_mapping(merge.ASSETS / (name + ".xml"), canonical)
        for source, destination in mapping.items():
            headers, rows = merge.transform(table(source, ["preserved text"]), mapping, canonical, "test")
            assert headers == [destination] and rows == [{destination: "preserved text"}]
    headers, rows = merge.transform(table("PASIENT_KJONN;PROVE_MANAD", ["value", "01"]), {}, canonical, "master")
    assert headers == ["Pasient_Kjnn", "Prove_M nad"]
    assert rows[0]["Prove_M nad"] == "01"


def test_smb_detects_nonparticipating_writer(tmp_path):
    master, metadata, report = inputs(tmp_path)
    remote = FakeSMB(master.read_bytes())
    original_command = remote.command
    reads = 0
    def command(verb, *args):
        nonlocal reads
        if verb == "get" and args[0] == "master.csv":
            reads += 1
            if reads == 2:
                remote.files["master.csv"] = b"Key;manual\nchanged;yes\n"
        original_command(verb, *args)
    remote.command = command
    config = {"master": None, "metadata": str(metadata), "master_smb": remote,
              "master_name": "master.csv", "nonempty_report": False}
    with pytest.raises(merge.MergeError, match="Remote master changed"):
        update.update(config, report, "fastq")
    assert remote.files["master.csv"] == b"Key;manual\nchanged;yes\n"
    assert not any(command[0] == "rename" for command in remote.commands)
    assert not remote.locked


def test_smb_metadata_download_with_mounted_master(tmp_path):
    master, metadata, report = inputs(tmp_path)
    remote = FakeSMB(b"unused")
    remote.files["metadata.csv"] = metadata.read_bytes()
    config = {"master": str(master), "metadata": None, "metadata_smb": remote,
              "metadata_name": "metadata.csv", "nonempty_report": False}
    assert update.update(config, report, "fastq")["added"] == 1


def test_backups_preserve_master_permissions(tmp_path):
    master, metadata, report = inputs(tmp_path)
    master.chmod(0o640)
    stats = merge.merge_files(master, metadata, report, "fastq")
    assert master.stat().st_mode & 0o777 == 0o640
    assert Path(stats["backup"]).stat().st_mode & 0o777 == 0o640


def test_new_fastq_without_metadata_has_blank_metadata_and_is_idempotent(tmp_path):
    master, metadata, report = inputs(tmp_path)
    report = save(report, table("Sample;RunID;New_report_column", ["SEQ-NEW", "RUN03", "0010"]))
    stats = merge.merge_files(master, metadata, report, "fastq")
    assert stats["added"] == 1 and stats["without_metadata"] == 1
    headers, rows = merge.read_table(master)
    assert len(rows) == 3
    added = rows[-1]
    assert added["Key"] == added["NGS_LW_ID"] == "SEQ-NEW"
    assert added[merge.KEY_SOURCE] == "Sample"
    assert added["NGS_Run_ID"] == "RUN03" and added["New_report_column"] == "0010"
    for column in ("Metadata_Key", "INF_Res", "Pasient_Kjnn", "Prove_M nad", "Added_metadata"):
        assert added[column] == ""
    before = master.read_bytes()
    repeated = merge.merge_files(master, metadata, report, "fastq")
    assert repeated["unchanged"] == 1 and not repeated["written"]
    assert master.read_bytes() == before


@pytest.mark.parametrize("workflow", ["fastq", "fasta"])
def test_later_metadata_enriches_provisional_row_without_duplicate(tmp_path, workflow):
    master = save(tmp_path / "master.csv", table("Key"))
    metadata = save(tmp_path / "metadata.csv", table("RekvNr;NGS_LW_ID;INF_Res"))
    report = save(tmp_path / "report.csv", table("Sample;RunID;DEPTH_HA", ["SEQ-NEW", "RUN03", "100"]))
    merge.merge_files(master, metadata, report, "fastq")
    headers, rows = merge.read_table(master)
    rows[0]["Comment_GC_Tessy"] = "retained comment"
    save(master, (headers, rows))
    save(metadata, table("RekvNr;NGS_LW_ID;INF_Res", ["0003", "SEQ-NEW", "new metadata"]))
    identifier = "Sample" if workflow == "fastq" else "OriginalName"
    value = "SEQ-NEW" if workflow == "fastq" else "0003"
    save(report, table(identifier + ";RunID", [value, "RUN04"]))
    stats = merge.merge_files(master, metadata, report, workflow)
    row, = merge.read_table(master)[1]
    assert stats["updated"] == 1 and stats["added"] == stats["without_metadata"] == 0
    assert row["Key"] == "0003" and row["NGS_LW_ID"] == "SEQ-NEW"
    assert row[merge.KEY_SOURCE] == "RekvNr"
    assert row["INF_Res"] == "new metadata" and row["NGS_Run_ID"] == "RUN04"
    assert row["Comment_GC_Tessy"] == "retained comment" and row["NGS_Reads_HA"] == "100"
    again = merge.merge_files(master, metadata, report, workflow)
    assert again["unchanged"] == 1 and not again["written"]


def test_unreported_provisional_sample_is_not_changed(tmp_path):
    master, metadata, report = inputs(tmp_path)
    save(report, table("Sample;RunID", ["SEQ-NEW", "RUN03"]))
    merge.merge_files(master, metadata, report, "fastq")
    before = merge.read_table(master)[1][-1]
    save(metadata, table("RekvNr;NGS_LW_ID;INF_Res", ["0003", "SEQ-NEW", "new metadata"]))
    save(report, table("Sample;RunID", ["SEQ-001", "RUN04"]))
    merge.merge_files(master, metadata, report, "fastq")
    assert merge.read_table(master)[1][-1] == before


@pytest.mark.parametrize("collision", ["master", "metadata", "resolved_key"])
def test_provisional_key_collisions_preserve_master(tmp_path, collision):
    master, metadata, report = inputs(tmp_path)
    if collision == "master":
        save(report, table("Sample;RunID", ["001", "RUN03"]))
    elif collision == "metadata":
        save(report, table("Sample;RunID", ["004", "RUN03"]))
    else:
        save(report, table("Sample;RunID", ["SEQ-NEW", "RUN03"]))
        merge.merge_files(master, metadata, report, "fastq")
        save(metadata, table("RekvNr;NGS_LW_ID", ["002", "SEQ-NEW"]))
    before = master.read_bytes()
    with pytest.raises(merge.MergeError, match="collides|existing master Key"):
        merge.merge_files(master, metadata, report, "fastq")
    assert master.read_bytes() == before


def test_missing_metadata_preserves_previously_imported_values(tmp_path):
    master, metadata, report = inputs(tmp_path)
    merge.merge_files(master, metadata, report, "fastq")
    first = merge.read_table(master)[1][0]
    save(metadata, table("RekvNr;NGS_LW_ID;INF_Res"))
    save(report, table("Sample;RunID", ["SEQ-001", "RUN04"]))
    stats = merge.merge_files(master, metadata, report, "fastq")
    row = merge.read_table(master)[1][0]
    assert stats["without_metadata"] == 1 and row["NGS_Run_ID"] == "RUN04"
    for column in ("Metadata_Key", "INF_Res", "Pasient_Kjnn", "Added_metadata"):
        assert row[column] == first[column]


def test_provisional_key_does_not_match_unrelated_metadata_key(tmp_path):
    master, metadata, report = inputs(tmp_path)
    save(report, table("Sample;RunID", ["SEQ-NEW", "RUN03"]))
    merge.merge_files(master, metadata, report, "fastq")
    save(metadata, table("RekvNr;NGS_LW_ID;INF_Res", ["SEQ-NEW", "OTHER-SEQUENCE", "unrelated"]))
    stats = merge.merge_files(master, metadata, report, "fastq")
    row = merge.read_table(master)[1][-1]
    assert stats["without_metadata"] == 1 and row["INF_Res"] == ""
    before = master.read_bytes()
    save(report, table("Sample;RunID", ["OTHER-SEQUENCE", "RUN04"]))
    with pytest.raises(merge.MergeError, match="unrelated provisional Key"):
        merge.merge_files(master, metadata, report, "fastq")
    assert master.read_bytes() == before

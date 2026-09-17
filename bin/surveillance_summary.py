#!/usr/bin/env python3
"""Create conservative, auditable influenza surveillance summary outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

SEGMENTS = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
EXPECTED_LENGTH = {
    "PB2": 2280,
    "PB1": 2274,
    "PA": 2151,
    "HA": 1701,
    "NP": 1497,
    "NA": 1410,
    "M": 982,
    "NS": 838,
}
SEGMENT_ALIASES = {"MP": "M", "M1": "M", "M2": "M"}
AMBIGUOUS_MIXED = set("RYSWKMBDHV")
EMPTY_VALUES = {"", "na", "nan", "none", "no matching mutations found", "not found"}
SCHEMA_VERSION = "1.0.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=("human-fastq", "human-fasta", "avian-fastq", "avian-fasta"))
    parser.add_argument("--fasta", nargs="*", default=[])
    parser.add_argument("--coverage", nargs="*", default=[])
    parser.add_argument("--subtype", nargs="*", default=[])
    parser.add_argument("--subtype-hits", nargs="*", default=[])
    parser.add_argument("--reassortment", nargs="*", default=[])
    parser.add_argument("--resistance", nargs="*", default=[])
    parser.add_argument("--resistance-databases", nargs="*", default=[])
    parser.add_argument("--mixed-min-count", type=int, default=5)
    parser.add_argument("--mixed-min-fraction", type=float, default=0.01)
    parser.add_argument("--outdir", type=Path, default=Path("."))
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def canonical_segment(value: str) -> str | None:
    upper = str(value).upper()
    tokens = [token for token in re.split(r"[^A-Z0-9]+", upper) if token]
    for token in tokens:
        token = SEGMENT_ALIASES.get(token, token)
        if token in SEGMENTS:
            return token
    match = re.search(r"(?:^|[_|.-])([1-8])(?:$|[_|.-])", upper)
    if match:
        return dict(zip("12345678", SEGMENTS))[match.group(1)]
    return None


def infer_sample(value: str) -> str:
    text = Path(str(value)).name
    text = re.sub(r"\.(?:fa|fasta|fna)(?:\.gz)?$", "", text, flags=re.I)
    if "|" in text:
        return text.split("|", 1)[0]
    match = re.search(r"(?i)(?:[_-])(?:PB2|PB1|PA|HA|NP|NA|MP|M|NS)(?:[_|.-]|$)", text)
    if match:
        return text[: match.start()]
    return re.sub(r"(?i)_(?:flumut|coverage|genin|genotyping).*$", "", text)


def read_fasta(paths: Iterable[str]) -> dict[tuple[str, str], str]:
    records: dict[tuple[str, str], str] = {}
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            continue
        header = None
        chunks: list[str] = []
        with open_text(path) as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        store_sequence(records, header, "".join(chunks), path)
                    header, chunks = line[1:].split()[0], []
                elif header is not None:
                    chunks.append(line)
        if header is not None:
            store_sequence(records, header, "".join(chunks), path)
    return records


def store_sequence(records: dict[tuple[str, str], str], header: str, sequence: str, path: Path) -> None:
    segment = canonical_segment(header) or canonical_segment(path.name)
    if not segment:
        return
    sample = infer_sample(header) or infer_sample(path.name)
    sequence = re.sub(r"\s+", "", sequence).upper()
    key = (sample, segment)
    if len(sequence) > len(records.get(key, "")):
        records[key] = sequence


def read_csv_rows(paths: Iterable[str], delimiter: str | None = None):
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            sep = delimiter
            if sep is None:
                try:
                    sep = csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
                except csv.Error:
                    sep = "\t" if "\t" in sample else ","
            for row in csv.DictReader(handle, delimiter=sep):
                yield path, {str(k).strip(): ("" if v is None else str(v).strip()) for k, v in row.items()}


def parse_coverage(paths: Iterable[str]) -> dict[tuple[str, str], float]:
    values: dict[tuple[str, str], float] = {}
    for path, row in read_csv_rows(paths):
        sample = row.get("Sample") or row.get("sample") or infer_sample(path.name)
        for column, raw in row.items():
            segment = canonical_segment(column)
            if not segment or not raw:
                continue
            try:
                value = float(raw)
            except ValueError:
                continue
            values[(sample, segment)] = max(values.get((sample, segment), 0.0), value)
    return values


def parse_subtypes(paths: Iterable[str]) -> dict[str, str]:
    calls: dict[str, str] = {}
    for path, row in read_csv_rows(paths):
        sample = row.get("Sample") or row.get("sample") or infer_sample(path.name)
        subtype = row.get("Subtype") or row.get("subtype") or ""
        if sample:
            calls[sample] = subtype
    return calls


def subtype_target(path: Path) -> str | None:
    name = path.name.lower()
    if name.startswith("ha_") or "_ha_" in name:
        return "HA"
    if name.startswith("na_") or "_na_" in name:
        return "NA"
    return None


def subtype_sample(path: Path) -> str:
    name = path.name
    name = re.sub(r"(?i)^(?:ha|na)_", "", name)
    name = re.sub(r"(?i)_(?:sorted|filtered).*", "", name)
    return re.sub(r"\.(?:tsv|txt)$", "", name, flags=re.I)


def parse_subtype_confidence(paths: Iterable[str]):
    hits: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for raw_path in paths:
        path = Path(raw_path)
        target = subtype_target(path)
        if not target or not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                try:
                    hit = {
                        "subject": fields[1],
                        "identity": float(fields[2]),
                        "length": int(float(fields[3])),
                        "bitscore": float(fields[11]),
                    }
                except ValueError:
                    continue
                hits[(subtype_sample(path), target)].append(hit)

    confidence: dict[tuple[str, str], dict[str, object]] = {}
    competing: dict[str, list[str]] = defaultdict(list)
    for key, entries in hits.items():
        unique = {}
        for hit in entries:
            subject = str(hit["subject"])
            if subject not in unique or float(hit["bitscore"]) > float(unique[subject]["bitscore"]):
                unique[subject] = hit
        ranked = sorted(
            unique.values(), key=lambda item: (float(item["bitscore"]), float(item["identity"])), reverse=True
        )
        top = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        separation = 1.0
        if second and float(top["bitscore"]) > 0:
            separation = max(0.0, 1.0 - float(second["bitscore"]) / float(top["bitscore"]))
            if float(second["bitscore"]) / float(top["bitscore"]) >= 0.95:
                competing[key[0]].append(f"{key[1]} competing subtype hits")
        coverage = min(1.0, int(top["length"]) / EXPECTED_LENGTH[key[1]])
        score = (float(top["identity"]) / 100.0) * coverage * (0.5 + 0.5 * separation)
        confidence[key] = {
            "subject": top["subject"],
            "score": round(score, 4),
            "identity": top["identity"],
            "query_coverage": round(coverage, 4),
            "bitscore_separation": round(separation, 4),
            "level": "high" if score >= 0.9 else "moderate" if score >= 0.75 else "low",
        }
    return confidence, competing


def qc_status(sequence: str | None, coverage: float | None, segment: str) -> tuple[str, str]:
    if sequence is None:
        return "missing", "no segment consensus"
    n_fraction = sequence.count("N") / len(sequence) if sequence else 1.0
    length_fraction = len(sequence) / EXPECTED_LENGTH[segment]
    effective_coverage = coverage if coverage is not None else 100.0 * (1.0 - n_fraction) * length_fraction
    if length_fraction >= 0.8 and n_fraction <= 0.05 and effective_coverage >= 80:
        return "pass", "length, ambiguity and coverage thresholds passed"
    if length_fraction >= 0.5 and n_fraction <= 0.2 and effective_coverage >= 50:
        return "warning", "partial or ambiguous segment consensus"
    return "fail", "insufficient segment consensus"


def database_metadata(paths: Iterable[str]) -> list[dict[str, str]]:
    metadata = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        version_match = re.search(r"(?<!\d)(20\d{2}(?:[-_]20?\d{2})?|\d{4})(?!\d)", path.stem)
        metadata.append(
            {
                "name": path.name,
                "version": version_match.group(1).replace("_", "-") if version_match else "unversioned",
                "sha256": digest.hexdigest(),
            }
        )
    return metadata


def resistance_findings(paths: Iterable[str]):
    findings = []
    samples = set()
    for path, row in read_csv_rows(paths):
        sample = row.get("Sample") or row.get("sample") or infer_sample(path.name)
        if not sample:
            continue
        samples.add(sample)
        for column, value in row.items():
            lowered = column.lower()
            if column.lower() == "sample" or not any(
                term in lowered for term in ("mutation", "flumut", "inhib", "resistance")
            ):
                continue
            if value.strip().lower() in EMPTY_VALUES:
                continue
            segment = canonical_segment(column) or canonical_segment(path.name) or "unknown"
            for finding in (item.strip() for item in re.split(r"[;|]", value) if item.strip()):
                findings.append(
                    {
                        "sample": sample,
                        "segment": segment,
                        "finding": finding,
                        "source_column": column,
                        "status": "detected",
                    }
                )
    return findings, samples


def parse_reassortment(paths: Iterable[str]):
    results = {}
    for path, row in read_csv_rows(paths):
        sample = row.get("Sample") or row.get("sample") or infer_sample(path.name)
        if not sample:
            continue
        status = row.get("Reassortment") or row.get("reassortment") or "Unknown"
        segment_calls = {
            segment: row.get(segment) or row.get("MP" if segment == "M" else segment) or "Missing"
            for segment in SEGMENTS
        }
        lineages = {
            re.sub(r"\([^)]*\)$", "", value)
            for value in segment_calls.values()
            if value and value != "Missing" and not value.startswith("TooLow")
        }
        results[sample] = {"status": status, "lineages": sorted(lineages), "segments": segment_calls}
    return results


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sequences = read_fasta(args.fasta)
    coverage = parse_coverage(args.coverage)
    subtype_calls = parse_subtypes(args.subtype)
    subtype_confidence, competing_hits = parse_subtype_confidence(args.subtype_hits)
    reassortment = parse_reassortment(args.reassortment)
    findings, resistance_samples = resistance_findings(args.resistance)
    databases = database_metadata(args.resistance_databases)

    samples = sorted(
        {sample for sample, _ in sequences}
        | {sample for sample, _ in coverage}
        | set(subtype_calls)
        | {sample for sample, _ in subtype_confidence}
        | set(reassortment)
        | resistance_samples
    )

    segment_rows = []
    mixed_flags = {}
    for sample in samples:
        mixed_evidence = list(competing_hits.get(sample, []))
        for segment in SEGMENTS:
            sequence = sequences.get((sample, segment))
            cov = coverage.get((sample, segment))
            ambiguous_count = sum(base in AMBIGUOUS_MIXED for base in sequence or "")
            ambiguous_fraction = ambiguous_count / len(sequence) if sequence else None
            if (
                sequence
                and ambiguous_count >= args.mixed_min_count
                and ambiguous_fraction is not None
                and ambiguous_fraction >= args.mixed_min_fraction
            ):
                mixed_evidence.append(f"{segment} has {ambiguous_count} mixed IUPAC sites ({ambiguous_fraction:.4f})")
            status, reason = qc_status(sequence, cov, segment)
            segment_rows.append(
                {
                    "sample": sample,
                    "segment": segment,
                    "qc_status": status,
                    "sequence_length": len(sequence) if sequence else 0,
                    "expected_length": EXPECTED_LENGTH[segment],
                    "coverage_percent": "" if cov is None else f"{cov:.2f}",
                    "n_fraction": "" if not sequence else f"{sequence.count('N') / len(sequence):.6f}",
                    "mixed_iupac_count": ambiguous_count,
                    "mixed_iupac_fraction": "" if ambiguous_fraction is None else f"{ambiguous_fraction:.6f}",
                    "qc_reason": reason,
                }
            )
        mixed_flags[sample] = {
            "status": (
                "flagged"
                if mixed_evidence
                else ("not_detected" if any(key[0] == sample for key in sequences) else "unknown")
            ),
            "evidence": sorted(set(mixed_evidence)),
        }

    write_tsv(
        args.outdir / "segment_qc.tsv",
        [
            "sample",
            "segment",
            "qc_status",
            "sequence_length",
            "expected_length",
            "coverage_percent",
            "n_fraction",
            "mixed_iupac_count",
            "mixed_iupac_fraction",
            "qc_reason",
        ],
        segment_rows,
    )

    reassortment_rows = []
    for sample in samples:
        result = reassortment.get(
            sample, {"status": "Unknown", "lineages": [], "segments": {segment: "Missing" for segment in SEGMENTS}}
        )
        reassortment_rows.append(
            {
                "sample": sample,
                "reassortment_status": result["status"],
                "lineage_count": len(result["lineages"]),
                "lineages": ";".join(result["lineages"]),
                **result["segments"],
                "method": "best BLAST identity per segment; review required",
            }
        )
    write_tsv(
        args.outdir / "reassortment_summary.tsv",
        ["sample", "reassortment_status", "lineage_count", "lineages", *SEGMENTS, "method"],
        reassortment_rows,
    )

    db_names = ";".join(item["name"] for item in databases) or "not_provided"
    db_versions = ";".join(item["version"] for item in databases) or "unknown"
    db_hashes = ";".join(item["sha256"] for item in databases) or "unknown"
    if not findings:
        findings = [
            {
                "sample": sample,
                "segment": "all",
                "finding": "",
                "source_column": "",
                "status": "none_detected",
            }
            for sample in samples
        ]
    resistance_rows = [
        {
            **finding,
            "database_name": db_names,
            "database_version": db_versions,
            "database_sha256": db_hashes,
        }
        for finding in findings
    ]
    write_tsv(
        args.outdir / "resistance_summary.tsv",
        [
            "sample",
            "segment",
            "status",
            "finding",
            "source_column",
            "database_name",
            "database_version",
            "database_sha256",
        ],
        resistance_rows,
    )

    qc_by_sample = defaultdict(dict)
    for row in segment_rows:
        qc_by_sample[row["sample"]][row["segment"]] = row["qc_status"]

    sample_summaries = []
    for sample in samples:
        confidence = {
            target: subtype_confidence.get((sample, target), {"score": None, "level": "unknown", "subject": None})
            for target in ("HA", "NA")
        }
        sample_summaries.append(
            {
                "sample": sample,
                "mode": args.mode,
                "subtype": subtype_calls.get(sample, ""),
                "subtype_confidence": confidence,
                "segment_qc": qc_by_sample[sample],
                "mixed_infection": mixed_flags[sample],
                "reassortment": reassortment.get(sample, {"status": "Unknown", "lineages": []}),
                "resistance_finding_count": sum(
                    row["sample"] == sample and row["status"] == "detected" for row in resistance_rows
                ),
            }
        )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "interpretation": "screening output for public-health surveillance; flagged results require review",
        "methods": {
            "segments": list(SEGMENTS),
            "qc": {
                "pass": "length >=80% expected, N <=5%, effective coverage >=80%",
                "warning": "length >=50% expected, N <=20%, effective coverage >=50%",
            },
            "mixed_infection": {
                "minimum_mixed_iupac_sites": args.mixed_min_count,
                "minimum_mixed_iupac_fraction": args.mixed_min_fraction,
                "competing_subtype_hit_ratio": 0.95,
            },
            "subtype_confidence": "identity x query coverage x (0.5 + 0.5 x top-hit bitscore separation); ranking score, not a calibrated probability",
        },
        "resistance_databases": databases,
        "samples": sample_summaries,
    }
    with (args.outdir / "sample_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()

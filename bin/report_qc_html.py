#!/usr/bin/env python3
"""Render a human FASTQ report as an offline QC overview; Python standard library only."""

import argparse
import csv
import hashlib
import html
import io
import json
import math
from collections import Counter
from pathlib import Path
from statistics import median

VERSION = "1.0.0"
SEGMENTS = ("HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2")
PROTEINS = ("HA1", "HA2", "NA", "M1", "M2", "NP", "NS", "PA", "PB1", "PB2")
MISSING = {"", "NA", "NAN", "NONE", "NULL", "N/A", "-"}
COVERAGE_THRESHOLD = 80
SUBTYPES = {
    "H3N2": "A/H3N2",
    "H1N1": "A/H1N1",
    "VIC": "B/Victoria",
    "VICVIC": "B/Victoria",
    "YAM": "B/Yamagata",
    "YAMYAM": "B/Yamagata",
}


def clean(value):
    value = (value or "").strip()
    return "" if value.upper() in MISSING else value


def number(value, maximum=None):
    try:
        result = float(clean(value))
    except ValueError:
        return None
    if not math.isfinite(result) or result < 0 or (maximum is not None and result > maximum):
        return None
    return result


def coverage(row, segment):
    aliases = ("Coverage-MP", "Coverage-M") if segment == "MP" else (f"Coverage-{segment}",)
    for column in aliases:
        if clean(row.get(column)):
            return number(row[column], 100)
    return None


def counts(values):
    return dict(sorted(Counter(values).items(), key=lambda item: (-item[1], item[0])))


def read_csv(path):
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    fields = reader.fieldnames
    if not fields or "Sample" not in fields:
        raise ValueError("The input CSV must contain a Sample column.")
    if len(set(fields)) != len(fields):
        raise ValueError("The input CSV contains duplicate column names.")
    rows = []
    identifiers = set()
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(
                f"Malformed CSV record ending at line {reader.line_num}: column count differs from header."
            )
        sample = clean(row["Sample"])
        if not sample:
            raise ValueError(f"Missing Sample at CSV line {reader.line_num}.")
        if sample in identifiers:
            raise ValueError(f"Duplicate Sample at CSV line {reader.line_num}; one row per sample is required.")
        identifiers.add(sample)
        rows.append(row)
    return fields, rows, hashlib.sha256(raw).hexdigest()


def summarise_sample(row, fields):
    cov = {seg: coverage(row, seg) for seg in SEGMENTS}
    raw_qc = clean(row.get("NGS_QC_Sum"))
    # Only an explicitly empty summary is a known no-flag result. NA is unavailable.
    qc_available = "NGS_QC_Sum" in fields and (not row["NGS_QC_Sum"].strip() or bool(raw_qc))
    reported_qc = "Review" if raw_qc else ("No flags" if qc_available else "Unavailable")
    reasons = []
    if raw_qc:
        reasons.append(f"Reported QC: {raw_qc}")
    if not qc_available:
        reasons.append("NGS QC summary unavailable")
    low = [seg for seg, value in cov.items() if value is not None and value < COVERAGE_THRESHOLD]
    missing = [seg for seg, value in cov.items() if value is None]
    if low:
        reasons.append("Coverage below 80%: " + ", ".join(low))
    if missing:
        reasons.append("Coverage missing or invalid: " + ", ".join(missing))
    nextclade = {protein: clean(row.get(f"Nextclade QC {protein}")) for protein in PROTEINS}
    concerning = [f"{protein}: {value}" for protein, value in nextclade.items() if value and value.lower() != "good"]
    absent = [protein for protein, value in nextclade.items() if not value]
    if concerning:
        reasons.append("Nextclade QC: " + "; ".join(concerning))
    if absent:
        reasons.append("Nextclade QC unavailable: " + ", ".join(absent))
    gisaid = clean(row.get("GISAID_Comment"))
    if gisaid and gisaid.lower() != "no flags":
        reasons.append("Reported GISAID comment: " + gisaid)
    status = "Review" if reasons else "No flags"
    if not qc_available and not raw_qc and all(value is None for value in cov.values()):
        status = "Unavailable"

    subtype_raw = clean(row.get("Subtype"))
    subtype = SUBTYPES.get(subtype_raw, subtype_raw) or "Unavailable"
    subclade = clean(row.get("Subclade_Nomenclature_Subclade")) or clean(row.get("subclade"))
    fraction = number(row.get("Subclade_Nomenclature_Subclade_Match_Fraction"), 1)
    char_fraction = number(row.get("Characterisation_Subclade_Match_Fraction"), 1)
    char_status = clean(row.get("Characterisation_Status"))
    provisional = bool(
        (fraction is not None and fraction < 1)
        or (char_fraction is not None and char_fraction < 1)
        or clean(row.get("Subclade_Nomenclature_Closest_Subclade_Missing_Mutations"))
        or clean(row.get("Characterisation_Missing_Subclade_Mutations"))
        or "review" in char_status.lower()
        or clean(row.get("Characterisation_Result")).lower().startswith("provisional:")
    )
    if subclade:
        call_status = "Provisional" if provisional else ("Exact rule match" if fraction == 1 else "Match not assessed")
    else:
        call_status = "Unavailable"
    conclusion = clean(row.get("Conclusion"))
    resistance = {
        key.removeprefix("DR_Res_"): clean(value) or "Unavailable"
        for key, value in row.items()
        if key.startswith("DR_Res_")
    }
    return {
        "sample": clean(row["Sample"]),
        "barcode": clean(row.get("Barcode")),
        "sequence_id": clean(row.get("SequenceID")),
        "plate_position": clean(row.get("PCR-PlatePosition")),
        "qc_status": status,
        "reported_qc_status": reported_qc,
        "reported_qc": raw_qc,
        "review_reasons": reasons,
        "coverage": cov,
        "segment_reads": {seg: number(row.get(f"DEPTH_{seg}")) for seg in SEGMENTS},
        "nextclade_qc": nextclade,
        "missing_nextclade_qc": absent,
        "subtype": subtype,
        "raw_subtype": subtype_raw,
        "sequence_result": clean(row.get("Sekvens_Resultat")),
        "subclade": subclade,
        "subclade_status": call_status,
        "subclade_match_fraction": fraction,
        "characterisation_status": char_status,
        "reference_virus": clean(row.get("Characterisation_Reference_Virus")),
        "characterisation_provisional": provisional,
        "reassortment_conclusion": conclusion,
        "resistance_codes": resistance,
        "irma": {key.removeprefix("IRMA_"): number(value) for key, value in row.items() if key.startswith("IRMA_")},
    }


def build_report(path):
    fields, rows, digest = read_csv(path)
    samples = [summarise_sample(row, fields) for row in rows]
    metadata = {
        column: sorted({clean(row.get(column)) for row in rows if clean(row.get(column))})
        for column in ("RunID", "Instrument ID", "Date", "Release Version")
    }
    warnings = []
    for column, values in metadata.items():
        if len(values) > 1:
            warnings.append(f"Multiple values in {column}: " + "; ".join(values))
        missing = sum(not clean(row.get(column)) for row in rows)
        if missing:
            warnings.append(f"{column} unavailable for {missing}/{len(rows)} samples.")
    if "NGS_QC_Sum" not in fields:
        warnings.append("NGS_QC_Sum column is absent; reported QC is unavailable.")
    if not samples:
        warnings.append("The CSV contains no sample records.")

    segment_summary = {}
    for seg in SEGMENTS:
        values = [s["coverage"][seg] for s in samples if s["coverage"][seg] is not None]
        reads = [s["segment_reads"][seg] for s in samples if s["segment_reads"][seg] is not None]
        segment_summary[seg] = {
            "at_least_80": sum(value >= COVERAGE_THRESHOLD for value in values),
            "below_80": sum(value < COVERAGE_THRESHOLD for value in values),
            "missing_or_invalid": len(samples) - len(values),
            "median_coverage": median(values) if values else None,
            "median_segment_reads": median(reads) if reads else None,
            "segment_reads_available": len(reads),
        }
    reported = counts(s["reported_qc_status"] for s in samples)
    assessed = counts(s["qc_status"] for s in samples)
    status = (
        "Review required" if warnings or any(s["qc_status"] != "No flags" for s in samples) else "No QC flags detected"
    )
    if not samples:
        status = "Not assessable"
    subtype_counts = counts(s["subtype"] for s in samples)
    provisional = sum(s["characterisation_provisional"] for s in samples)
    text = [
        f"{len(samples)} samples are represented in the CSV: {reported.get('No flags', 0)} have no reported NGS QC flags, "
        f"{reported.get('Review', 0)} have reported flags, and {reported.get('Unavailable', 0)} have unavailable NGS QC summaries.",
        "Reported subtypes: "
        + (", ".join(f"{key} ({value})" for key, value in subtype_counts.items()) or "none")
        + ".",
        f"{provisional} samples have provisional or incomplete characterisation evidence. "
        f"{sum(bool(s['missing_nextclade_qc']) for s in samples)} samples lack at least one expected Nextclade QC value.",
    ]
    return {
        "schema_version": "1.0",
        "generator_version": VERSION,
        "source": {"filename": path.name, "sha256": digest},
        "metadata": metadata,
        "assessment": status,
        "warnings": warnings,
        "summary_text": text,
        "rules": {
            "coverage_threshold_percent": COVERAGE_THRESHOLD,
            "reported_qc": "NGS_QC_Sum from the source CSV; no mutation calls are recalculated",
            "review": "Reported flags, unavailable NGS summary, low/missing/invalid coverage, "
            "non-good/missing Nextclade QC, or a reported GISAID comment",
            "run_acceptance": "Control outcomes, expected sample count and laboratory acceptance rules are not supplied",
        },
        "counts": {
            "samples": len(samples),
            "reported_qc": reported,
            "assessment": assessed,
            "subtypes": subtype_counts,
            "provisional_characterisation": provisional,
            "characterisation_status": counts(s["characterisation_status"] or "Unavailable" for s in samples),
            "subclades": counts(
                f"{s['subtype']} · {s['subclade'] or 'Unavailable'} · {s['subclade_status']}" for s in samples
            ),
        },
        "segments": segment_summary,
        "samples": samples,
    }


def esc(value):
    return html.escape(str(value), quote=True)


def fmt(value, suffix="", precision=2):
    return "—" if value is None else f"{value:,.{precision}f}".rstrip("0").rstrip(".") + suffix


def badge(value):
    color = {"No flags": "good", "Review": "review", "Unavailable": "missing"}.get(value, "missing")
    return f'<span class="badge {color}">{esc(value)}</span>'


def distribution(items):
    total = sum(items.values()) or 1
    return (
        "".join(
            f'<div class="distribution"><span>{esc(key)}</span><b>{value}</b>'
            f'<progress value="{value}" max="{total}">{value}/{total}</progress></div>'
            for key, value in items.items()
        )
        or "<p>No records.</p>"
    )


def render(report, template):
    samples = report["samples"]
    run = ", ".join(report["metadata"]["RunID"]) or Path(report["source"]["filename"]).stem
    rows = []
    for sample in samples:
        cells = []
        for seg in SEGMENTS:
            value = sample["coverage"][seg]
            level = "missing" if value is None else ("good" if value >= COVERAGE_THRESHOLD else "review")
            title = (
                f"{seg}: reported coverage {fmt(value, '%')}; IRMA segment reads {fmt(sample['segment_reads'][seg])}"
            )
            cells.append(
                f'<td class="metric {level}" title="{esc(title)}" data-coverage="{esc(fmt(value, "%"))}" '
                f'data-reads="{esc(fmt(sample["segment_reads"][seg]))}">{esc(fmt(value, "%"))}</td>'
            )
        reasons = "".join(f"<li>{esc(reason)}</li>" for reason in sample["review_reasons"])
        info = {
            "Reported NGS QC": sample["reported_qc"] or sample["reported_qc_status"],
            "Sequence result": sample["sequence_result"] or "Unavailable",
            "Barcode / plate": " / ".join(filter(None, [sample["barcode"], sample["plate_position"]])) or "Unavailable",
            "Sequence ID": sample["sequence_id"] or "Unavailable",
            "Characterisation": sample["characterisation_status"] or "Unavailable",
            "Reference-virus category": sample["reference_virus"] or "Unavailable",
            "Subclade evidence": sample["subclade_status"],
            "Reassortment screen (reported)": sample["reassortment_conclusion"] or "Unavailable",
            "Drug lookup codes (reported)": "; ".join(
                f"{key}: {value}" for key, value in sample["resistance_codes"].items()
            )
            or "Unavailable",
            "IRMA read flow": "; ".join(f"{key}: {fmt(value, precision=5)}" for key, value in sample["irma"].items())
            or "Unavailable",
        }
        detail = "".join(f"<dt>{esc(key)}</dt><dd>{esc(value)}</dd>" for key, value in info.items())
        subclade = esc(sample["subclade"] or "Unavailable")
        if sample["subclade"]:
            subclade += f'<small>{esc(sample["subclade_status"])}</small>'
        rows.append(
            f'<tr data-qc="{esc(sample["qc_status"])}" data-subtype="{esc(sample["subtype"])}">'
            f'<th scope="row">{esc(sample["sample"])}</th><td>{badge(sample["qc_status"])}</td>'
            f'<td>{esc(sample["subtype"])}</td><td>{subclade}</td>{"".join(cells)}'
            f'<td><details><summary>Details</summary><div class="sample-details">'
            f'<strong>{esc(sample["sample"])}</strong><ul>{reasons}</ul><dl>{detail}</dl></div></details></td></tr>'
        )
    segment_rows = "".join(
        f'<tr><th scope="row">{seg}</th><td>{values["at_least_80"]}</td><td>{values["below_80"]}</td>'
        f'<td>{values["missing_or_invalid"]}</td><td>{fmt(values["median_coverage"], "%")}</td>'
        f'<td>{fmt(values["median_segment_reads"])}</td><td>{values["segment_reads_available"]}/{len(samples)}</td></tr>'
        for seg, values in report["segments"].items()
    )
    reported = report["counts"]["reported_qc"]
    replacements = {
        "TITLE": esc(run),
        "ASSESSMENT": esc(report["assessment"]),
        "RUN_META": esc(
            " · ".join(
                " / ".join(report["metadata"][key]) or f"{key} unavailable"
                for key in ("Instrument ID", "Date", "Release Version")
            )
        ),
        "TOTAL": str(len(samples)),
        "NO_FLAGS": str(reported.get("No flags", 0)),
        "FLAGGED": str(reported.get("Review", 0)),
        "UNAVAILABLE": str(reported.get("Unavailable", 0)),
        "SUMMARY": "".join(f"<li>{esc(line)}</li>" for line in report["summary_text"]),
        "WARNINGS": "".join(f'<p class="notice">{esc(value)}</p>' for value in report["warnings"]),
        "SUBTYPES": distribution(report["counts"]["subtypes"]),
        "CHARACTERISATION": distribution(report["counts"]["characterisation_status"]),
        "SUBCLADES": distribution(report["counts"]["subclades"]),
        "SEGMENT_ROWS": segment_rows,
        "SAMPLE_ROWS": "".join(rows),
        "SEGMENT_HEADERS": "".join(f'<th scope="col">{seg}</th>' for seg in SEGMENTS),
        "SUBTYPE_OPTIONS": "".join(
            f'<option value="{esc(key)}">{esc(key)}</option>' for key in report["counts"]["subtypes"]
        ),
        "ASSESSED_COUNTS": esc(
            ", ".join(f"{value} {key.lower()}" for key, value in report["counts"]["assessment"].items())
        ),
        "SOURCE": esc(report["source"]["filename"]),
        "SHA256": report["source"]["sha256"],
        "VERSION": VERSION,
    }
    # Substitute once, so user-supplied text resembling a placeholder stays literal.
    import re

    return re.sub(r"@@([A-Z_0-9]+)@@", lambda match: replacements[match[1]], template)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Final human FASTQ report CSV")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--template", type=Path, default=Path(__file__).parent / "templates" / "report_qc.html")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    try:
        report = build_report(args.input)
        page = render(report, args.template.read_text(encoding="utf-8"))
        args.outdir.mkdir(parents=True, exist_ok=True)
        prefix = args.outdir / f"{args.input.stem}_qc"
        Path(f"{prefix}.html").write_text(page, encoding="utf-8")
        Path(f"{prefix}.json").write_text(
            json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
        )
    except (OSError, UnicodeError, ValueError, csv.Error) as error:
        parser.exit(2, f"QC report error: {error}\n")
    print(f"Created {prefix}.html and {prefix}.json ({report['counts']['samples']} samples)")


if __name__ == "__main__":
    main()

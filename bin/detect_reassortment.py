#!/usr/bin/env python3
"""Summarise influenza segment BLAST hits and screen for reassortment.

The preferred subject identifier is:

    ORIGIN|SUBTYPE|STRAIN|SEGMENT|ACCESSION

SEGMENT may be omitted because the query identifier determines the segment.
Legacy STRAIN|LINEAGE|SEGMENT|ACCESSION, STRAIN|SEGMENT|ACCESSION, and
STRAIN|ACCESSION identifiers remain accepted, but unavailable metadata is
reported as UNKNOWN.

This is a similarity-based screen, not a phylogenetic reassortment analysis.
"""

import argparse
import re
from typing import Dict, Iterable, List

import pandas as pd


EXPECTED_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
IDENTITY_THRESHOLD = 80.0  # %
UNKNOWN = "UNKNOWN"

SEGMENT_ALIASES = {
    "M": "MP",
    "M1": "MP",
    "M2": "MP",
    "HA1": "HA",
    "HA2": "HA",
    "NS1": "NS",
    "NS2": "NS",
}

FALLBACK_SEGMENT_PATTERNS = [
    "PB2",
    "PB1",
    "PA",
    "HA",
    "NP",
    "NA",
    "MP",
    "NS",
    "M2",
    "M1",
    "M",
    "HA2",
    "HA1",
    "NS2",
    "NS1",
]


def canonicalize_segment(token) -> str:
    """Normalize observed segment labels to report schema labels."""
    if token is None or pd.isna(token):
        return pd.NA
    value = str(token).strip().upper()
    if not value:
        return pd.NA
    if value in EXPECTED_SEGMENTS:
        return value
    if value in SEGMENT_ALIASES:
        return SEGMENT_ALIASES[value]
    return pd.NA


def infer_segment_from_qseqid(qseqid) -> str:
    """Infer a segment from common query-identifier formats."""
    if qseqid is None or pd.isna(qseqid):
        return pd.NA

    raw = str(qseqid).strip()
    if not raw:
        return pd.NA
    upper = raw.upper()

    legacy = re.search(r"_([^_]+)$", upper)
    if legacy:
        segment = canonicalize_segment(legacy.group(1))
        if not pd.isna(segment):
            return segment

    tokens = [token for token in re.split(r"[^A-Z0-9]+", upper) if token]
    for token in tokens:
        segment = canonicalize_segment(token)
        if not pd.isna(segment):
            return segment

    for pattern in FALLBACK_SEGMENT_PATTERNS:
        if re.search(rf"(?<![A-Z0-9]){pattern}(?![A-Z0-9])", upper):
            segment = canonicalize_segment(pattern)
            if not pd.isna(segment):
                return segment

    return pd.NA


def normalize_origin(value: str) -> str:
    """Normalize common host/origin labels used by reference databases."""
    origin = value.strip().upper()
    aliases = {
        "BIRD": "AVIAN",
        "BIRDS": "AVIAN",
        "HOMO SAPIENS": "HUMAN",
        "PIG": "SWINE",
        "PIGS": "SWINE",
    }
    return aliases.get(origin, origin) if origin else UNKNOWN


def normalize_subtype(value: str) -> str:
    """Normalize common subtype and influenza B lineage labels."""
    subtype = value.strip()
    if not subtype:
        return UNKNOWN
    aliases = {
        "PDM09": "H1N1",
        "H1N1PDM": "H1N1",
        "H1N1PDM09": "H1N1",
        "VIC": "B/Victoria",
        "VICTORIA": "B/Victoria",
        "B/VIC": "B/Victoria",
        "B/VICTORIA": "B/Victoria",
        "YAM": "B/Yamagata",
        "YAMAGATA": "B/Yamagata",
        "B/YAM": "B/Yamagata",
        "B/YAMAGATA": "B/Yamagata",
    }
    return aliases.get(subtype.upper(), subtype.upper())


def looks_like_enriched_header(parts: List[str]) -> bool:
    """Distinguish ORIGIN|SUBTYPE|... from legacy STRAIN|LINEAGE|...."""
    if len(parts) < 3:
        return False
    subtype = normalize_subtype(parts[1])
    subtype_like = bool(
        re.fullmatch(r"H\d+N\d+(?:PDM09)?", subtype, flags=re.IGNORECASE)
        or subtype.lower() in {"b/victoria", "b/yamagata"}
        or subtype == UNKNOWN
    )
    strain_like = bool(re.match(r"^[A-D]/", parts[0], flags=re.IGNORECASE))
    return subtype_like and not strain_like


def parse_subject_id(subject_id) -> Dict[str, str]:
    """Parse enriched and legacy reassortment-database identifiers."""
    parsed = {
        "origin": UNKNOWN,
        "subtype": UNKNOWN,
        "strain": UNKNOWN,
        "subject_segment": UNKNOWN,
        "accession": UNKNOWN,
    }
    if subject_id is None or pd.isna(subject_id):
        return parsed

    parts = [part.strip() for part in str(subject_id).split("|")]
    if parts and parts[0]:
        parsed["strain"] = parts[0]

    if len(parts) >= 5 and looks_like_enriched_header(parts):
        segment = canonicalize_segment(parts[3])
        parsed.update(
            origin=normalize_origin(parts[0]),
            subtype=normalize_subtype(parts[1]),
            strain=parts[2] or UNKNOWN,
            subject_segment=segment if not pd.isna(segment) else UNKNOWN,
            accession=parts[4] or UNKNOWN,
        )
    elif len(parts) >= 4 and looks_like_enriched_header(parts):
        # Enriched header without subject segment; the query ID supplies it.
        parsed.update(
            origin=normalize_origin(parts[0]),
            subtype=normalize_subtype(parts[1]),
            strain=parts[2] or UNKNOWN,
            accession=parts[3] or UNKNOWN,
        )
    elif len(parts) == 4:
        # Legacy STRAIN|LINEAGE|SEGMENT|ACCESSION.
        segment = canonicalize_segment(parts[2])
        parsed.update(
            subtype=normalize_subtype(parts[1]),
            subject_segment=segment if not pd.isna(segment) else UNKNOWN,
            accession=parts[3] or UNKNOWN,
        )
    elif len(parts) == 3:
        # Legacy STRAIN|SEGMENT|ACCESSION.
        segment = canonicalize_segment(parts[1])
        parsed.update(
            subject_segment=segment if not pd.isna(segment) else UNKNOWN,
            accession=parts[2] or UNKNOWN,
        )
    elif len(parts) == 2:
        parsed["accession"] = parts[1] or UNKNOWN

    return parsed


def load_reference_metadata(path: str) -> Dict[str, Dict[str, str]]:
    """Load accession annotations used to enrich legacy subject identifiers."""
    if not path:
        return {}

    table = pd.read_csv(path, dtype=str).fillna("")
    required = {"accession", "origin", "subtype", "strain"}
    missing_columns = required.difference(table.columns)
    if missing_columns:
        raise ValueError(
            "Reference metadata is missing column(s): "
            + ", ".join(sorted(missing_columns))
        )
    if table["accession"].duplicated().any():
        duplicates = sorted(table.loc[table["accession"].duplicated(), "accession"].unique())
        raise ValueError("Duplicate reference metadata accession(s): " + ", ".join(duplicates))

    lookup = {}
    for row in table.to_dict("records"):
        accession = row["accession"].strip()
        if not accession:
            continue
        lookup[accession] = {
            "origin": normalize_origin(row["origin"]),
            "subtype": normalize_subtype(row["subtype"]),
            "strain": row["strain"].strip() or UNKNOWN,
        }
    return lookup


def apply_reference_metadata(
    parsed: Dict[str, str], lookup: Dict[str, Dict[str, str]]
) -> Dict[str, str]:
    """Fill unavailable legacy-header fields from accession metadata."""
    reference = lookup.get(parsed["accession"])
    if not reference:
        return parsed
    for field in ("origin", "subtype", "strain"):
        if parsed[field] == UNKNOWN and reference[field] != UNKNOWN:
            parsed[field] = reference[field]
    return parsed


def format_reference(hit) -> str:
    """Return the compact ORIGIN|SUBTYPE|STRAIN segment result."""
    return f"{hit['origin']}|{hit['subtype']}|{hit['strain']}"


def join_values(values: Iterable[str]) -> str:
    """Produce a deterministic, spreadsheet-friendly summary value."""
    cleaned = sorted({str(value) for value in values if value and value != UNKNOWN})
    return ";".join(cleaned) if cleaned else UNKNOWN


def build_conclusion(accepted: List[Dict[str, str]], missing: List[str], low: List[str]) -> str:
    """Create an actionable conclusion from all eight segment calls."""
    origins = {hit["origin"] for hit in accepted if hit["origin"] != UNKNOWN}
    subtypes = {hit["subtype"] for hit in accepted if hit["subtype"] != UNKNOWN}
    strains = {hit["strain"] for hit in accepted if hit["strain"] != UNKNOWN}
    metadata_missing = [
        hit["segment"]
        for hit in accepted
        if UNKNOWN in (hit["origin"], hit["subtype"], hit["strain"])
    ]
    non_human = sorted(origin for origin in origins if origin != "HUMAN")

    if missing or low or metadata_missing:
        reasons = []
        if missing:
            reasons.append(f"missing segments: {','.join(missing)}")
        if low:
            reasons.append(f"low-identity segments: {','.join(low)}")
        if metadata_missing:
            reasons.append(f"reference metadata missing for: {','.join(metadata_missing)}")
        if non_human:
            reasons.append(f"observed non-human origin(s): {','.join(non_human)}")
        if len(origins) > 1:
            reasons.append(f"observed mixed origins: {','.join(sorted(origins))}")
        if len(subtypes) > 1:
            reasons.append(f"observed subtype discordance: {','.join(sorted(subtypes))}")
        return "INCONCLUSIVE - " + "; ".join(reasons)

    origin_text = ",".join(sorted(origins))
    subtype_text = ",".join(sorted(subtypes))

    if len(origins) > 1:
        detail = f"ALERT - mixed origins ({origin_text})"
        if len(subtypes) > 1:
            detail += f" and subtypes ({subtype_text})"
        return detail

    origin = next(iter(origins))
    if origin != "HUMAN":
        if len(subtypes) > 1:
            return f"ALERT - non-human origin {origin} with mixed subtypes ({subtype_text})"
        subtype = next(iter(subtypes))
        if len(strains) > 1:
            return f"ALERT - possible {origin} {subtype} reassortment; multiple reference strains"
        return f"FLAG - all segments match one non-human {origin} {subtype} reference strain"

    if len(subtypes) > 1:
        return f"ALERT - human subtype discordance ({subtype_text})"

    subtype = next(iter(subtypes))
    if len(strains) > 1:
        return f"REVIEW - possible within-subtype reassortment; multiple HUMAN {subtype} reference strains"
    return f"CONSISTENT - all segments match one HUMAN {subtype} reference strain"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, help="BLAST outfmt 6 table")
    parser.add_argument("--output", required=True, help="Single-line CSV result")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument(
        "--metadata",
        help="Optional accession/origin/subtype/strain CSV for legacy subject IDs",
    )
    return parser.parse_args()


def read_blast(path: str) -> pd.DataFrame:
    """Read an outfmt 6 file, including a valid empty result."""
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    try:
        return pd.read_csv(path, sep="\t", names=columns)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)


def main() -> None:
    args = parse_args()
    blast = read_blast(args.blast)
    blast["segment"] = blast["qseqid"].apply(infer_segment_from_qseqid)

    metadata_lookup = load_reference_metadata(args.metadata)
    parsed_metadata = (
        blast["sseqid"]
        .apply(lambda subject_id: apply_reference_metadata(parse_subject_id(subject_id), metadata_lookup))
        .apply(pd.Series)
    )
    blast = pd.concat([blast, parsed_metadata], axis=1)

    best = (
        blast.dropna(subset=["segment"])
        .sort_values(["pident", "bitscore", "length"], ascending=False)
        .drop_duplicates("segment")
    )

    row = {"Sample": args.sample}
    accepted = []
    missing = []
    low = []

    for segment in EXPECTED_SEGMENTS:
        candidates = best[best["segment"] == segment]
        if candidates.empty:
            row[segment] = "Missing"
            missing.append(segment)
            continue

        hit = candidates.iloc[0]
        identity = round(float(hit["pident"]), 1)
        reference = format_reference(hit)
        if identity < IDENTITY_THRESHOLD:
            row[segment] = f"TooLow({identity}%):{reference}"
            low.append(segment)
            continue

        row[segment] = f"{reference}({identity}%)"
        accepted.append(
            {
                "segment": segment,
                "origin": hit["origin"],
                "subtype": hit["subtype"],
                "strain": hit["strain"],
            }
        )

    metadata_complete = all(
        UNKNOWN not in (hit["origin"], hit["subtype"], hit["strain"])
        for hit in accepted
    )
    reference_profiles = {
        (hit["origin"], hit["subtype"], hit["strain"])
        for hit in accepted
    }
    if missing or low or not metadata_complete:
        row["Reassortment"] = "Unknown"
    else:
        row["Reassortment"] = "No" if len(reference_profiles) == 1 else "Yes"

    row["Conclusion"] = build_conclusion(accepted, missing, low)
    row["Origins"] = join_values(hit["origin"] for hit in accepted)
    row["Subtypes"] = join_values(hit["subtype"] for hit in accepted)
    row["ReferenceStrains"] = join_values(hit["strain"] for hit in accepted)

    pd.DataFrame([row]).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()

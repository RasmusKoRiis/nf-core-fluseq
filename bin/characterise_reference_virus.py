#!/usr/bin/env python3
"""Classify seasonal influenza samples against reporting reference viruses.

The classifier consumes the existing subclade-nomenclature CSV and the annual
characterisation guideline tables. Only rows explicitly marked
``reporting_category=yes`` can become the primary category. A sample in a
descendant, non-reporting lineage is expressed as the nearest reporting
reference virus plus the intervening guideline mutations and any sample-specific
amino-acid mutations reported by the subclade caller.
"""

import argparse
import csv
import os
import re
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


MISSING = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "unassigned",
    "no ha sequence",
    "no fasta input",
    "unsupported subtype",
}

PROFILE_DEFINITIONS = {
    "H1N1": {
        "filename": "H1_characterisation_guidelines_NH_2025_2026.csv",
        "tokens": ("H1", "H1N1", "H1N1PDM", "H1N1PDM09"),
        "nomenclature_profiles": ("H1N1PDM_HA",),
    },
    "H3N2": {
        "filename": "H3_characterisation_guidelines_NH_2025_2026.csv",
        "tokens": ("H3", "H3N2"),
        "nomenclature_profiles": ("H3N2_HA",),
    },
    "B/Victoria": {
        "filename": "BVIC_characterisation_guidelines_NH_2025_2026.csv",
        "tokens": ("VIC", "VICVIC", "BVIC", "BVICTORIA"),
        "nomenclature_profiles": ("BVIC_HA",),
    },
    "B/Yamagata": {
        "filename": "BYAM_characterisation_guidelines_NH_2025_2026.csv",
        "tokens": ("YAM", "YAMYAM", "BYAM", "BYAMAGATA"),
        "nomenclature_profiles": ("BYAM_HA",),
    },
}

GUIDELINE_COLUMNS = {
    "reference_virus",
    "reference_role",
    "clade",
    "subclade",
    "signature_amino_acid_substitutions",
    "ancestor",
    "reporting_category",
}

OUTPUT_COLUMNS = [
    "Characterisation_Profile",
    "Characterisation_Status",
    "Characterisation_Reference_Virus",
    "Characterisation_Reference_Role",
    "Characterisation_Reporting_Category_Subclade",
    "Characterisation_Closest_Guideline_Reference",
    "Characterisation_Closest_Guideline_Subclade",
    "Characterisation_Guideline_Extra_Mutations",
    "Characterisation_Sample_Extra_Mutations",
    "Characterisation_All_Extra_Mutations",
    "Characterisation_Result",
    "Characterisation_Subclade_Match_Fraction",
    "Characterisation_Missing_Subclade_Mutations",
    "Characterisation_Guideline_Source",
]


def is_missing(value) -> bool:
    return str(value or "").strip().lower() in MISSING


def normalize_token(value) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(value or "").upper())


def label_aliases(value) -> Set[str]:
    """Return equivalent labels, including K/former J.2.4.1-style aliases."""
    raw = str(value or "").strip()
    if is_missing(raw):
        return set()
    aliases = {raw.casefold()}
    leading = re.match(r"^([A-Za-z0-9.]+)", raw)
    if leading:
        aliases.add(leading.group(1).casefold())
    for former in re.findall(r"former\s+([A-Za-z0-9.]+)", raw, flags=re.IGNORECASE):
        aliases.add(former.casefold())
    return aliases


def labels_equal(left, right) -> bool:
    return bool(label_aliases(left).intersection(label_aliases(right)))


def is_descendant(child, parent) -> bool:
    child_aliases = label_aliases(child)
    parent_aliases = label_aliases(parent)
    for child_alias in child_aliases:
        for parent_alias in parent_aliases:
            if child_alias == parent_alias or child_alias.startswith(parent_alias + "."):
                return True
    return False


def label_depth(value) -> int:
    aliases = label_aliases(value)
    return max((alias.count(".") + 1 for alias in aliases), default=0)


def row_label(row: Dict[str, str]) -> str:
    return row.get("subclade", "").strip() or row.get("clade", "").strip()


def valid_reference(row: Optional[Dict[str, str]]) -> bool:
    if not row:
        return False
    reference = row.get("reference_virus", "").strip()
    return not is_missing(reference) and reference.casefold() != "none assigned"


def is_reporting_category(row: Dict[str, str]) -> bool:
    return row.get("reporting_category", "").strip().lower() == "yes" and valid_reference(row)


def choose_profile(subtype: str, row: Dict[str, str]) -> Optional[Tuple[str, Dict[str, object]]]:
    haystacks = {
        normalize_token(subtype),
        normalize_token(row.get("Subclade_Nomenclature_Profile", "")),
    }
    for profile_name, definition in PROFILE_DEFINITIONS.items():
        profile_tokens = {
            normalize_token(token) for token in tuple(definition["tokens"]) + tuple(definition["nomenclature_profiles"])
        }
        if any(any(token and token in haystack for token in profile_tokens) for haystack in haystacks):
            return profile_name, definition
    return None


def read_guidelines(path: str) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing_columns = GUIDELINE_COLUMNS.difference(columns)
        if missing_columns:
            raise ValueError(f"{path} is missing guideline column(s): {', '.join(sorted(missing_columns))}")
        rows = []
        for index, raw_row in enumerate(reader):
            row = {key: str(value or "").strip() for key, value in raw_row.items()}
            row["_index"] = str(index)
            rows.append(row)
        return rows


def match_rank(guideline: Dict[str, str], sample_subclade: str, sample_clade: str) -> Optional[Tuple[int, int, int]]:
    guideline_subclade = guideline.get("subclade", "")
    guideline_clade = guideline.get("clade", "")
    index = int(guideline.get("_index", 0))

    if not is_missing(guideline_subclade) and not is_missing(sample_subclade):
        if labels_equal(sample_subclade, guideline_subclade):
            return 3, label_depth(guideline_subclade), -index
        if is_descendant(sample_subclade, guideline_subclade):
            return 2, label_depth(guideline_subclade), -index

    if not is_missing(guideline_clade) and not is_missing(sample_clade):
        if labels_equal(sample_clade, guideline_clade) or is_descendant(sample_clade, guideline_clade):
            return 1, label_depth(guideline_clade), -index
    return None


def nearest_guideline_row(
    guidelines: Sequence[Dict[str, str]], sample_subclade: str, sample_clade: str
) -> Optional[Dict[str, str]]:
    ranked = [
        (rank, row)
        for row in guidelines
        for rank in [match_rank(row, sample_subclade, sample_clade)]
        if rank is not None
    ]
    return max(ranked, key=lambda item: item[0])[1] if ranked else None


def reporting_base(
    guidelines: Sequence[Dict[str, str]], sample_subclade: str, sample_clade: str
) -> Optional[Dict[str, str]]:
    candidates = []
    for row in guidelines:
        if not is_reporting_category(row):
            continue
        label = row_label(row)
        if not is_missing(sample_subclade) and not is_missing(row.get("subclade", "")):
            if is_descendant(sample_subclade, row["subclade"]):
                candidates.append((label_depth(label), -int(row["_index"]), row))
        elif not is_missing(sample_clade) and not is_missing(row.get("clade", "")):
            if is_descendant(sample_clade, row["clade"]):
                candidates.append((label_depth(label), -int(row["_index"]), row))
    return max(candidates, key=lambda item: item[:2])[2] if candidates else None


def closest_named_reference(
    guidelines: Sequence[Dict[str, str]],
    sample_subclade: str,
    sample_clade: str,
    matched: Optional[Dict[str, str]],
    base: Optional[Dict[str, str]],
) -> Optional[Dict[str, str]]:
    if valid_reference(matched):
        return matched

    candidates = []
    for row in guidelines:
        if not valid_reference(row):
            continue
        rank = match_rank(row, sample_subclade, sample_clade)
        if rank is not None:
            candidates.append((rank, row))
    if candidates:
        return max(candidates, key=lambda item: item[0])[1]
    return base if valid_reference(base) else None


def mutation_key(value: str) -> str:
    match = re.search(r"del\([^)]+\)|[A-Z*]\d+[A-Z*]", value, flags=re.IGNORECASE)
    return (match.group(0) if match else value).casefold()


def unique_values(values: Iterable[str]) -> List[str]:
    output = []
    seen = set()
    for value in values:
        cleaned = str(value or "").strip()
        if is_missing(cleaned):
            continue
        key = mutation_key(cleaned)
        if key not in seen:
            seen.add(key)
            output.append(cleaned)
    return output


def signature_mutations(signature: str) -> List[str]:
    pieces = re.split(r"\s*;\s*|\s+\+\s+", str(signature or "").strip())
    mutations = []
    for piece in pieces:
        cleaned = piece.strip()
        if not cleaned or cleaned.lower().startswith("clade "):
            continue
        cleaned = re.sub(r"^usually\s+", "", cleaned, flags=re.IGNORECASE)
        mutations.append(cleaned)
    return unique_values(mutations)


def guideline_extra_mutations(
    guidelines: Sequence[Dict[str, str]],
    base: Optional[Dict[str, str]],
    sample_subclade: str,
) -> List[str]:
    if not base or is_missing(sample_subclade):
        return []
    base_label = row_label(base)
    base_keys = {
        mutation_key(value) for value in signature_mutations(base.get("signature_amino_acid_substitutions", ""))
    }
    path_rows = []
    for row in guidelines:
        guideline_subclade = row.get("subclade", "")
        if is_missing(guideline_subclade) or labels_equal(guideline_subclade, base_label):
            continue
        if is_descendant(sample_subclade, guideline_subclade) and is_descendant(guideline_subclade, base_label):
            path_rows.append(row)
    path_rows.sort(key=lambda row: (label_depth(row["subclade"]), int(row["_index"])))

    extras = []
    for row in path_rows:
        for mutation in signature_mutations(row.get("signature_amino_acid_substitutions", "")):
            if mutation_key(mutation) not in base_keys:
                extras.append(mutation)
    return unique_values(extras)


def sample_extra_mutations(value: str) -> List[str]:
    mutations = []
    for mutation in str(value or "").split(";"):
        cleaned = mutation.strip()
        if is_missing(cleaned) or cleaned.lower().startswith("nuc:"):
            continue
        mutations.append(cleaned)
    return unique_values(mutations)


def incomplete_subclade_call(row: Dict[str, str]) -> bool:
    missing = row.get("Subclade_Nomenclature_Closest_Subclade_Missing_Mutations", "")
    if not is_missing(missing):
        return True
    fraction = row.get("Subclade_Nomenclature_Subclade_Match_Fraction", "")
    try:
        return float(fraction) < 1.0
    except (TypeError, ValueError):
        return False


def empty_characterisation() -> Dict[str, str]:
    return {column: "NA" for column in OUTPUT_COLUMNS}


def characterise_row(row: Dict[str, str], subtype: str, guidelines_dir: str) -> Dict[str, str]:
    output = empty_characterisation()
    selected = choose_profile(subtype, row)
    if selected is None:
        output["Characterisation_Status"] = "Unsupported subtype"
        output["Characterisation_Result"] = "Unsupported subtype"
        return output

    profile_name, definition = selected
    guideline_path = os.path.join(guidelines_dir, str(definition["filename"]))
    guidelines = read_guidelines(guideline_path)
    sample_subclade = row.get("Subclade_Nomenclature_Subclade", "")
    sample_clade = row.get("Subclade_Nomenclature_Clade", "")
    matched = nearest_guideline_row(guidelines, sample_subclade, sample_clade)
    base = reporting_base(guidelines, sample_subclade, sample_clade)
    closest = closest_named_reference(guidelines, sample_subclade, sample_clade, matched, base)
    guideline_extras = guideline_extra_mutations(guidelines, base, sample_subclade)
    sample_extras = sample_extra_mutations(row.get("Subclade_Nomenclature_Unique_Mutations", ""))
    all_extras = unique_values(guideline_extras + sample_extras)

    output.update(
        {
            "Characterisation_Profile": profile_name,
            "Characterisation_Subclade_Match_Fraction": row.get("Subclade_Nomenclature_Subclade_Match_Fraction", "NA")
            or "NA",
            "Characterisation_Missing_Subclade_Mutations": row.get(
                "Subclade_Nomenclature_Closest_Subclade_Missing_Mutations", "NA"
            )
            or "NA",
            "Characterisation_Guideline_Source": os.path.basename(guideline_path),
            "Characterisation_Guideline_Extra_Mutations": "; ".join(guideline_extras) or "NA",
            "Characterisation_Sample_Extra_Mutations": "; ".join(sample_extras) or "NA",
            "Characterisation_All_Extra_Mutations": "; ".join(all_extras) or "NA",
        }
    )

    if base:
        output["Characterisation_Reference_Virus"] = base["reference_virus"]
        output["Characterisation_Reference_Role"] = base.get("reference_role", "") or "NA"
        output["Characterisation_Reporting_Category_Subclade"] = row_label(base) or "NA"
    if closest:
        output["Characterisation_Closest_Guideline_Reference"] = closest["reference_virus"]
        output["Characterisation_Closest_Guideline_Subclade"] = row_label(closest) or "NA"

    if incomplete_subclade_call(row):
        status = "Review - incomplete subclade call"
    elif base:
        exact = matched is base and labels_equal(sample_subclade or sample_clade, row_label(base))
        status = "Exact reporting category" if exact else "Derived reporting category"
    elif not any(is_reporting_category(item) for item in guidelines):
        status = "No reporting categories in guideline"
    elif matched:
        status = "No reporting category"
    else:
        status = "Unassigned"
    output["Characterisation_Status"] = status

    if base:
        result = f"{base['reference_virus']}-like"
        if all_extras:
            result += " + " + "; ".join(all_extras)
        if status.startswith("Review"):
            result = "Provisional: " + result
    elif closest:
        result = f"No reporting category; closest guideline reference {closest['reference_virus']}-like"
        if sample_extras:
            result += " + " + "; ".join(sample_extras)
    elif status == "No reporting categories in guideline":
        result = "No reporting categories in current guideline"
    else:
        result = status
    output["Characterisation_Result"] = result
    return output


def read_subtype(path: str) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return stripped.split()[0]
    return ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Subclade nomenclature CSV")
    parser.add_argument("--subtype-file", required=True)
    parser.add_argument("--guidelines-dir", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    subtype = read_subtype(args.subtype_file)
    with open(args.input, "r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        input_columns = list(reader.fieldnames or [])
        rows = list(reader)

    output_columns = input_columns + [column for column in OUTPUT_COLUMNS if column not in input_columns]
    with open(args.output, "w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=output_columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            row.update(characterise_row(row, subtype, args.guidelines_dir))
            writer.writerow(row)


if __name__ == "__main__":
    main()

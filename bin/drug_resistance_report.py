#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


RAW_COLUMNS = [
    "M2 inhibtion mutations",
    "NA inhibtion mutations",
    "PA inhibtion mutations",
]
RESULT_COLUMNS = [
    "DR_Res_Adamantine",
    "DR_Res_Oseltamivir",
    "DR_Res_Zanamivir",
    "DR_Res_Peramivir",
    "DR_Res_Laninamivir",
    "DR_Res_Baloxavir",
]
DETAIL_COLUMNS = ["DR_M2_Mut", "DR_NA_Mut", "DR_PA_Mut"]
OUTPUT_COLUMNS = ["SampleID", "OriginalName", "Subtype", *RAW_COLUMNS, *RESULT_COLUMNS, *DETAIL_COLUMNS]
MISSING_VALUES = {"", "NA", "NAN", "NONE"}


def parse_args():
    parser = argparse.ArgumentParser(description="Build the human FASTA drug-resistance-only report.")
    parser.add_argument("--id-map", required=True, type=Path)
    parser.add_argument("--subtype", nargs="*", type=Path, default=[])
    parser.add_argument("--resistance", nargs="*", type=Path, default=[])
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_rows(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def is_missing(value):
    return str(value or "").strip().upper() in MISSING_VALUES


def clean_value(value):
    value = str(value or "").strip()
    return "NA" if is_missing(value) else value


def normalize_resistance_column(column):
    # Some NA translations are named NA1/NA2/etc.; present them in one NA field.
    if re.fullmatch(r"NA\d+ inhibtion mutations", column, flags=re.IGNORECASE):
        return "NA inhibtion mutations"
    return column


def ensure_sample(report, order, sample_id, original_name=None):
    sample_id = str(sample_id or "").strip()
    if not sample_id:
        return None
    if sample_id not in report:
        report[sample_id] = {
            "SampleID": sample_id,
            "OriginalName": str(original_name or sample_id).strip(),
            "Subtype": "NA",
            **{column: "NA" for column in RAW_COLUMNS},
        }
        order.append(sample_id)
    elif original_name and report[sample_id]["OriginalName"] == sample_id:
        report[sample_id]["OriginalName"] = str(original_name).strip()
    return report[sample_id]


def merge_value(target, column, value):
    value = clean_value(value)
    if is_missing(target.get(column)) and not is_missing(value):
        target[column] = value


def classify(value, no_mutation_code):
    value = clean_value(value)
    if value == "NA":
        return "NA"
    if "no matching mutations" in value.lower():
        return no_mutation_code
    return "Review"


def mutation_detail(value, classification):
    value = clean_value(value)
    if value == "NA":
        return "NA"
    return value if classification == "Review" else "No Mutations"


def build_report(id_map, subtype_files, resistance_files):
    report = {}
    order = []

    for row in read_rows(id_map, delimiter="\t"):
        ensure_sample(report, order, row.get("SampleID"), row.get("OriginalName"))

    for path in subtype_files:
        for row in read_rows(path):
            target = ensure_sample(report, order, row.get("Sample") or row.get("SampleID"))
            if target is not None:
                merge_value(target, "Subtype", row.get("Subtype"))

    for path in resistance_files:
        for row in read_rows(path):
            target = ensure_sample(report, order, row.get("Sample") or row.get("SampleID"))
            if target is None:
                continue
            for column, value in row.items():
                normalized_column = normalize_resistance_column(str(column).strip())
                if normalized_column in RAW_COLUMNS:
                    merge_value(target, normalized_column, value)

    rows = []
    for sample_id in order:
        row = report[sample_id]
        for column in RAW_COLUMNS:
            row[column] = clean_value(row.get(column))

        row["DR_Res_Adamantine"] = classify(row["M2 inhibtion mutations"], "AANI")
        for drug in ("Oseltamivir", "Zanamivir", "Peramivir", "Laninamivir"):
            row[f"DR_Res_{drug}"] = classify(row["NA inhibtion mutations"], "AANI")
        row["DR_Res_Baloxavir"] = classify(row["PA inhibtion mutations"], "AANS")

        row["DR_M2_Mut"] = mutation_detail(row["M2 inhibtion mutations"], row["DR_Res_Adamantine"])
        row["DR_NA_Mut"] = mutation_detail(row["NA inhibtion mutations"], row["DR_Res_Oseltamivir"])
        row["DR_PA_Mut"] = mutation_detail(row["PA inhibtion mutations"], row["DR_Res_Baloxavir"])
        rows.append(row)

    return rows


def write_report(rows, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    rows = build_report(args.id_map, args.subtype, args.resistance)
    write_report(rows, args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


SEGMENT_SUFFIXES = {"HA", "NA", "PB1", "PB2", "PA", "NP", "NS", "MP", "M"}


def normalize_sample_name(value: str) -> str:
    sample = str(value).strip()
    if not sample:
        return "NA"
    if "|" in sample:
        sample = sample.split("|", 1)[0].strip()
    if "_" in sample:
        left, right = sample.rsplit("_", 1)
        if left and right.upper() in SEGMENT_SUFFIXES:
            return left.strip() or "NA"
    return sample


def slim_report(input_path: Path, output_path: Path) -> None:
    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    if "Genotype_Genin2" not in frame.columns and "Genotype" in frame.columns:
        frame = frame.rename(columns={"Genotype": "Genotype_Genin2"})
    if "Genotype_Genin2" not in frame.columns:
        frame["Genotype_Genin2"] = "NA"

    sample_column = next(
        (
            candidate
            for candidate in [
                "Sample Name",
                "Sample",
                "SampleID",
                "SequenceID",
                "sample_id",
                "id",
                "ID",
                "Name",
            ]
            if candidate in frame.columns
        ),
        None,
    )
    frame["Sample Name"] = frame[sample_column] if sample_column else input_path.stem
    frame = frame[["Sample Name", "Genotype_Genin2"]]
    frame = frame.applymap(lambda value: str(value).strip() or "NA")
    frame["Sample Name"] = frame["Sample Name"].map(normalize_sample_name)
    frame = frame.drop_duplicates(subset=["Sample Name"], keep="first")

    temporary = output_path.with_suffix(".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the reduced GenIn2 report.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    slim_report(args.input, args.output)


if __name__ == "__main__":
    main()

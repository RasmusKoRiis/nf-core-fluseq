#!/usr/bin/env python3
import argparse
import csv
import os
import re


PROFILES = {
    "H3N2_HA": {
        "subtype_tokens": ("H3", "H3N2"),
        "rule_dir": "H3N2_HA",
        "source": "influenza-clade-nomenclature/seasonal_A-H3N2_HA",
        "features": {"SigPep": 18, "HA1": 66, "HA2": 1053},
    },
    "H1N1pdm_HA": {
        "subtype_tokens": ("H1", "H1N1", "H1N1PDM"),
        "rule_dir": "H1N1pdm_HA",
        "source": "influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA",
        "features": {"SigPep": 21, "HA1": 72, "HA2": 1053},
    },
    "B-Vic_HA": {
        "subtype_tokens": ("VIC", "VICVIC", "B/VIC", "B/VICTORIA"),
        "rule_dir": "B-Vic_HA",
        "source": "influenza-clade-nomenclature/seasonal_B-Vic_HA",
        "features": {"SigPep": 34, "HA1": 79, "HA2": 1120},
    },
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

OUTPUT_COLUMNS = [
    "Sample",
    "Subclade_Nomenclature_Profile",
    "Subclade_Nomenclature_Clade",
    "Subclade_Nomenclature_Subclade",
    "Subclade_Nomenclature_Key_Mutations",
    "Subclade_Nomenclature_Clade_Key_Mutations",
    "Subclade_Nomenclature_Subclade_Match_Fraction",
    "Subclade_Nomenclature_Clade_Match_Fraction",
    "Subclade_Nomenclature_Source",
]


def read_fasta(paths):
    records = []
    for path in paths:
        header = None
        seq = []
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        records.append((path, header, "".join(seq)))
                    header = line[1:].strip()
                    seq = []
                else:
                    seq.append(line)
        if header is not None:
            records.append((path, header, "".join(seq)))
    return records


def read_first_token(path):
    if not path:
        return ""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line:
                    return line.split()[0]
    except FileNotFoundError:
        return ""
    return ""


def normalize_token(value):
    return re.sub(r"[^A-Z0-9/]+", "", str(value).upper())


def choose_profile(subtype, records):
    haystacks = [normalize_token(subtype)]
    for path, header, _seq in records:
        haystacks.append(normalize_token(os.path.basename(path)))
        haystacks.append(normalize_token(header))

    for profile_name, profile in PROFILES.items():
        for haystack in haystacks:
            for token in profile["subtype_tokens"]:
                if normalize_token(token) in haystack:
                    return profile_name, profile
    return "", None


def is_ha_record(path, header):
    text = f"{os.path.basename(path)} {header}".upper()
    if re.search(r"(^|[^A-Z0-9])HA([^A-Z0-9]|$)", text):
        return True
    if re.search(r"(^|[^A-Z0-9])01[-_]?HA([^A-Z0-9]|$)", text):
        return True
    return False


def choose_ha_record(records):
    for record in records:
        if is_ha_record(record[0], record[1]):
            return record
    return records[0] if records else None


def clean_seq(seq):
    return re.sub(r"\s+", "", seq).upper().replace("U", "T")


def translate_codon(codon):
    codon = codon.upper()
    if codon == "---":
        return "-"
    if "-" in codon or len(codon) != 3:
        return "X"
    if re.search(r"[^ACGT]", codon):
        return "X"
    return CODON_TABLE.get(codon, "X")


def observed_state(seq, locus, position, features):
    try:
        pos = int(position)
    except (TypeError, ValueError):
        return ""

    if locus == "nuc":
        idx = pos - 1
        if 0 <= idx < len(seq):
            return seq[idx]
        return ""

    start = features.get(locus)
    if start is None:
        return ""
    idx = start - 1 + (pos - 1) * 3
    return translate_codon(seq[idx:idx + 3])


def format_mutation(rule, observed=None):
    expected = f"{rule['locus']}:{rule['site']}{rule['alt']}"
    if observed is None or observed == rule["alt"]:
        return expected
    return f"{expected}(observed={observed or 'NA'})"


def parse_scalar(value):
    value = value.strip()
    if value in ("[]", "none", "None", "null", "Null", "~"):
        return "" if value != "[]" else []
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    return value


def read_yaml_record(path):
    record = {}
    current_list = None
    current_item = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            indent = len(raw_line) - len(raw_line.lstrip(" "))
            line = raw_line.strip()

            if indent == 0 and not line.startswith("- ") and ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                current_item = None
                if value:
                    record[key] = parse_scalar(value)
                    current_list = None
                else:
                    record[key] = []
                    current_list = key
                continue

            if current_list and line.startswith("- "):
                item_text = line[2:].strip()
                if ":" in item_text:
                    key, value = item_text.split(":", 1)
                    current_item = {key.strip(): parse_scalar(value)}
                    record[current_list].append(current_item)
                else:
                    current_item = parse_scalar(item_text)
                    record[current_list].append(current_item)
                continue

            if current_item is not None and isinstance(current_item, dict) and ":" in line:
                key, value = line.split(":", 1)
                current_item[key.strip()] = parse_scalar(value)

    return record


def read_yaml_rules(path):
    definitions = {}
    parents = {}
    clade_links = {}
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Missing nomenclature YAML directory: {path}")
    for filename in sorted(os.listdir(path)):
        if not filename.endswith((".yml", ".yaml")):
            continue
        record = read_yaml_record(os.path.join(path, filename))
        name = str(record.get("name") or os.path.splitext(filename)[0]).strip()
        if not name:
            continue
        rules = []
        for mutation in record.get("defining_mutations") or []:
            if not isinstance(mutation, dict):
                continue
            locus = str(mutation.get("locus") or "").strip()
            site = str(mutation.get("position") or "").strip()
            alt = str(mutation.get("state") or "").strip()
            if locus and site and alt:
                rules.append({"locus": locus, "site": site, "alt": alt})
        definitions[name] = rules
        parent = str(record.get("parent") or "").strip()
        parents[name] = "" if parent.lower() == "none" else parent
        clade = str(record.get("clade") or "").strip()
        if clade and clade.lower() != "none":
            clade_links[name] = clade
    return definitions, parents, clade_links


def read_clade_yaml_rules(path, subclade_definitions):
    definitions = {}
    parents = {}
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Missing clade YAML directory: {path}")
    for filename in sorted(os.listdir(path)):
        if not filename.endswith((".yml", ".yaml")):
            continue
        record = read_yaml_record(os.path.join(path, filename))
        name = str(record.get("name") or os.path.splitext(filename)[0]).strip()
        if not name:
            continue
        display_name = str(record.get("short_name") or name).strip()
        rules = []
        for mutation in record.get("defining_mutations") or []:
            if not isinstance(mutation, dict):
                continue
            locus = str(mutation.get("locus") or "").strip()
            site = str(mutation.get("position") or "").strip()
            alt = str(mutation.get("state") or "").strip()
            if locus and site and alt:
                rules.append({"locus": locus, "site": site, "alt": alt})
        alias_of = str(record.get("alias_of") or "").strip()
        if not rules and alias_of in subclade_definitions:
            rules = subclade_definitions[alias_of]
        definitions[display_name] = rules
        parent = str(record.get("parent") or "").strip()
        parents[display_name] = "" if parent.lower() == "none" else parent
    return definitions, parents


def evaluate_rule_set(seq, rules, features):
    details = []
    matched = []
    for rule in rules:
        observed = observed_state(seq, rule["locus"], rule["site"], features)
        ok = observed == rule["alt"]
        detail = {
            "rule": rule,
            "observed": observed,
            "matched": ok,
        }
        details.append(detail)
        if ok:
            matched.append(format_mutation(rule, observed))
    total = len(rules)
    nmatch = len(matched)
    fraction = (nmatch / total) if total else 0.0
    return {
        "total": total,
        "matched": nmatch,
        "fraction": fraction,
        "details": details,
        "matched_mutations": matched,
        "exact": total > 0 and nmatch == total,
    }


def parent_depth(name, parents):
    depth = 0
    seen = set()
    cur = name
    while parents.get(cur) and cur not in seen:
        seen.add(cur)
        cur = parents[cur]
        depth += 1
    return depth


def call_hierarchical(seq, definitions, parents, features):
    evaluations = {
        name: evaluate_rule_set(seq, rules, features)
        for name, rules in definitions.items()
    }
    exact = [name for name, ev in evaluations.items() if ev["exact"]]
    if exact:
        name = sorted(
            exact,
            key=lambda n: (parent_depth(n, parents), len(definitions[n]), n),
            reverse=True,
        )[0]
        return name, evaluations[name], True

    if evaluations:
        name = sorted(
            evaluations,
            key=lambda n: (evaluations[n]["fraction"], evaluations[n]["matched"], parent_depth(n, parents), n),
            reverse=True,
        )[0]
        return "Unassigned", evaluations[name], False

    return "Unassigned", {"fraction": 0.0, "matched_mutations": []}, False


def call_flat(seq, definitions, features):
    evaluations = {
        name: evaluate_rule_set(seq, rules, features)
        for name, rules in definitions.items()
    }
    exact = [name for name, ev in evaluations.items() if ev["exact"]]
    if exact:
        name = sorted(exact, key=lambda n: (n.count("."), len(definitions[n]), len(n), n), reverse=True)[0]
        return name, evaluations[name], True
    if evaluations:
        name = sorted(
            evaluations,
            key=lambda n: (evaluations[n]["fraction"], evaluations[n]["matched"], n.count("."), len(n), n),
            reverse=True,
        )[0]
        return "Unassigned", evaluations[name], False
    return "Unassigned", {"fraction": 0.0, "matched_mutations": []}, False


def call_sample(sample_id, subtype, fasta_paths, rules_dir):
    records = read_fasta(fasta_paths)
    profile_name, profile = choose_profile(subtype, records)
    row = {column: "NA" for column in OUTPUT_COLUMNS}
    row["Sample"] = sample_id

    if not records:
        row["Subclade_Nomenclature_Profile"] = "No FASTA input"
        return row
    if profile is None:
        row["Subclade_Nomenclature_Profile"] = "Unsupported subtype"
        return row

    ha_record = choose_ha_record(records)
    if ha_record is None:
        row["Subclade_Nomenclature_Profile"] = profile_name
        row["Subclade_Nomenclature_Clade"] = "No HA sequence"
        row["Subclade_Nomenclature_Subclade"] = "No HA sequence"
        return row

    seq = clean_seq(ha_record[2])
    rule_path = os.path.join(rules_dir, profile["rule_dir"])
    sub_defs, sub_parents, subclade_to_clade = read_yaml_rules(os.path.join(rule_path, "subclades"))
    clade_defs, _clade_parents = read_clade_yaml_rules(os.path.join(rule_path, "clades"), sub_defs)

    subclade, sub_eval, _sub_exact = call_hierarchical(seq, sub_defs, sub_parents, profile["features"])
    clade, clade_eval, _clade_exact = call_flat(seq, clade_defs, profile["features"])
    if subclade != "Unassigned" and subclade in subclade_to_clade:
        linked_clade = subclade_to_clade[subclade]
        if linked_clade in clade_defs:
            clade = linked_clade
            clade_eval = evaluate_rule_set(seq, clade_defs[linked_clade], profile["features"])
    subclade_mutations = ";".join(sub_eval.get("matched_mutations", [])) if subclade != "Unassigned" else ""
    clade_mutations = ";".join(clade_eval.get("matched_mutations", [])) if clade != "Unassigned" else ""

    row.update({
        "Subclade_Nomenclature_Profile": profile_name,
        "Subclade_Nomenclature_Clade": clade,
        "Subclade_Nomenclature_Subclade": subclade,
        "Subclade_Nomenclature_Key_Mutations": subclade_mutations or "NA",
        "Subclade_Nomenclature_Clade_Key_Mutations": clade_mutations or "NA",
        "Subclade_Nomenclature_Subclade_Match_Fraction": f"{sub_eval.get('fraction', 0.0):.3f}",
        "Subclade_Nomenclature_Clade_Match_Fraction": f"{clade_eval.get('fraction', 0.0):.3f}",
        "Subclade_Nomenclature_Source": profile["source"],
    })
    return row


def main():
    parser = argparse.ArgumentParser(description="Call seasonal influenza HA clade/subclade nomenclature rules.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--subtype-file", required=True)
    parser.add_argument("--rules-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("fasta", nargs="+")
    args = parser.parse_args()

    subtype = read_first_token(args.subtype_file)
    row = call_sample(args.sample_id, subtype, args.fasta, args.rules_dir)

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()

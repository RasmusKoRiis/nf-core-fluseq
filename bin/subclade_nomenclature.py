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
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

OUTPUT_COLUMNS = [
    "Sample",
    "Subclade_Nomenclature_Profile",
    "Subclade_Nomenclature_Clade",
    "Subclade_Nomenclature_Clade_Long",
    "Subclade_Nomenclature_Subclade",
    "Subclade_Nomenclature_Lineage_Path",
    "Subclade_Nomenclature_Key_Mutations",
    "Subclade_Nomenclature_Lineage_Additive_Mutations",
    "Subclade_Nomenclature_Lineage_Key_Mutations",
    "Subclade_Nomenclature_Clade_Key_Mutations",
    "Subclade_Nomenclature_Closest_Subclade",
    "Subclade_Nomenclature_Closest_Subclade_Missing_Mutations",
    "Subclade_Nomenclature_Unique_Mutations",
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


def read_single_fasta(path):
    records = read_fasta([path])
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return clean_seq(records[0][2]).replace("-", "")


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


def build_reference_to_query_map(reference, query):
    reference = clean_seq(reference).replace("-", "")
    query = clean_seq(query).replace("-", "")
    n = len(reference)
    m = len(query)
    gap = -5
    match = 2
    mismatch = -1

    prev = [j * gap for j in range(m + 1)]
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        trace[0][j] = 2
    for i in range(1, n + 1):
        curr = [i * gap] + [0] * m
        trace[i][0] = 1
        ref_base = reference[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ref_base == query[j - 1] else mismatch)
            up = prev[j] + gap
            left = curr[j - 1] + gap
            if diag >= up and diag >= left:
                curr[j] = diag
                trace[i][j] = 0
            elif up >= left:
                curr[j] = up
                trace[i][j] = 1
            else:
                curr[j] = left
                trace[i][j] = 2
        prev = curr

    ref_to_query = [None] * (n + 1)
    i = n
    j = m
    while i > 0 or j > 0:
        step = trace[i][j]
        if i > 0 and j > 0 and step == 0:
            ref_to_query[i] = j
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or step == 1):
            ref_to_query[i] = None
            i -= 1
        else:
            j -= 1
    return ref_to_query, query


def translate_codon(codon):
    codon = codon.upper()
    if codon == "---":
        return "-"
    if "-" in codon or len(codon) != 3:
        return "X"
    if re.search(r"[^ACGT]", codon):
        return "X"
    return CODON_TABLE.get(codon, "X")


def observed_state(seq, locus, position, features, ref_to_query=None):
    try:
        pos = int(position)
    except (TypeError, ValueError):
        return ""

    if locus == "nuc":
        query_pos = ref_to_query[pos] if ref_to_query and pos < len(ref_to_query) else pos
        if query_pos is None:
            return "-"
        idx = query_pos - 1
        if 0 <= idx < len(seq):
            return seq[idx]
        return ""

    start = features.get(locus)
    if start is None:
        return ""
    ref_positions = [start + (pos - 1) * 3 + offset for offset in range(3)]
    codon = []
    for ref_pos in ref_positions:
        query_pos = ref_to_query[ref_pos] if ref_to_query and ref_pos < len(ref_to_query) else ref_pos
        if query_pos is None:
            return "-"
        idx = query_pos - 1
        if not (0 <= idx < len(seq)):
            return ""
        codon.append(seq[idx])
    return translate_codon("".join(codon))


def feature_lengths(reference, features):
    reference = clean_seq(reference).replace("-", "")
    ordered = sorted(features.items(), key=lambda item: item[1])
    lengths = {}
    for index, (locus, start) in enumerate(ordered):
        end = ordered[index + 1][1] - 1 if index + 1 < len(ordered) else len(reference)
        lengths[locus] = max(0, (end - start + 1) // 3)
    return lengths


def observed_mutations(seq, reference, features, ref_to_query=None):
    reference = clean_seq(reference).replace("-", "")
    mutations = []
    seen = set()
    aligned_positions = []
    if ref_to_query:
        aligned_positions = [index for index in range(1, len(ref_to_query)) if ref_to_query[index] is not None]
    first_aligned = min(aligned_positions) if aligned_positions else 1
    last_aligned = max(aligned_positions) if aligned_positions else len(reference)

    for pos, ref_base in enumerate(reference, start=1):
        if pos < first_aligned or pos > last_aligned:
            continue
        observed = observed_state(seq, "nuc", pos, features, ref_to_query)
        if not observed or observed in {"N"} or observed == ref_base:
            continue
        key = ("nuc", str(pos), observed)
        if key not in seen:
            seen.add(key)
            mutations.append({"locus": "nuc", "site": str(pos), "alt": observed})

    for locus, length in feature_lengths(reference, features).items():
        for pos in range(1, length + 1):
            start = features[locus] + (pos - 1) * 3
            end = start + 2
            if start < first_aligned or end > last_aligned:
                continue
            ref_aa = observed_state(reference, locus, pos, features)
            observed = observed_state(seq, locus, pos, features, ref_to_query)
            if not observed or observed in {"X"} or observed == ref_aa:
                continue
            key = (locus, str(pos), observed)
            if key not in seen:
                seen.add(key)
                mutations.append({"locus": locus, "site": str(pos), "alt": observed})

    return mutations


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


def mutation_key(rule):
    return (rule["locus"], str(rule["site"]), rule["alt"])


def mutation_site_key(rule):
    return (rule["locus"], str(rule["site"]))


def merge_rules(rule_sets):
    merged = []
    indexes = {}
    for rules in rule_sets:
        for rule in rules:
            key = mutation_site_key(rule)
            if key in indexes:
                merged.pop(indexes[key])
                merged.append(rule)
                indexes = {mutation_site_key(item): idx for idx, item in enumerate(merged)}
            else:
                indexes[key] = len(merged)
                merged.append(rule)
    return merged


def lineage_names(name, parents):
    names = []
    seen = set()
    current = name
    while current and current not in seen:
        seen.add(current)
        names.append(current)
        current = parents.get(current, "")
    return list(reversed(names))


def lineage_rule_set(name, definitions, parents):
    return merge_rules(definitions.get(item, []) for item in lineage_names(name, parents))


def format_rule_list(rules):
    return ";".join(format_mutation(rule) for rule in rules) or "NA"


def lineage_additive_mutations(name, definitions, parents):
    names = lineage_names(name, parents)
    final_rules = {}
    owners = {}
    order = []
    for lineage_name in names:
        for rule in definitions.get(lineage_name, []):
            key = mutation_site_key(rule)
            if key in final_rules:
                order.remove(key)
            order.append(key)
            final_rules[key] = rule
            owners[key] = lineage_name

    grouped = {lineage_name: [] for lineage_name in names}
    for key in order:
        grouped[owners[key]].append(final_rules[key])

    return (
        " | ".join(
            f"{lineage_name}:{format_rule_list(grouped[lineage_name])}"
            for lineage_name in names
            if grouped[lineage_name]
        )
        or "NA"
    )


def read_clade_yaml_rules(path, subclade_definitions, subclade_parents):
    definitions = {}
    parents = {}
    long_names = {}
    raw_rules = {}
    raw_parents = {}
    raw_display_names = {}
    raw_aliases = {}
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
        parent = str(record.get("parent") or "").strip()
        raw_rules[name] = rules
        raw_parents[name] = "" if parent.lower() == "none" else parent
        raw_display_names[name] = display_name
        raw_aliases[name] = alias_of

    memo = {}

    def clade_lineage_rules(name):
        if name in memo:
            return memo[name]
        alias_of = raw_aliases.get(name, "")
        if alias_of in subclade_definitions:
            rules = merge_rules(
                [
                    lineage_rule_set(alias_of, subclade_definitions, subclade_parents),
                    raw_rules.get(name, []),
                ]
            )
        else:
            rule_sets = []
            parent = raw_parents.get(name, "")
            if parent in raw_rules:
                rule_sets.append(clade_lineage_rules(parent))
            rule_sets.append(raw_rules.get(name, []))
            rules = merge_rules(rule_sets)
        memo[name] = rules
        return rules

    for name in sorted(raw_rules):
        display_name = raw_display_names[name]
        rules = clade_lineage_rules(name)
        definitions[display_name] = rules
        long_names[display_name] = name
        parents[display_name] = raw_parents.get(name, "")
        if name != display_name:
            definitions[name] = rules
            long_names[name] = name
            parents[name] = raw_parents.get(name, "")
    return definitions, parents, long_names


def evaluate_rule_set(seq, rules, features, ref_to_query=None):
    details = []
    matched = []
    for rule in rules:
        observed = observed_state(seq, rule["locus"], rule["site"], features, ref_to_query)
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


def missing_mutations(evaluation):
    missing = []
    for detail in evaluation.get("details", []):
        if not detail.get("matched"):
            missing.append(format_mutation(detail["rule"], detail.get("observed")))
    return missing


def parent_depth(name, parents):
    depth = 0
    seen = set()
    cur = name
    while parents.get(cur) and cur not in seen:
        seen.add(cur)
        cur = parents[cur]
        depth += 1
    return depth


def call_hierarchical(seq, definitions, parents, features, ref_to_query=None):
    evaluations = {name: evaluate_rule_set(seq, rules, features, ref_to_query) for name, rules in definitions.items()}
    lineage_evaluations = {
        name: evaluate_rule_set(seq, lineage_rule_set(name, definitions, parents), features, ref_to_query)
        for name in definitions
    }

    if lineage_evaluations:
        name = sorted(
            lineage_evaluations,
            key=lambda n: (
                lineage_evaluations[n]["matched"],
                lineage_evaluations[n]["fraction"],
                parent_depth(n, parents),
                len(lineage_rule_set(n, definitions, parents)),
                n,
            ),
            reverse=True,
        )[0]
        evaluations[name]["candidate"] = name
        evaluations[name]["lineage_evaluation"] = lineage_evaluations[name]
        return (
            (name if lineage_evaluations[name]["exact"] else "Unassigned"),
            evaluations[name],
            lineage_evaluations[name]["exact"],
        )

    return "Unassigned", {"fraction": 0.0, "matched_mutations": [], "candidate": ""}, False


def linked_clade_for_subclade(subclade, parents, clade_links):
    seen = set()
    current = subclade
    while current and current not in seen:
        seen.add(current)
        clade = clade_links.get(current)
        if clade:
            return clade
        current = parents.get(current, "")
    return ""


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

    seq = clean_seq(ha_record[2]).replace("-", "")
    rule_path = os.path.join(rules_dir, profile["rule_dir"])
    reference = read_single_fasta(os.path.join(rule_path, "reference.fasta"))
    ref_to_query, seq = build_reference_to_query_map(reference, seq)
    sub_defs, sub_parents, subclade_to_clade = read_yaml_rules(os.path.join(rule_path, "subclades"))
    clade_defs, _clade_parents, clade_long_names = read_clade_yaml_rules(
        os.path.join(rule_path, "clades"),
        sub_defs,
        sub_parents,
    )

    subclade, sub_eval, sub_exact = call_hierarchical(seq, sub_defs, sub_parents, profile["features"], ref_to_query)
    candidate_subclade = sub_eval.get("candidate", "")
    lineage_target = subclade if subclade != "Unassigned" else candidate_subclade
    lineage_path = " -> ".join(lineage_names(lineage_target, sub_parents)) if lineage_target else ""
    lineage_additions = lineage_additive_mutations(lineage_target, sub_defs, sub_parents) if lineage_target else ""
    lineage_eval = {"matched_mutations": [], "fraction": 0.0}
    lineage_rules = []
    if lineage_target:
        lineage_rules = lineage_rule_set(lineage_target, sub_defs, sub_parents)
        lineage_eval = evaluate_rule_set(seq, lineage_rules, profile["features"], ref_to_query)
    clade = "Unassigned"
    clade_long = "Unassigned"
    clade_eval = {"fraction": 0.0, "matched_mutations": []}
    if lineage_target:
        linked_clade = linked_clade_for_subclade(lineage_target, sub_parents, subclade_to_clade)
        if linked_clade:
            clade = linked_clade
            clade_long = clade_long_names.get(linked_clade, linked_clade)
            if linked_clade in clade_defs:
                clade_eval = evaluate_rule_set(seq, clade_defs[linked_clade], profile["features"], ref_to_query)
    subclade_mutations = ";".join(sub_eval.get("matched_mutations", [])) if lineage_target else ""
    lineage_mutations = ";".join(lineage_eval.get("matched_mutations", [])) if lineage_target else ""
    clade_mutations = ";".join(clade_eval.get("matched_mutations", [])) if clade != "Unassigned" else ""
    sub_eval_for_report = sub_eval.get("lineage_evaluation", lineage_eval) if lineage_target else sub_eval
    missing = ";".join(missing_mutations(sub_eval_for_report)) if lineage_target and not sub_exact else ""
    lineage_keys = {mutation_key(rule) for rule in lineage_rules}
    unique_mutations = ";".join(
        format_mutation(rule)
        for rule in observed_mutations(seq, reference, profile["features"], ref_to_query)
        if mutation_key(rule) not in lineage_keys
    )

    row.update(
        {
            "Subclade_Nomenclature_Profile": profile_name,
            "Subclade_Nomenclature_Clade": clade,
            "Subclade_Nomenclature_Clade_Long": clade_long,
            "Subclade_Nomenclature_Subclade": lineage_target or "Unassigned",
            "Subclade_Nomenclature_Lineage_Path": lineage_path or "NA",
            "Subclade_Nomenclature_Key_Mutations": subclade_mutations or "NA",
            "Subclade_Nomenclature_Lineage_Additive_Mutations": lineage_additions or "NA",
            "Subclade_Nomenclature_Lineage_Key_Mutations": lineage_mutations or "NA",
            "Subclade_Nomenclature_Clade_Key_Mutations": clade_mutations or "NA",
            "Subclade_Nomenclature_Closest_Subclade": candidate_subclade or "NA",
            "Subclade_Nomenclature_Closest_Subclade_Missing_Mutations": missing or "NA",
            "Subclade_Nomenclature_Unique_Mutations": unique_mutations or "NA",
            "Subclade_Nomenclature_Subclade_Match_Fraction": f"{sub_eval_for_report.get('fraction', 0.0):.3f}",
            "Subclade_Nomenclature_Clade_Match_Fraction": f"{clade_eval.get('fraction', 0.0):.3f}",
            "Subclade_Nomenclature_Source": profile["source"],
        }
    )
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

"""Check that local module documentation exists and describes real ports."""

import json
import re
from pathlib import Path

import jsonschema
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
MODULES = sorted((ROOT / "modules/local").glob("*/main.nf"))
SCHEMA = json.loads((ROOT / "tests/schemas/module-meta.schema.json").read_text())


def declarations(source, section):
    """Read the single-line channel declarations used by these local modules."""
    source = re.sub(r"/\*.*?\*/", "", source, flags=re.S)
    source = re.sub(r"(?m)^\s*//.*$", "", source)
    start = re.search(rf"(?m)^\s*{section}:\s*$", source)
    if start is None:
        return []
    block = re.split(
        r"(?m)^\s*(?:input|output|when|script|stub|shell|exec):",
        source[start.end() :],
        maxsplit=1,
    )[0]
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    assert all(
        re.match(r"(?:tuple|val|path|file)\b", line) for line in lines
    ), "Update the metadata contract check to support new declaration syntax"
    return lines


@pytest.mark.parametrize("main", MODULES, ids=lambda p: p.parent.name)
def test_local_module_metadata(main):
    metadata_path = main.with_name("meta.yml")
    assert metadata_path.is_file(), f"Missing module metadata: {metadata_path}"
    metadata = yaml.safe_load(metadata_path.read_text())
    jsonschema.Draft202012Validator(SCHEMA).validate(metadata)
    source = main.read_text()
    processes = re.findall(r"(?m)^process\s+(\w+)\s*\{", source)
    assert processes == [metadata["name"].upper()]

    inputs = declarations(source, "input")
    assert len(metadata["input"]) == len(inputs)
    for declaration, documented in zip(inputs, metadata["input"]):
        is_tuple = declaration.startswith("tuple ")
        assert isinstance(documented, list) == is_tuple
        elements = documented if is_tuple else [documented]
        names = re.findall(r"\b(?:val|path|file)\s*\(?\s*(\w+)", declaration)
        assert [next(iter(item)) for item in elements] == names

    outputs = declarations(source, "output")
    emits = [re.search(r"emit:\s*(\w+)", line)[1] for line in outputs]
    assert list(metadata["output"]) == emits
    for declaration, emit in zip(outputs, emits):
        documented = metadata["output"][emit]
        is_tuple = declaration.startswith("tuple ")
        assert isinstance(documented[0], list) == is_tuple
        elements = documented[0] if is_tuple else documented
        assert len(elements) == len(re.findall(r"\b(?:val|path|file)\b", declaration))
        patterns = re.findall(r"""\bpath\s*\(?\s*["']([^"']+)["']""", declaration)
        assert [
            details["pattern"] for item in elements for details in item.values() if "pattern" in details
        ] == patterns

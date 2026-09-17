"""Protect runtime database updates and workflow failure handling."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_all_processes_default_to_ignored_errors():
    base_config = (ROOT / "conf" / "base.config").read_text(encoding="utf-8")
    assert "errorStrategy = 'ignore'" in base_config

    for module in (ROOT / "modules").glob("**/main.nf"):
        source = module.read_text(encoding="utf-8")
        strategies = re.findall(r"errorStrategy\s*(?:=\s*)?['\"]([^'\"]+)['\"]", source)
        assert all(strategy == "ignore" for strategy in strategies), module


def test_flumut_attempts_update_with_bundled_database_fallback():
    source = (ROOT / "modules" / "local" / "flumut" / "main.nf").read_text(encoding="utf-8")
    update = source.index("if ! flumut --update")
    analysis = source.index("flumut -m")

    assert update < analysis
    assert "using the database bundled with the container" in source


def test_subclade_rules_always_download_latest_main_branches():
    source = (ROOT / "modules" / "local" / "subclade_nomenclature_rules" / "main.nf").read_text(
        encoding="utf-8"
    )

    assert "cache false" in source
    assert source.count("archive/refs/heads/main.tar.gz") == 3
    assert not re.search(r"archive/[0-9a-f]{40}\.tar\.gz", source)


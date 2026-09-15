"""Check that server resource caps remain process settings rather than pipeline parameters."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_server_resource_limits_are_not_schema_validated_parameters():
    source = (ROOT / "conf" / "server.config").read_text(encoding="utf-8")
    params_block = source.split("params {", 1)[1].split("}", 1)[0]

    assert "max_cpus" not in params_block
    assert "max_memory" not in params_block
    assert "max_time" not in params_block
    assert "process.resourceLimits" in source
    assert "cpus: 16" in source
    assert "memory: '256.GB'" in source
    assert "time: '20.h'" in source

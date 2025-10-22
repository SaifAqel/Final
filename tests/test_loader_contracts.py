from io_loader import load_config
from units import Q_

def test_required_stage_fields_present():
    stages, gas, water, _ = load_config("config/stages.yaml", "config/streams.yaml")
    for s in stages:
        assert "inner_diameter" in s.spec
        assert "inner_length" in s.spec
        assert s.spec["inner_diameter"].check("[length]")
        assert s.spec["inner_length"].check("[length]")

def test_stream_types_and_units():
    _, gas, water, _ = load_config("config/stages.yaml", "config/streams.yaml")
    assert gas.mass_flow.check("[mass]/[time]")
    assert gas.T.check("[temperature]")
    assert gas.P.check("[pressure]")
    assert water.mass_flow.check("[mass]/[time]")
    assert water.h.check("[energy]/[mass]")
    assert water.P.check("[pressure]")

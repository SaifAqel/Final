from io_loader import load_config, _q
from units import Q_
import yaml, pathlib

def test__q_parses_quantity():
    q = _q({"value": 5, "unit": "m"})
    assert isinstance(q, Q_)
    assert q.to("m").magnitude == 5

def test_load_config_minimal_paths_ok():
    stages, gas, water, drum = load_config("config/stages.yaml", "config/streams.yaml", "config/drum.yaml")
    assert len(stages) >= 1
    assert gas.mass_flow.check("[mass] / [time]")
    assert water.h.check("[energy] / [mass]")
    assert drum.Di.check("[length]")

import pytest
import logging
from logging_utils import setup_logging
from io_loader import load_config
from pipeline import SixStageCounterflow

@pytest.mark.usefixtures("fake_gas")
def test_pipeline_runs_small_and_outputs_order():
    setup_logging(logging.ERROR)
    stages, gas, water, _ = load_config("config/stages.yaml", "config/streams.yaml")
    # shrink length to speed up and reduce external dependence
    for s in stages:
        if "inner_length" in s.spec:
            s.spec["inner_length"] *= 0.1
    gh, wh = SixStageCounterflow(stages).run(gas, water)
    assert len(gh) > 0 and len(wh) > 0
    # outlet indexing contract: hot is last, cold is first (counterflow)
    assert gh[-1].T.to("K").magnitude < gh[0].T.to("K").magnitude
    assert wh[0].h.to("J/kg").magnitude > wh[-1].h.to("J/kg").magnitude

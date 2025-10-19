# ======================== tests/test_smoke.py ========================
# Optional smoke test. Runs without external deps.
from logging_utils import setup_logging
from io_loader import load_config
from pipeline import SixStageCounterflow
import logging

def test_pipeline_smoke():
    setup_logging(logging.ERROR)
    stages, gas, water = load_config()
    pipe = SixStageCounterflow(stages)
    g_hist, w_hist = pipe.run(gas, water)
    assert len(g_hist) > 0
    assert len(w_hist) > 0

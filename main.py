import argparse, logging
from logging_utils import setup_logging
from io_loader import load_config
from pipeline import SixStageCounterflow

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="config/stages.yaml")
    ap.add_argument("--streams", default="config/streams.yaml")
    ap.add_argument("--log", default="DEBUG")
    args = ap.parse_args()

    setup_logging(args.log)

    stages, gas, water = load_config(args.stages, args.streams)
    pipe = SixStageCounterflow(stages)
    gh, wh = pipe.run(gas, water)

    g_out, w_out = gh[-1], wh[-1]
    logging.getLogger("main").info("done", extra={"stage": stages[-1].name, "step": "outlet"})
    print(f"Stages: {len(stages)}")
    print(f"Gas T_out: {g_out.T:.3f~P}")
    print(f"Gas P_out: {g_out.P:.1f~P}")
    print(f"Water h_out: {w_out.h:.2f~P}")

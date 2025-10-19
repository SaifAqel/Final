import argparse, logging
from logging_utils import setup_logging
from io_loader import load_config
from pipeline import SixStageCounterflow
from units import Q_

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="config/stages.yaml")
    ap.add_argument("--streams", default="config/streams.yaml")
    ap.add_argument("--log", default="TRACE")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()

    setup_logging(args.log)

    stages, gas, water = load_config(args.stages, args.streams)
    if args.fast:
        for s in stages:
            s.L = Q_(0.2, "m")
            if "UA_per_m" in s.spec: s.spec["UA_per_m"] *= 0.8

    pipe = SixStageCounterflow(stages)
    gh, wh = pipe.run(gas, water)
    g_out, w_out = gh[-1], wh[-1]
    print(f"OK | stages={len(stages)} | Gas Tout={g_out.T:.2f~P} | Water h_out={w_out.h:.2f~P}")

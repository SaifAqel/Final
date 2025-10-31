# main.py
import argparse, logging
from logging_utils import setup_logging
from io_loader import load_config
from pipeline import SixStageCounterflow
from geometry import GeometryBuilder

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="config/stages.yaml")
    ap.add_argument("--streams", default="config/streams.yaml")
    ap.add_argument("--drum", default="config/drum.yaml")
    ap.add_argument("--log", default="DEBUG")
    args = ap.parse_args()

    setup_logging(args.log)

    geom, gas, water, drum = load_config(args.stages, args.streams, args.drum)
    stages = GeometryBuilder(drum).enrich(geom)

    pipe = SixStageCounterflow(stages)
    stage_results = pipe.run(gas, water)

    last = stage_results[-1]
    g_out = last.steps[-1].gas
    # cold outlet is the side opposite hot outlet; that corresponds to the final water
    w_out = last.steps[0].water

    logging.getLogger("main").info("done", extra={"stage": stages[-1].name, "step": "outlet"})
    print(f"Stages: {len(stages)}")
    print(f"Gas T_out: {g_out.T:.3f~P}")
    print(f"Gas P_out: {g_out.P:.1f~P}")
    print(f"Water h_out: {w_out.h:.2f~P}")
    print(f"Stage Q: {last.Q_stage:.3f~P}, UA: {last.UA_stage:.3f~P}")
    print(stage_results)

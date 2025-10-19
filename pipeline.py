# pipeline.py
from dataclasses import replace
from typing import List, Tuple
import logging
from logging_utils import trace_calls
from solver import StageSolver
from models import GasStream, WaterStream, HXStage
log = logging.getLogger("pipeline")

class SixStageCounterflow:
    def __init__(self, stages: List[HXStage]):
        self.stages = stages

    @trace_calls()
    def run(self, gas: GasStream, water: WaterStream) -> Tuple[List[GasStream], List[WaterStream]]:
        gas_hist: List[GasStream] = []     # was: gas_hist, water_hist = []
        water_hist: List[WaterStream] = []

        g_cur = gas
        w_cur = water
        for stg in self.stages:
            log.info("stage-start", extra={"stage": stg.name, "step": "start"})
            g_list, w_list = StageSolver(stg, g_cur, w_cur).solve()

            # append without duplicating the shared boundary between stages
            if gas_hist:   gas_hist.extend(g_list[1:])
            else:          gas_hist.extend(g_list)

            if water_hist: water_hist.extend(w_list[1:])
            else:          water_hist.extend(w_list)

            # next-stage inlets: hot takes last of g_list; cold takes index 0 of w_list (counterflow)
            g_cur = replace(g_list[-1])
            w_cur = replace(w_list[0])
            log.info("stage-end", extra={"stage": stg.name, "step": "end"})
        return gas_hist, water_hist

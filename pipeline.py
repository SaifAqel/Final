from typing import List, Tuple
from dataclasses import replace
import logging
from logging_utils import trace_calls
from models import GasStream, WaterStream, HXStage
from solver import StageSolver
log = logging.getLogger("pipeline")

class SixStageCounterflow:
    def __init__(self, stages: List[HXStage]):
        self.stages = stages

    @trace_calls()
    def run(self, gas: GasStream, water: WaterStream) -> Tuple[List[GasStream], List[WaterStream]]:
        gas_hist: List[GasStream] = []
        water_hist: List[WaterStream] = []
        for stg in self.stages:
            log.info("stage-start", extra={"stage": stg.name, "step": "start"})
            solver = StageSolver(stg, gas, water)
            g_list, w_list = solver.solve()
            gas_hist.extend(g_list)
            water_hist.extend(w_list)
            gas = replace(g_list[-1])
            water = replace(w_list[-1])
            log.info("stage-end", extra={"stage": stg.name, "step": "end"})
        return gas_hist, water_hist

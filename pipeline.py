# pipeline.py
from dataclasses import replace
from typing import List
import logging
from logging_utils import trace_calls
from solver import StageSolver
from models import GasStream, WaterStream, HXStage
from results import StageResult

log = logging.getLogger("pipeline")

class SixStageCounterflow:
    def __init__(self, stages: List[HXStage]):
        self.stages = stages

    def run(self, gas: GasStream, water: WaterStream) -> List[StageResult]:
        results: List[StageResult] = []

        g_cur = gas
        w_cur = water
        for stg in self.stages:
            log.info("stage-start", extra={"stage": stg.name, "step": "start"})
            res: StageResult = StageSolver(stg, g_cur, w_cur).solve()
            results.append(res)

            # next-stage inlets:
            # hot takes last gas state at x = L
            g_cur = replace(res.steps[-1].gas)
            # cold takes first recorded water state; we keep counterflow across stages
            w_cur = replace(res.steps[0].water)

            log.info("stage-end", extra={"stage": stg.name, "step": "end"})

        return results

from dataclasses import replace
from typing import List, Tuple
import logging
from logging_utils import trace_calls
from units import Q_, ureg
from models import HXStage, GasStream, WaterStream
from physics import heat_rate_per_length, rhs

log = logging.getLogger("solver")

class StageSolver:
    def __init__(self, stage: HXStage, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = replace(gas, stage=stage.name)
        self.water = replace(water, stage=stage.name)

    @trace_calls()
    def solve(self, dx: Q_ = Q_(0.1,"m"), tol_dT: Q_ = Q_(2.0,"K")) -> Tuple[List[GasStream], List[WaterStream]]:
        L: Q_ = self.stage.L
        x: Q_ = Q_(0.0,"m")
        gas_hist: List[GasStream] = []
        water_hist: List[WaterStream] = []
        step_i = 0
        while x < L:
            step_i += 1
            hr = heat_rate_per_length(self.gas, self.water, self.stage)
            derivs = rhs(self.gas, self.water, hr["qprime"], self.stage)
            dT_est: Q_ = abs((derivs["dTgdx"] * dx).to("K"))
            if dT_est > tol_dT:
                dx *= 0.5
                log.trace("cut step", extra={"stage": self.stage.name, "step": step_i})
                continue
            if dT_est < 0.25 * tol_dT:
                dx *= 1.2

            g_next = replace(
                self.gas,
                T=(self.gas.T + derivs["dTgdx"] * dx).to("K"),
                P=(self.gas.P + derivs["dpgdx"] * dx).to("Pa"),
            )
            w_next = replace(
                self.water,
                h=(self.water.h + derivs["dhwdx"] * dx).to("J/kg"),
            )
            gas_hist.append(self.gas)
            water_hist.append(self.water)
            self.gas, self.water = g_next, w_next
            x += dx
            if step_i % 10 == 0:
                log.info("advancing", extra={"stage": self.stage.name, "step": step_i})
        return gas_hist, water_hist

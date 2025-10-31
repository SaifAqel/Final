# solver.py
from dataclasses import replace
from typing import List, Optional, Tuple
from results import StepResult, StageResult
from units import Q_
from models import HXStage, GasStream, WaterStream
from physics import wall_resistance, fouling_resistances
from water_htc import water_htc
from props import WaterProps
from gas_htc import gas_htc, cp_gas
import numpy as np
import logging

log = logging.getLogger("solver")

class StageSolver:
    def __init__(self, stage: HXStage, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = replace(gas, stage=stage.name)
        self.water = replace(water, stage=stage.name)

    @staticmethod
    def solve_step(g, w, stage, Tgw_guess, Tww_guess, qprime_guess, i: Optional[int]=None, x: Optional[Q_]=None, dx: Optional[Q_]=None) -> StepResult:
        spec = stage.spec
        Pg = spec["hot_wet_P"]
        Pw = spec["cold_wet_P"]
        Tg = g.T
        Tw = WaterProps.T_from_Ph(w.P, w.h)
        Tgw = Tgw_guess
        Tww = Tww_guess
        qprime = qprime_guess
        alpha = 0.25
        tolT = Q_(1e-3,"K"); tolq = Q_(1e-3,"W/m"); maxit = 10

        for _ in range(maxit):
            h_g = gas_htc(g, spec, Tgw)
            qpp_cold = (qprime / Pw).to("W/m^2")
            h_c, boiling = water_htc(w, stage, Tww, qpp_cold)
            Rfg, Rfc = fouling_resistances(spec)
            Rw = wall_resistance(spec)
            Rg = (1/(h_g*Pg)).to("K*m/W")
            Rc = (1/(h_c*Pw)).to("K*m/W")

            UA_prime = (1/(Rg + Rfg + Rw + Rfc + Rc)).to("W/K/m")

            qprime_new = (UA_prime * (Tg - Tw)).to("W/m")

            qpp_hot  = (qprime_new / Pg).to("W/m^2")
            qpp_cold = (qprime_new / Pw).to("W/m^2")

            Tgw_new = (Tg - qpp_hot/h_g - qpp_hot*Rfg*Pg).to("K")
            Tww_new = (Tw + qpp_cold*Rw*Pw + qpp_cold*Rfc*Pw + qpp_cold/h_c).to("K")

            dTgw = abs(Tgw_new - Tgw); dTww = abs(Tww_new - Tww); dq = abs(qprime_new - qprime)
            if dTgw < tolT and dTww < tolT and dq < tolq:
                Tgw, Tww, qprime = Tgw_new, Tww_new, qprime_new
                break

            Tgw = (alpha*Tgw_new + (1-alpha)*Tgw).to("K")
            Tww = (alpha*Tww_new + (1-alpha)*Tww).to("K")
            qprime = (alpha*qprime_new + (1-alpha)*qprime).to("W/m")

        return StepResult(
            i=i, x=x, dx=dx,
            gas=g, water=w,
            Tgw=Tgw, Tww=Tww,
            UA_prime=UA_prime,
            qprime=qprime,
            boiling=boiling
        )

    def _n_segments(self) -> int:
        # simple rule: ~0.25 m segments, min 10, max 400
        L = self.stage.spec["inner_length"].to("m").magnitude
        N = int(np.clip(np.ceil(L / 0.25), 10, 400))
        return max(1, N)

    def solve(self) -> StageResult:
        stg = self.stage
        L = stg.spec["inner_length"].to("m")
        N = self._n_segments()
        dx = (L / N).to("m")

        steps: List[StepResult] = []

        # marching variables
        g = self.gas            # hot stream starts at x=0
        w = self.water          # cold stream starts at x=L (counterflow)

        # initial guesses
        Tgw_guess = g.T
        Tww_guess = WaterProps.T_from_Ph(w.P, w.h)
        qprime_guess = Q_(0.0, "W/m")

        Q_stage = Q_(0.0, "W")
        UA_stage = Q_(0.0, "W/K")

        # march left->right for gas, right->left for water
        for i in range(N):
            xL = (i * dx).to("m")

            # solve local film temps and q' with current bulk states
            step = self.solve_step(
                g, w, stg,
                Tgw_guess=Tgw_guess, Tww_guess=Tww_guess, qprime_guess=qprime_guess,
                i=i, x=xL, dx=dx
            )

            # heat over segment
            dq = (step.qprime * dx).to("W")      # positive if gas hotter than water

            # update bulk states to segment outlet positions
            # gas: T decreases
            cpg = cp_gas(g)                      # J/kg/K
            dTg = (dq / (g.mass_flow * cpg)).to("K")
            g = replace(g, T=(g.T - dTg).to("K"))

            # water: h increases
            dhw = (dq / w.mass_flow).to("J/kg")
            w = replace(w, h=(w.h + dhw).to("J/kg"))

            # accumulate integrals
            Q_stage = (Q_stage + dq).to("W")
            UA_stage = (UA_stage + (step.UA_prime * dx)).to("W/K")

            # store step at the segment exit location x+dx for gas
            steps.append(replace(step, gas=g, water=w))

            # propagate guesses
            Tgw_guess = step.Tgw
            Tww_guess = step.Tww
            qprime_guess = step.qprime

        return StageResult(
            stage_name=stg.name,
            stage_kind=stg.kind,
            steps=steps,
            Q_stage=Q_stage,
            UA_stage=UA_stage
        )

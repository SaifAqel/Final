# solver.py
from dataclasses import replace
from typing import List, Tuple
import logging
from logging_utils import trace_calls
from units import Q_, ureg
from models import HXStage, GasStream, WaterStream
from physics import cp_gas, ua_per_m
from props import WaterProps

log = logging.getLogger("solver")

class StageSolver:
    def __init__(self, stage: HXStage, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = replace(gas, stage=stage.name)     # hot enters at x=0
        self.water = replace(water, stage=stage.name) # cold enters at x=L

    @trace_calls()
    def solve(self, N: int = 60) -> Tuple[List[GasStream], List[WaterStream]]:
        L: Q_ = self.stage.L
        dx: Q_ = (L / N).to("m")

        # Initialize local marching states
        g = self.gas                      # at x=0  (hot inlet)
        w = self.water                    # at x=L  (cold inlet)

        # Recover initial water temperature safely once
        try:
            Tw = WaterProps.T_from_Ph(w.P, w.h)
        except Exception:
            # fall back near saturated liquid temp at P
            Tw = (WaterProps.Tsat(w.P) - Q_(5.0, "K")).to("K")

        gas_hist: List[GasStream] = [g]
        # For counterflow, we’ll build water from the hot-inlet end by prepending
        water_hist: List[WaterStream] = [w]

        UA_per_m: Q_ = ua_per_m(self.stage)

        for i in range(N):
            # Heat rate per length using current interfacial ΔT
            qprime = UA_per_m * (g.T - Tw)              # W/m

            # Hot side derivative
            cpg = cp_gas(g)
            dTgdx = - qprime / (g.mass_flow * cpg)      # K/m
            dpgdx = Q_(0.0, "Pa/m")

            # Cold side derivative in +x coordinates:
            # dh/dx = + q'/m ; but cold marches from x=L to 0, so use -dx step
            cpw = WaterProps.cp_from_PT(w.P, Tw)
            dTwdx = - qprime / (w.mass_flow * cpw)      # K/m

            # Advance hot forward (+dx)
            g_next = replace(
                g,
                T=(g.T + dTgdx * dx).to("K"),
                P=(g.P + dpgdx * dx).to("Pa"),
            )

            # Advance cold backward (−dx): temperature decreases by dTwdx*dx in x sense,
            # but we are stepping toward smaller x, so subtract
            Tw_next = (Tw - dTwdx * dx).to("K")
            # keep Tw in a sane band to avoid IAPWS blowups
            if Tw_next.magnitude < 273.16: Tw_next = Q_(273.16, "K")
            if Tw_next.magnitude > 1073.0: Tw_next = Q_(1073.0, "K")

            # Water enthalpy from energy balance (exact, no inversion):
            h_next = (w.h + (qprime / w.mass_flow) * dx).to("J/kg")
            w_next = replace(w, h=h_next)

            gas_hist.append(g_next)
            # Build water history from hot-inlet end: prepend so index 0 is at x=0
            water_hist.insert(0, w_next)

            g, w, Tw = g_next, w_next, Tw_next

            if (i+1) % 10 == 0:
                log.info("advancing", extra={"stage": self.stage.name, "step": i+1})

        # include final states
        return gas_hist, water_hist

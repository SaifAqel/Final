from dataclasses import replace
from typing import List, Tuple
import logging
from logging_utils import _fmt
from units import Q_
from models import HXStage, GasStream, WaterStream
from physics import cp_gas, ua_per_m
from props import WaterProps

log = logging.getLogger("solver")

class StageSolver:
    def __init__(self, stage: HXStage, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = replace(gas, stage=stage.name)
        self.water = replace(water, stage=stage.name)

    def solve(self, N: int = 60) -> Tuple[List[GasStream], List[WaterStream]]:
        if "inner_length" not in self.stage.spec:
            raise KeyError(f"{self.stage.name}: 'inner_length' required in spec")
        L: Q_ = self.stage.spec["inner_length"].to("m")
        dx: Q_ = (L / N).to("m")

        g = self.gas
        w = self.water

        try:
            Tw = WaterProps.T_from_Ph(w.P, w.h)
        except Exception as e:
            raise RuntimeError(f"{self.stage.name}: invalid water state for T_from_Ph") from e

        gas_hist: List[GasStream] = [g]
        water_hist: List[WaterStream] = [w]
        qpp_hist = []

        for i in range(N):
            if qpp_hist:
                self.stage.spec["qpp_guess"] = qpp_hist[-1]
            UA_per_m_val, qpp = ua_per_m(g, w, self.stage)
            qpp_hist.append(qpp)
            self.stage.spec["qpp_hist"] = qpp_hist
            UA_per_m_val, qpp = ua_per_m(g, w, self.stage)
            qprime = UA_per_m_val * (g.T - Tw)

            cpg = cp_gas(g)
            dTgdx = - qprime / (g.mass_flow * cpg)
            dpgdx = Q_(0.0, "Pa/m")  # unchanged model assumption

            cpw = WaterProps.cp_from_PT(w.P, Tw)
            dTwdx = - qprime / (w.mass_flow * cpw)

            g_next = replace(g, T=(g.T + dTgdx * dx).to("K"), P=(g.P + dpgdx * dx).to("Pa"))
            Tw_next = (Tw - dTwdx * dx).to("K")

            h_next = (w.h + (qprime / w.mass_flow) * dx).to("J/kg")
            w_next = replace(w, h=h_next)

            gas_hist.append(g_next)
            water_hist.insert(0, w_next)

            g, w, Tw = g_next, w_next, Tw_next

            if (i + 1) % 10 == 0:
                log.info(
                    "g: m=%s T=%s P=%s | w: m=%s h=%s T=%s P=%s | UA'=%s q'=%s x=%s/%s"
                    % (_fmt(g.mass_flow), _fmt(g.T), _fmt(g.P),
                       _fmt(w.mass_flow), _fmt(w.h), _fmt(Tw), _fmt(w.P),
                       _fmt(UA_per_m_val), _fmt(qprime),
                       _fmt((i + 1) * dx), _fmt(L)),
                    extra={"stage": self.stage.name, "step": i + 1},
                )

        return gas_hist, water_hist

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
    @staticmethod
    def solve_wall_state(g, w, stage, Tg, Tw, qpp):
        A = stage.spec["inner_perimeter"]
        n = 100
        tq = Q_(1e2, "watt / meter**2")   # heat-flux atol
        tT = Q_(1e-3, "kelvin") 
        T_wg = Tg
        T_ww = Tw
        for _ in range(n):
            UA_per_m, qpp, _ = ua_per_m(g, w, stage, T_wg, T_ww, qpp)
            qpp_new = (UA_per_m * (Tg - Tw) / A)
            if abs(qpp_new - qpp) < tq and abs(T_wg - Tg) < tT and abs(T_ww - Tw) < tT:
                return T_wg, T_ww, UA_per_m, qpp_new
            qpp = qpp_new
        return T_wg, T_ww, UA_per_m, qpp


    def solve(self, N: int = 60) -> Tuple[List[GasStream], List[WaterStream]]:
        if "inner_length" not in self.stage.spec:
            raise KeyError(f"{self.stage.name}: 'inner_length' required in spec")
        L: Q_ = self.stage.spec["inner_length"].to("m")
        dx: Q_ = (L / N).to("m")

        g = self.gas
        w = self.water

        Tw = WaterProps.T_from_Ph(w.P, w.h)

        gas_hist: List[GasStream] = [g]
        water_hist: List[WaterStream] = [w]
        qpp_hist = []

        for i in range(N):

            if qpp_hist:
                qpp_guess = qpp_hist[-1]
            else:
                # one dry call to get UA0 using trivial wall temps and dummy qpp
                UA0, _, boiling = ua_per_m(g, w, self.stage, g.T, Tw, Q_(0.0, "W/m^2"))
                A = self.stage.spec["inner_perimeter"]
                qpp_guess = (UA0 * (g.T - Tw) / A).to("W/m^2")

            T_wg, T_ww, UA_per_m_val, qpp = self.solve_wall_state(g, w, self.stage, g.T, Tw, qpp_guess)
            qpp_hist.append(qpp)
            self.stage.spec["qpp_hist"] = qpp_hist
            self.stage.spec["qpp_guess"] = qpp

            A = self.stage.spec["inner_diameter"]
            qprime = UA_per_m_val * (g.T - Tw)

            cpg = cp_gas(g)
            dTgdx = - qprime / (g.mass_flow * cpg)
            dpgdx = Q_(0.0, "Pa/m")

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
                    "g: m=%s T=%s P=%s | w: m=%s h=%s T=%s P=%s | T_wg=%s T_ww=%s UA'=%s q''=%s x=%s/%s boiling=%s"
                    % (_fmt(g.mass_flow), _fmt(g.T), _fmt(g.P),
                    _fmt(w.mass_flow), _fmt(w.h), _fmt(Tw), _fmt(w.P),
                    _fmt(T_wg), _fmt(T_ww), _fmt(UA_per_m_val), _fmt(qpp),
                    _fmt((i + 1) * dx), _fmt(L), _fmt(boiling)),
                    extra={"stage": self.stage.name, "step": i + 1},
                )

        return gas_hist, water_hist


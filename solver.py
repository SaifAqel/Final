from dataclasses import replace
from typing import List, Tuple
import logging
from logging_utils import _fmt
from units import Q_
from models import HXStage, GasStream, WaterStream
from physics import cp_gas, gas_htc, wall_resistance, fouling_resistances
from water_htc import water_htc
from props import WaterProps

log = logging.getLogger("solver")

class StageSolver:
    def __init__(self, stage: HXStage, gas: GasStream, water: WaterStream):
        self.stage = stage
        self.gas = replace(gas, stage=stage.name)
        self.water = replace(water, stage=stage.name)

    @staticmethod
    def solve_step(g, w, stage, Tgw_guess, Tww_guess, qprime_guess):
        spec = stage.spec
        Pg = spec["inner_perimeter"]
        Pw = spec["cold_wet_P"]
        Tg = g.T
        Tw = WaterProps.T_from_Ph(w.P, w.h)
        Tgw = Tgw_guess
        Tww = Tww_guess
        qprime = qprime_guess
        alpha = 0.5
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
                Tgw, Tww, qprime = Tgw_new, Tww_new, qprime_new   # use actual converged values
                break

            Tgw = (alpha*Tgw_new + (1-alpha)*Tgw).to("K")
            Tww = (alpha*Tww_new + (1-alpha)*Tww).to("K")
            qprime = (alpha*qprime_new + (1-alpha)*qprime).to("W/m")
        
        h_g = gas_htc(g, spec, Tgw)
        qpp_cold = (qprime / Pw).to("W/m^2")
        h_c, boiling = water_htc(w, stage, Tww, qpp_cold)
        Rfg, Rfc = fouling_resistances(spec)
        Rw = wall_resistance(spec)
        Rg = (1/(h_g*Pg)).to("K*m/W")
        Rc = (1/(h_c*Pw)).to("K*m/W")
        UA_prime = (1/(Rg + Rfg + Rw + Rfc + Rc)).to("W/K/m")
        qprime = (UA_prime * (Tg - Tw)).to("W/m")

        return Tgw, Tww, UA_prime, qprime, boiling

    def solve(self, N: int = 60) -> Tuple[List[GasStream], List[WaterStream]]:
        L: Q_ = self.stage.spec["inner_length"].to("m")
        dx: Q_ = (L / N).to("m")

        g = self.gas
        w = self.water

        Tw = WaterProps.T_from_Ph(w.P, w.h)

        gas_hist: List[GasStream] = [g]
        water_hist: List[WaterStream] = [w]
        qprime_hist = []
        Tgw_hist = []
        Tww_hist = []

        for i in range(N):
            if qprime_hist:
                qprime_guess = qprime_hist[-1]
            if Tgw_hist:
                Tgw_guess = Tgw_hist[-1]
            if Tww_hist:
                Tww_guess = Tww_hist[-1]
            else:
                Tgw_guess = g.T
                Tww_guess = Tw
                qprime_guess = Q_(1e-9, "W/m")
                
            Tgw, Tww, UA_prime, qprime, boiling = self.solve_step(g, w, self.stage, Tgw_guess, Tww_guess, qprime_guess)
            qprime_hist.append(qprime)
            Tgw_hist.append(Tgw)
            Tww_hist.append(Tww)

            cpg = cp_gas(g)
            dTgdx = - qprime / (g.mass_flow * cpg)
            dpgdx = Q_(0.0, "Pa/m")
            dhwdx = (qprime / w.mass_flow).to("J/kg/m")

            T_next = (g.T + dTgdx * dx).to("K")
            P_next = (g.P + dpgdx * dx).to("Pa")
            h_next = (w.h + dhwdx * dx).to("J/kg")

            g_next = replace(g, T=T_next, P=P_next)
            w_next = replace(w, h=h_next)

            gas_hist.append(g_next)
            water_hist.insert(0, w_next)

            g, w = g_next, w_next

            if (i + 1) % 10 == 0:
                log.info(
                    "g: m=%s T=%s P=%s | w: m=%s h=%s T=%s P=%s | T_wg=%s T_ww=%s UA'=%s q''=%s x=%s/%s boiling=%s"
                    % (_fmt(g.mass_flow), _fmt(g.T), _fmt(g.P),
                    _fmt(w.mass_flow), _fmt(w.h), _fmt(Tw), _fmt(w.P),
                    _fmt(Tgw), _fmt(Tww), _fmt(UA_prime), _fmt(qprime),
                    _fmt((i + 1) * dx), _fmt(L), _fmt(boiling)),
                    extra={"stage": self.stage.name, "step": i + 1},
                )

        return gas_hist, water_hist


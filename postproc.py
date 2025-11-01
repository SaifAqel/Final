# postproc.py
from __future__ import annotations
from typing import List
import pandas as pd
from results import GlobalProfile
from props import WaterProps, GasProps

_gas = GasProps()  # reuse for all rows

def _mag_or_nan(q, unit):
    return q.to(unit).magnitude if q is not None else float("nan")

def profile_to_dataframe(gp: "GlobalProfile") -> "pd.DataFrame":
    rows = []
    for i in range(len(gp.x)):
        g = gp.gas[i]
        w = gp.water[i]

        # water side primitives
        Tw = WaterProps.T_from_Ph(w.P, w.h)
        xq = WaterProps.quality_from_Ph(w.P, w.h)
        Two_phase = xq is not None

        # water temperatures
        if Two_phase:
            Tw = WaterProps.Tsat(w.P)
        else:
            Tw = WaterProps.T_from_Ph(w.P, w.h)

        # water props
        if Two_phase:
            w_cp  = None
            w_mu  = None
            w_k   = None
            w_rho = WaterProps.rho_from_Px(w.P, xq) if xq is not None else None
        else:
            w_cp  = WaterProps.cp_from_Ph(w.P, w.h)
            w_mu  = WaterProps.mu_from_Ph(w.P, w.h)
            w_k   = WaterProps.k_from_Ph(w.P, w.h)
            w_rho = WaterProps.rho_from_Ph(w.P, w.h)

        # gas props (need composition from StepResult; already present in gp.gas[i])
        g_h   = _gas.h_sensible(g.T, g.P, g.comp)
        g_cp  = _gas.cp(g.T, g.P, g.comp)
        g_mu  = _gas.mu(g.T, g.P, g.comp)
        g_k   = _gas.k(g.T, g.P, g.comp)
        g_rho = _gas.rho(g.T, g.P, g.comp)

        row = {
            "stage_index": gp.stage_index[i],
            "stage_name": gp.stage_name[i],
            "i": i,
            "x[m]": gp.x[i].to("m").magnitude,
            "dx[m]": gp.dx[i].to("m").magnitude,
            "qprime[W/m]": gp.qprime[i].to("W/m").magnitude,
            "UA_prime[W/K/m]": gp.UA_prime[i].to("W/K/m").magnitude,

            # gas stream state + props
            "gas_T[K]": g.T.to("K").magnitude,
            "gas_P[Pa]": g.P.to("Pa").magnitude,
            "gas_h[J/kg]": g_h.to("J/kg").magnitude,
            "gas_cp[J/kg/K]": g_cp.to("J/kg/K").magnitude,
            "gas_mu[Pa*s]": g_mu.to("Pa*s").magnitude,
            "gas_k[W/m/K]": g_k.to("W/m/K").magnitude,
            "gas_rho[kg/m^3]": g_rho.to("kg/m^3").magnitude,

            # water stream state + props
            "water_h[J/kg]": w.h.to("J/kg").magnitude,
            "water_P[Pa]": w.P.to("Pa").magnitude,
            "water_T[K]": Tw.to("K").magnitude,
            "water_cp[J/kg/K]": _mag_or_nan(w_cp, "J/kg/K"),
            "water_mu[Pa*s]": _mag_or_nan(w_mu, "Pa*s"),
            "water_k[W/m/K]": _mag_or_nan(w_k, "W/m/K"),
            "water_rho[kg/m^3]": _mag_or_nan(w_rho, "kg/m^3"),
            "water_x[-]": _mag_or_nan(xq, ""),
            "h_gas[W/m^2/K]": gp.h_g[i].to("W/m^2/K").magnitude,   # <-- add
            "h_water[W/m^2/K]": gp.h_c[i].to("W/m^2/K").magnitude, # <-- add
            "boiling": bool(xq is not None),
        }

        # optional: write gas composition so CSV alone is self-sufficient
        for sp, q in (g.comp or {}).items():
            row[f"y_{sp}[-]"] = q.to("").magnitude

        rows.append(row)

    return pd.DataFrame(rows)


def summary_from_profile(gp: "GlobalProfile") -> tuple[list[dict], float, float]:
    rows = []
    Q_total = 0.0
    UA_total = 0.0
    # group consecutive rows by stage_index
    import itertools
    for k, grp in itertools.groupby(range(len(gp.x)), key=lambda i: gp.stage_index[i]):
        idxs = list(grp)
        name = gp.stage_name[idxs[0]]
        # integrals
        Q_stage = sum((gp.qprime[i] * gp.dx[i]).to("W").magnitude for i in idxs)
        UA_stage = sum((gp.UA_prime[i] * gp.dx[i]).to("W/K").magnitude for i in idxs)
        # endpoints along gas x in this stage
        gas_in_T = gp.gas[idxs[0]].T.to("K").magnitude
        gas_out_T = gp.gas[idxs[-1]].T.to("K").magnitude
        water_in_h = gp.water[idxs[-1]].h.to("J/kg").magnitude   # counter-current at x=L
        water_out_h = gp.water[idxs[0]].h.to("J/kg").magnitude   # counter-current at x=0
        row = {
            "stage_index": k,
            "stage_name": name,
            "stage_kind": "",   # fill if you extend GlobalProfile to include kind per point or a map
            "Q_stage[W]": Q_stage,
            "UA_stage[W/K]": UA_stage,
            "gas_in_T[K]": gas_in_T,
            "gas_out_T[K]": gas_out_T,
            "water_in_h[J/kg]": water_in_h,
            "water_out_h[J/kg]": water_out_h,
        }
        rows.append(row)
        Q_total += Q_stage
        UA_total += UA_stage
    return rows, Q_total, UA_total

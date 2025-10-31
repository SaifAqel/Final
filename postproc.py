# =========================================================
# FILE: postproc.py  (optional)
# =========================================================
from __future__ import annotations
from typing import List
import pandas as pd

from results import GlobalProfile
from props import WaterProps

def profile_to_dataframe(gp: "GlobalProfile") -> "pd.DataFrame":
    rows = []
    for i in range(len(gp.x)):
        w = gp.water[i]
        Tw = WaterProps.T_from_Ph(w.P, w.h)
        xq = WaterProps.quality_from_Ph(w.P, w.h)
        rows.append({
            "stage_index": gp.stage_index[i],
            "stage_name": gp.stage_name[i],
            "i": i,
            "x[m]": gp.x[i].to("m").magnitude,
            "dx[m]": gp.dx[i].to("m").magnitude,
            "qprime[W/m]": gp.qprime[i].to("W/m").magnitude,
            "UA_prime[W/K/m]": gp.UA_prime[i].to("W/K/m").magnitude,
            # wall temperatures are not in GlobalProfile; omit Tgw/Tww or add to profile if needed
            "boiling": bool(WaterProps.quality_from_Ph(w.P, w.h) is not None),
            "gas_T[K]": gp.gas[i].T.to("K").magnitude,
            "gas_P[Pa]": gp.gas[i].P.to("Pa").magnitude,
            "water_h[J/kg]": w.h.to("J/kg").magnitude,
            "water_P[Pa]": w.P.to("Pa").magnitude,
            "water_T[K]": Tw.to("K").magnitude,
            "water_x[-]": (xq.to("").magnitude if xq is not None else float("nan")),
        })
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

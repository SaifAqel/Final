# =========================================================
# FILE: postproc.py  (optional)
# =========================================================
from __future__ import annotations
from typing import List
import pandas as pd

from results import StageResult

def results_to_dataframe(stage_results: List[StageResult]) -> "pd.DataFrame":
    """
    Build a tidy DataFrame of per-step results across all stages.
    Columns: stage_index, stage_name, i, x[m], dx[m], qprime[W/m], UA_prime[W/K/m],
             Tgw[K], Tww[K], boiling[bool],
             gas_T[K], gas_P[Pa], water_h[J/kg], water_P[Pa], water_T[K], water_x[- or NaN]
    """
    rows = []
    from props import WaterProps
    for st in stage_results:
        for s in st.steps:
            Tw = WaterProps.T_from_Ph(s.water.P, s.water.h)
            xq = WaterProps.quality_from_Ph(s.water.P, s.water.h)
            rows.append({
                "stage_index": s.stage_index,
                "stage_name": s.stage_name,
                "i": s.i,
                "x[m]": s.x.to("m").magnitude if s.x is not None else None,
                "dx[m]": s.dx.to("m").magnitude if s.dx is not None else None,
                "qprime[W/m]": s.qprime.to("W/m").magnitude,
                "UA_prime[W/K/m]": s.UA_prime.to("W/K/m").magnitude,
                "Tgw[K]": s.Tgw.to("K").magnitude,
                "Tww[K]": s.Tww.to("K").magnitude,
                "boiling": bool(s.boiling),
                "gas_T[K]": s.gas.T.to("K").magnitude,
                "gas_P[Pa]": s.gas.P.to("Pa").magnitude,
                "water_h[J/kg]": s.water.h.to("J/kg").magnitude,
                "water_P[Pa]": s.water.P.to("Pa").magnitude,
                "water_T[K]": Tw.to("K").magnitude,
                "water_x[-]": (xq.to("").magnitude if xq is not None else float("nan")),
            })
    return pd.DataFrame(rows)

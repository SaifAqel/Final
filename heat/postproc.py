# postproc.py
from __future__ import annotations
import pandas as pd
from common.results import GlobalProfile
from common.props import WaterProps, GasProps
from common.units import Q_
from heat.gas_htc import emissivity 

_gas = GasProps()  # reuse for all rows

def _mag_or_nan(q, unit):
    return q.to(unit).magnitude if q is not None else float("nan")

def profile_to_dataframe(gp: "GlobalProfile", *, remap_water: bool = True) -> "pd.DataFrame":
    # Build per-stage [first,last] index ranges once
    stage_ranges: dict[int, tuple[int, int]] = {}
    for i in range(len(gp.x)):
        k = gp.stage_index[i]
        if k not in stage_ranges:
            stage_ranges[k] = [i, i]
        else:
            stage_ranges[k][1] = i
    stage_ranges = {k: (v[0], v[1]) for k, v in stage_ranges.items()}

    rows = []
    for i in range(len(gp.x)):
        g = gp.gas[i]

        # stage index and geometry (same for gas + water in this row)
        k_stage = gp.stage_index[i]
        sr_stage = gp.stage_results[k_stage]
        A_hot = sr_stage.hot_flow_A   # gas-side flow area
        A_cold = sr_stage.cold_flow_A # water-side flow area
        Dh_hot = sr_stage.hot_Dh      # gas-side hydraulic diameter
        Dh_cold = sr_stage.cold_Dh    # water-side hydraulic diameter

        # mirrored index for water within the same stage block
        if remap_water:
            i0, iN = stage_ranges[gp.stage_index[i]]
            j = i0 + (iN - i)          # mirror: i -> j
        else:
            j = i

        w = gp.water[j]

        # water side primitives (from mirrored j)
        xq = WaterProps.quality_from_Ph(w.P, w.h)
        Two_phase = xq is not None

        if Two_phase:
            Tw = WaterProps.Tsat(w.P)
            w_cp = w_mu = w_k = None
            w_rho = WaterProps.rho_from_Px(w.P, xq) if xq is not None else None
        else:
            Tw   = WaterProps.T_from_Ph(w.P, w.h)
            w_cp = WaterProps.cp_from_Ph(w.P, w.h)
            w_mu = WaterProps.mu_from_Ph(w.P, w.h)
            w_k  = WaterProps.k_from_Ph(w.P, w.h)
            w_rho = WaterProps.rho_from_Ph(w.P, w.h)

        # gas props (local i)
        g_h   = _gas.h_sensible(g.T, g.P, g.comp)
        g_cp  = _gas.cp(g.T, g.P, g.comp)
        g_mu  = _gas.mu(g.T, g.P, g.comp)
        g_k   = _gas.k(g.T, g.P, g.comp)
        g_rho = _gas.rho(g.T, g.P, g.comp)

        # NEW: gas velocity and Reynolds number
        gas_V = (g.mass_flow / (g_rho * A_hot)).to("m/s")
        Re_gas = (g_rho * gas_V * Dh_hot / g_mu).to("").magnitude

        # NEW: water velocity and Reynolds number
        if w_rho is not None and A_cold is not None:
            water_V = (w.mass_flow / (w_rho * A_cold)).to("m/s")
        else:
            water_V = None

        if w_rho is not None and w_mu is not None and water_V is not None:
            Re_water = (w_rho * water_V * Dh_cold / w_mu).to("").magnitude
        else:
            Re_water = float("nan")

        # NEW: gas emissivity (H2O + CO2) using same model as gas_htc
        y = g.comp or {}
        yH2O = y.get("H2O", Q_(0.0, "")).to("").magnitude
        yCO2 = y.get("CO2", Q_(0.0, "")).to("").magnitude
        P_Pa = g.P.to("Pa").magnitude
        pH2O = yH2O * P_Pa
        pCO2 = yCO2 * P_Pa
        # mean beam length same logic as _mean_beam_length (0.9 * Dh)
        Lb_m = (0.9 * Dh_hot).to("m").magnitude
        gas_eps = emissivity(
            g.T.to("K").magnitude,
            pH2O,
            pCO2,
            Lb_m,
        )

        row = {
            "stage_index": gp.stage_index[i],
            "stage_name": gp.stage_name[i],
            "i": i,
            "x[m]": gp.x[i].to("m").magnitude,
            "dx[m]": gp.dx[i].to("m").magnitude,
            "qprime[W/m]": gp.qprime[i].to("W/m").magnitude,
            "UA_prime[W/K/m]": gp.UA_prime[i].to("W/K/m").magnitude,

            # gas stream state + props (local i)
            "gas_T[K]": g.T.to("K").magnitude,
            "gas_P[Pa]": g.P.to("Pa").magnitude,
            "gas_h[J/kg]": g_h.to("J/kg").magnitude,
            "gas_cp[J/kg/K]": g_cp.to("J/kg/K").magnitude,
            "gas_mu[Pa*s]": g_mu.to("Pa*s").magnitude,
            "gas_k[W/m/K]": g_k.to("W/m/K").magnitude,
            "gas_rho[kg/m^3]": g_rho.to("kg/m^3").magnitude,

            # NEW: gas velocity & Reynolds number & emissivity
            "gas_V[m/s]": gas_V.to("m/s").magnitude,
            "Re_gas[-]": Re_gas,
            "gas_eps[-]": gas_eps,

            # water stream state + props (from mirrored j)
            "water_h[J/kg]": w.h.to("J/kg").magnitude,
            "water_P[Pa]": w.P.to("Pa").magnitude,
            "water_T[K]": Tw.to("K").magnitude,
            "water_cp[J/kg/K]": _mag_or_nan(w_cp, "J/kg/K"),
            "water_mu[Pa*s]": _mag_or_nan(w_mu, "Pa*s"),
            "water_k[W/m/K]": _mag_or_nan(w_k, "W/m/K"),
            "water_rho[kg/m^3]": _mag_or_nan(w_rho, "kg/m^3"),
            "water_x[-]": _mag_or_nan(xq, ""),

            # NEW: water velocity and Reynolds (NaN where undefined, e.g. 2-phase)
            "water_V[m/s]": _mag_or_nan(water_V, "m/s") if isinstance(water_V, Q_) else float("nan"),
            "Re_water[-]": Re_water,

            # local HTCs and ΔP stay aligned with gas x
            "h_gas[W/m^2/K]": gp.h_g[i].to("W/m^2/K").magnitude,
            "h_water[W/m^2/K]": gp.h_c[i].to("W/m^2/K").magnitude,

            "dP_fric[Pa]": gp.dP_fric[i].to("Pa").magnitude,
            "dP_minor[Pa]": gp.dP_minor[i].to("Pa").magnitude,
            "dP_total[Pa]": gp.dP_total[i].to("Pa").magnitude,

            "boiling": bool(xq is not None),
        }

        # optional gas composition
        for sp, q in (g.comp or {}).items():
            row[f"y_{sp}[-]"] = q.to("").magnitude

        rows.append(row)


    return pd.DataFrame(rows)



def summary_from_profile(gp: "GlobalProfile") -> tuple[list[dict], float, float]:
    rows = []
    Q_total = 0.0
    UA_total = 0.0

    import itertools
    for k, grp in itertools.groupby(range(len(gp.x)), key=lambda i: gp.stage_index[i]):
        idxs = list(grp)
        name = gp.stage_name[idxs[0]]

        # integrals
        Q_stage = sum((gp.qprime[i] * gp.dx[i]).to("W").magnitude for i in idxs)
        UA_stage = sum((gp.UA_prime[i] * gp.dx[i]).to("W/K").magnitude for i in idxs)

        # ΔP sums
        dP_fric = sum(gp.dP_fric[i].to("Pa").magnitude for i in idxs)
        dP_minor = sum(gp.dP_minor[i].to("Pa").magnitude for i in idxs)
        dP_total = sum(gp.dP_total[i].to("Pa").magnitude for i in idxs)

        # endpoints along gas x in this stage
        gas_in_T = gp.gas[idxs[0]].T.to("K").magnitude
        gas_out_T = gp.gas[idxs[-1]].T.to("K").magnitude
        water_in_h = gp.water[idxs[-1]].h.to("J/kg").magnitude   # counter-current at x=L
        water_out_h = gp.water[idxs[0]].h.to("J/kg").magnitude   # counter-current at x=0

        row = {
            "stage_index": k,
            "stage_name": name,
            "stage_kind": "",
            "Q_stage[W]": Q_stage,
            "UA_stage[W/K]": UA_stage,
            "gas_in_T[K]": gas_in_T,
            "gas_out_T[K]": gas_out_T,
            "water_in_h[J/kg]": water_in_h,
            "water_out_h[J/kg]": water_out_h,
            # new ΔP stage totals
            "ΔP_stage_fric[Pa]": dP_fric,
            "ΔP_stage_minor[Pa]": dP_minor,
            "ΔP_stage_total[Pa]": dP_total,
        }
        rows.append(row)
        Q_total += Q_stage
        UA_total += UA_stage

    return rows, Q_total, UA_total

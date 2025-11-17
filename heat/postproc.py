# postproc.py
from __future__ import annotations
import pandas as pd
from common.results import GlobalProfile, CombustionResult
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

    stage_offsets: dict[int, Q_] = {}
    offset = Q_(0.0, "m")
    for k, sr in enumerate(gp.stage_results):
        stage_offsets[k] = offset
        if sr.steps:
            last = sr.steps[-1]
            offset = (offset + last.x + last.dx).to("m")

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

        # --- x global: offset per stage + local x ---
        x_local = gp.x[i].to("m")
        x_global = (stage_offsets[k_stage] + x_local).to("m")

        # NOTE: all units kept as in original code, only ordering/selection changed
        row = {
            # 1) identifiers
            "stage_name": gp.stage_name[i],                      # stage name
            "i": i,                                              # global step index

            # 2) geometry / per-step
            "x[m]": x_global.magnitude,                          # GLOBAL x
            "dx[m]": gp.dx[i].to("m").magnitude,
            "qprime[MW/m]": gp.qprime[i].to("MW/m").magnitude,
            "UA_prime[MW/K/m]": gp.UA_prime[i].to("MW/K/m").magnitude,

            # 3) gas state
            "gas_P[Pa]": g.P.to("Pa").magnitude,
            "gas_T[°C]": g.T.to("degC").magnitude,
            "gas_h[kJ/kg]": g_h.to("kJ/kg").magnitude,

            # 4) water state (mirrored j)
            "water_P[Pa]": w.P.to("Pa").magnitude,
            "water_T[°C]": Tw.to("degC").magnitude,
            "water_h[kJ/kg]": w.h.to("kJ/kg").magnitude,

            # 5) phase / radiation flags
            "gas_eps[-]": gas_eps,
            "water_x[-]": _mag_or_nan(xq, ""),
            "boiling": bool(xq is not None),

            # 6) gas hydraulics & HTC
            "gas_V[m/s]": gas_V.to("m/s").magnitude,
            "Re_gas[-]": Re_gas,
            "h_gas[W/m^2/K]": gp.h_g[i].to("W/m^2/K").magnitude,

            # 7) water hydraulics & HTC
            "water_V[m/s]": (
                _mag_or_nan(water_V, "m/s") if isinstance(water_V, Q_) else float("nan")
            ),
            "Re_water[-]": Re_water,
            "h_water[W/m^2/K]": gp.h_c[i].to("W/m^2/K").magnitude,

            # 8) pressure drops (gas side, per step)
            "dP_fric[Pa]": gp.dP_fric[i].to("Pa").magnitude,
            "dP_minor[Pa]": gp.dP_minor[i].to("Pa").magnitude,
            "dP_total[Pa]": gp.dP_total[i].to("Pa").magnitude,

            # 9) gas transport / thermophysical
            "gas_cp[kJ/kg/K]": g_cp.to("kJ/kg/K").magnitude,
            "gas_mu[Pa*s]": g_mu.to("Pa*s").magnitude,
            "gas_k[W/m/K]": g_k.to("W/m/K").magnitude,
            "gas_rho[kg/m^3]": g_rho.to("kg/m^3").magnitude,

            # 10) water transport / thermophysical
            "water_cp[kJ/kg/K]": _mag_or_nan(w_cp, "kJ/kg/K"),
            "water_mu[Pa*s]": _mag_or_nan(w_mu, "Pa*s"),
            "water_k[W/m/K]": _mag_or_nan(w_k, "W/m/K"),
            "water_rho[kg/m^3]": _mag_or_nan(w_rho, "kg/m^3"),
        }

        rows.append(row)

    return pd.DataFrame(rows)


def summary_from_profile(gp: "GlobalProfile", combustion: CombustionResult | None = None,) -> tuple[list[dict], float, float]:
    rows = []
    Q_total = 0.0
    UA_total = 0.0
    Q_total_conv = 0.0
    Q_total_rad  = 0.0 
    dP_total_fric = 0.0
    dP_total_minor = 0.0
    dP_total_total = 0.0
    stack_T_C = None

    import itertools
    for k, grp in itertools.groupby(range(len(gp.x)), key=lambda i: gp.stage_index[i]):
        idxs = list(grp)
        name = gp.stage_name[idxs[0]]

        # integrals
        Q_stage = sum((gp.qprime[i] * gp.dx[i]).to("MW").magnitude for i in idxs)
        UA_stage = sum((gp.UA_prime[i] * gp.dx[i]).to("MW/K").magnitude for i in idxs)

        sr_stage = gp.stage_results[k]
        Q_stage_conv = sum((st.qprime_conv * st.dx).to("MW").magnitude for st in sr_stage.steps)
        Q_stage_rad  = sum((st.qprime_rad  * st.dx).to("MW").magnitude for st in sr_stage.steps)

        # ΔP sums
        dP_fric = sum(gp.dP_fric[i].to("Pa").magnitude for i in idxs)
        dP_minor = sum(gp.dP_minor[i].to("Pa").magnitude for i in idxs)
        dP_total = sum(gp.dP_total[i].to("Pa").magnitude for i in idxs)

        # endpoints along gas x in this stage
        g_in  = gp.gas[idxs[0]]
        g_out = gp.gas[idxs[-1]]

        # gas endpoints: P, T, h (sensible), in/out
        gas_in_T  = g_in.T.to("degC").magnitude
        gas_out_T = g_out.T.to("degC").magnitude

        gas_in_P  = g_in.P.to("Pa").magnitude
        gas_out_P = g_out.P.to("Pa").magnitude

        gas_in_h  = _gas.h_sensible(g_in.T,  g_in.P,  g_in.comp).to("kJ/kg").magnitude
        gas_out_h = _gas.h_sensible(g_out.T, g_out.P, g_out.comp).to("kJ/kg").magnitude

        # water endpoints: counter-current (water inlet at gas x=L)
        w_in  = gp.water[idxs[-1]]   # water inlet
        w_out = gp.water[idxs[0]]    # water outlet

        water_in_h  = w_in.h.to("kJ/kg").magnitude
        water_out_h = w_out.h.to("kJ/kg").magnitude

        water_in_P  = w_in.P.to("Pa").magnitude
        water_out_P = w_out.P.to("Pa").magnitude

        water_in_T  = WaterProps.T_from_Ph(w_in.P,  w_in.h).to("degC").magnitude
        water_out_T = WaterProps.T_from_Ph(w_out.P, w_out.h).to("degC").magnitude



        row = {
            "stage_index": k,
            "stage_name": name,
            "stage_kind": gp.stage_results[k].stage_kind,

            "Q_stage[MW]": Q_stage,
            "UA_stage[MW/K]": UA_stage,

            # gas endpoints
            "gas_in_P[Pa]": gas_in_P,
            "gas_in_T[°C]": gas_in_T,
            "gas_in_h[kJ/kg]": gas_in_h,
            "gas_out_P[Pa]": gas_out_P,
            "gas_out_T[°C]": gas_out_T,
            "gas_out_h[kJ/kg]": gas_out_h,

            # water endpoints (counter-current)
            "water_in_P[Pa]": water_in_P,
            "water_in_T[°C]": water_in_T,
            "water_in_h[kJ/kg]": water_in_h,
            "water_out_P[Pa]": water_out_P,
            "water_out_T[°C]": water_out_T,
            "water_out_h[kJ/kg]": water_out_h,

            # ΔP stage totals
            "ΔP_stage_fric[Pa]": dP_fric,
            "ΔP_stage_minor[Pa]": dP_minor,
            "ΔP_stage_total[Pa]": dP_total,

            # heat splits
            "Q_conv_stage[MW]": Q_stage_conv,
            "Q_rad_stage[MW]": Q_stage_rad,

            # boiler-level fields (blank at stage level)
            "η_direct[-]": "",
            "η_indirect[-]": "",
            "Q_total_useful[MW]": "",
            "Q_in_total[MW]": "",
            "P_LHV[MW]": "",
            "LHV_mass[kJ/kg]": "",
        }
        
        rows.append(row)

        Q_total += Q_stage
        UA_total += UA_stage
        Q_total_conv += Q_stage_conv
        Q_total_rad  += Q_stage_rad
        dP_total_fric  += dP_fric
        dP_total_minor += dP_minor
        dP_total_total += dP_total
        stack_T_C = gas_out_T

    # Global boiler useful duty (W)
    Q_useful = Q_total

    # Default values if no combustion info
    Q_in_total = None
    P_LHV_W = None
    LHV_mass_kJkg = None
    eta_direct = None
    eta_indirect = None

    if combustion is not None:
        # Q_in is currently in kW
        Q_in_total = combustion.Q_in.to("MW").magnitude

        # firing capacity based on LHV (kW) – try fuel_P_LHV if present, else LHV field
        if combustion.fuel_P_LHV is not None:
            P_LHV_W = combustion.fuel_P_LHV.to("MW").magnitude
        else:
            P_LHV_W = combustion.LHV.to("MW").magnitude

        # mass-based LHV if available
        if combustion.fuel_LHV_mass is not None:
            LHV_mass_kJkg = combustion.fuel_LHV_mass.to("kJ/kg").magnitude

        # Direct efficiency: useful duty / firing capacity (LHV basis)
        if P_LHV_W and P_LHV_W > 0.0:
            eta_direct = Q_useful / P_LHV_W

        # Approximate indirect efficiency: useful duty / total input heat
        if Q_in_total and Q_in_total > 0.0:
            eta_indirect = Q_useful / Q_in_total

    total_row = {
        "stage_index": "",
        "stage_name": "TOTAL_BOILER",
        "stage_kind": "",

        "Q_stage[MW]": Q_useful,
        "UA_stage[MW/K]": UA_total,

        # endpoints not meaningful at TOTAL_BOILER level
        "gas_in_P[Pa]": "",
        "gas_in_T[°C]": "",
        "gas_in_h[kJ/kg]": "",
        "gas_out_P[Pa]": "",
        "gas_out_T[°C]": "",
        "gas_out_h[kJ/kg]": "",
        "water_in_P[Pa]": "",
        "water_in_T[°C]": "",
        "water_in_h[kJ/kg]": "",
        "water_out_P[Pa]": "",
        "water_out_T[°C]": "",
        "water_out_h[kJ/kg]": "",

        "ΔP_stage_fric[Pa]": dP_total_fric,
        "ΔP_stage_minor[Pa]": dP_total_minor,
        "ΔP_stage_total[Pa]": dP_total_total,

        "stack_temperature[°C]": stack_T_C,

        "Q_conv_stage[MW]": Q_total_conv,
        "Q_rad_stage[MW]": Q_total_rad,

        "η_direct[-]": eta_direct if eta_direct is not None else "",
        "η_indirect[-]": eta_indirect if eta_indirect is not None else "",
        "Q_total_useful[MW]": Q_useful,
        "Q_in_total[MW]": Q_in_total if Q_in_total is not None else "",
        "P_LHV[MW]": P_LHV_W if P_LHV_W is not None else "",
        "LHV_mass[kJ/kg]": LHV_mass_kJkg if LHV_mass_kJkg is not None else "",
    }


    rows.append(total_row)

    return rows, Q_total, UA_total

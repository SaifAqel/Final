from math import pi, log10
from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage
from props import GasProps, WaterProps

sigma = Q_(5.670374419e-8, "W/m^2/K^4")

P_CRIT_WATER = Q_(22.064, "MPa")
MW_WATER = 18.01528

_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")

def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})
def cp_water(w: WaterStream) -> Q_: return WaterProps.cp_from_Ph(w.P, w.h)
def T_water(w: WaterStream) -> Q_:  return WaterProps.T_from_Ph(w.P, w.h)

def _Aprime_hot(spec: dict) -> Q_:
    if "inner_diameter" not in spec:
        raise KeyError("spec['inner_diameter'] required")
    Di = spec["inner_diameter"].to("m")
    nt = int(spec.get("tubes_number", Q_(1.0, "dimensionless")).to("dimensionless").magnitude)
    return (pi * Di * nt).to("m")

def _Aprime_cold(spec: dict) -> Q_:
    if "cold_Pi" in spec: return spec["cold_Pi"].to("m")
    if "cold_Di" in spec: return (pi * spec["cold_Di"].to("m")).to("m")
    raise KeyError("cold-side geometry required: 'cold_Pi' or 'cold_Di'")

def _gas_velocity(g: GasStream, spec: dict) -> Q_:
    Di = spec["inner_diameter"].to("m")
    nt = int(spec.get("tubes_number", Q_(1.0, "dimensionless")).to("dimensionless").magnitude)
    Aflow = nt * (pi * (Di/2)**2)
    rho = _gas.rho(g.T, g.P, g.comp or {})
    return (g.mass_flow / (rho * Aflow)).to("m/s")

def _cold_velocity(w: WaterStream, spec: dict, Tw: Q_) -> Q_:
    if "cold_Ai" in spec:
        Aflow = spec["cold_Ai"].to("m^2")
    elif "cold_Di" in spec:
        Di = spec["cold_Di"].to("m"); Aflow = pi * (Di/2)**2
    else:
        raise KeyError("cold-side area required: 'cold_Ai' or 'cold_Di'")
    rho = WaterProps.rho_from_PT(w.P, Tw)
    return (w.mass_flow / (rho * Aflow)).to("m/s")

def _h_gas_conv(g: GasStream, spec: dict) -> Q_:
    Di = spec["inner_diameter"].to("m")
    V = _gas_velocity(g, spec)
    rho = _gas.rho(g.T, g.P, g.comp or {})
    mu  = _gas.mu(g.T, g.P, g.comp or {})
    k   = _gas.k(g.T, g.P, g.comp or {})
    Re = (rho * V * Di / mu).to("dimensionless").magnitude
    Pr = (cp_gas(g) * mu / k).to("dimensionless").magnitude
    Nu = 3.66 if Re < 2300 else 0.023 * (Re**0.8) * (Pr**0.4)
    return (Nu * k / Di).to("W/m^2/K")

def _h_gas_rad(g: GasStream, spec: dict) -> Q_:
    if "eps_in" not in spec:
        raise KeyError("spec['eps_in'] required for radiation")
    eps = spec["eps_in"].to("dimensionless").magnitude
    Tf = g.T.to("K")
    return (4.0 * eps * sigma * Tf**3).to("W/m^2/K")

def _is_boiling(P: Q_, h: Q_) -> bool:
    hf = WaterProps.h_f(P).to("J/kg")
    hg = WaterProps.h_g(P).to("J/kg")
    hJ = h.to("J/kg")
    return hf.magnitude <= hJ.magnitude <= hg.magnitude

def _h_water_singlephase(w: WaterStream, spec: dict, Tw: Q_) -> Q_:
    if "cold_Di" in spec:
        Dh = spec["cold_Di"].to("m")
    elif "cold_Ai" in spec and "cold_Pi" in spec:
        Dh = (4 * spec["cold_Ai"] / spec["cold_Pi"]).to("m")
    else:
        raise KeyError("cold-side hydraulic diameter info required")

    V   = _cold_velocity(w, spec, Tw)
    rho = WaterProps.rho_from_PT(w.P, Tw)
    mu  = WaterProps.mu_from_PT(w.P, Tw)
    k   = WaterProps.k_from_PT(w.P, Tw)
    cpw = WaterProps.cp_from_PT(w.P, Tw)

    Re = (rho * V * Dh / mu).to("dimensionless").magnitude
    Pr = (cpw * mu / k).to("dimensionless").magnitude
    Nu = 3.66 if Re < 2300 else 0.023 * (Re**0.8) * (Pr**0.3)
    return (Nu * k / Dh).to("W/m^2/K")

def _h_water_boil_cooper(P: Q_, qpp: Q_, Rp: Q_) -> Q_:
    p_r = (P.to("MPa") / P_CRIT_WATER).magnitude
    Rp_um = Rp.to("micrometer").magnitude
    if Rp_um <= 0:
        raise ValueError("surface roughness must be positive")
    q_kWm2 = qpp.to("kW/m^2").magnitude
    if q_kWm2 <= 0:
        raise ValueError("heat flux must be positive")
    h_kWm2K = 55.0 * (p_r**0.12) * ((-log10(Rp_um))**-0.55) * (MW_WATER**-0.5) * (q_kWm2**0.67)
    return Q_(h_kWm2K, "kW/m^2/K").to("W/m^2/K")

@trace_calls(values=True)
def ua_per_m(g, w, stage):
    spec = stage.spec
    Agh = _Aprime_hot(spec)
    Ac  = _Aprime_cold(spec)
    Tw  = WaterProps.T_from_Ph(w.P, w.h)

    h_g = _h_gas_conv(g, spec) + _h_gas_rad(g, spec)

    if "boiling_mode" in spec:
        boiling = str(spec["boiling_mode"]).lower() in {"1","true","yes","boil","nukiyama"}
    else:
        boiling = _is_boiling(w.P, w.h)

    if boiling:
        if "roughness_out" not in spec:
            raise KeyError("spec['roughness_out'] required for boiling")
        qpp_guess = spec.get("qpp_guess", Q_(1e5, "W/m^2"))  # first cell only; solver will overwrite
        h_c = _h_water_boil_cooper(w.P, qpp_guess, spec["roughness_out"])
    else:
        h_c = _h_water_singlephase(w, spec, Tw)

    # thermal resistances and UA'
    t_fg = spec.get("foul_t_in", Q_(0.0, "m"))
    k_fg = spec.get("foul_k_in", Q_(1e9, "W/m/K"))
    t_fc = spec.get("foul_t_out", Q_(0.0, "m"))
    k_fc = spec.get("foul_k_out", Q_(1e9, "W/m/K"))
    if "wall_t" not in spec or "wall_k" not in spec:
        raise KeyError("spec['wall_t'] and spec['wall_k'] required")
    Awm = Agh

    Rg  = (1.0 / (h_g * Agh)).to("K*m/W")
    Rfg = (t_fg / (k_fg * Agh)).to("K*m/W")
    Rw  = (spec["wall_t"].to("m") / (spec["wall_k"].to("W/m/K") * Awm)).to("K*m/W")
    Rfc = (t_fc / (k_fc * Ac)).to("K*m/W")
    Rc  = (1.0 / (h_c * Ac)).to("K*m/W")
    UA_prime = (1.0 / (Rg + Rfg + Rw + Rfc + Rc)).to("W/K/m")

    # define physical heat flux ALWAYS (prevents UnboundLocalError)
    dT = (g.T - Tw).to("K")
    qprime = (UA_prime * dT).to("W/m")   # per length
    qpp_phys = (qprime / Ac).to("W/m^2")

    return UA_prime, qpp_phys

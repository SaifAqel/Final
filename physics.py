# --- inside physics.py ---
from math import pi
from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage
from props import GasProps, WaterProps

sigma = Q_(5.670374419e-8, "W/m^2/K^4")  # Stefan–Boltzmann

_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")

def cp_gas(g: GasStream) -> Q_:
    return _gas.cp(g.T, g.P, g.comp or {})

def cp_water(w: WaterStream) -> Q_:
    return WaterProps.cp_from_Ph(w.P, w.h)

def T_water(w: WaterStream) -> Q_:
    return WaterProps.T_from_Ph(w.P, w.h)

# ---------- geometry helpers (per-length areas) ----------
def _Aprime_hot(spec: dict) -> Q_:
    if "hot_Di" in spec:
        Di = spec["hot_Di"].to("m")
        nt = int(spec.get("hot_ntubes", Q_(1.0, "dimensionless")).magnitude)
        return (pi * Di * nt).to("m")  # m^2 per m length
    # fallback: if only perimeter is implicit, assume 1 m perimeter
    return Q_(1.0, "m")

def _Aprime_cold(spec: dict) -> Q_:
    if "cold_Pi" in spec:
        return spec["cold_Pi"].to("m")
    if "cold_Di" in spec:
        return (pi * spec["cold_Di"].to("m")).to("m")
    return Q_(1.0, "m")

def _Aprime_wall_mean(spec: dict) -> Q_:
    # crude mean between sides
    return 0.5 * (_Aprime_hot(spec) + _Aprime_cold(spec))

# ---------- velocities and properties ----------
def _gas_velocity(g: GasStream, spec: dict) -> Q_:
    # Flow area on hot side
    if "hot_Di" in spec:
        Di = spec["hot_Di"].to("m")
        nt = int(spec.get("hot_ntubes", Q_(1.0, "dimensionless")).magnitude)
        Aflow = nt * (pi * (Di/2)**2)
    else:
        # assume 1 m^2 if unknown; keeps units consistent
        Aflow = Q_(1.0, "m^2")
    rho = _gas.rho(g.T, g.P, g.comp or {})
    return (g.mass_flow / (rho * Aflow)).to("m/s")

def _cold_velocity(w: WaterStream, spec: dict, Tw: Q_) -> Q_:
    if "cold_Ai" in spec:
        Aflow = spec["cold_Ai"].to("m^2")
    elif "cold_Di" in spec:
        Di = spec["cold_Di"].to("m"); Aflow = pi * (Di/2)**2
    else:
        Aflow = Q_(1.0, "m^2")
    rho = WaterProps.rho_from_PT(w.P, Tw)  # was rho_from_Ph
    return (w.mass_flow / (rho * Aflow)).to("m/s")


# ---------- h correlations ----------
def _h_gas_conv(g: GasStream, spec: dict) -> Q_:
    # Dittus–Boelter for turbulent internal flow
    Di = spec.get("hot_Di")
    if not Di:
        return Q_(10.0, "W/m^2/K")  # benign fallback
    Di = Di.to("m")
    V = _gas_velocity(g, spec)
    rho = _gas.rho(g.T, g.P, g.comp or {})
    mu = _gas.mu(g.T, g.P, g.comp or {})
    k  = _gas.k(g.T, g.P, g.comp or {})
    Re = (rho * V * Di / mu).to("dimensionless").magnitude
    Pr = (cp_gas(g) * mu / k).to("dimensionless").magnitude
    if Re < 2300:
        # simple laminar tube correlation Nu=3.66
        Nu = 3.66
    else:
        # Dittus–Boelter: Nu = 0.023 Re^0.8 Pr^n (n=0.4 heating)
        Nu = 0.023 * (Re**0.8) * (Pr**0.4)
    h = (Nu * k / Di).to("W/m^2/K")
    return h

def _h_gas_rad(g: GasStream, spec: dict) -> Q_:
    eps = spec.get("hot_eps", Q_(0.8, "dimensionless")).to("dimensionless").magnitude
    Tf = g.T.to("K")
    return (4.0 * eps * sigma * Tf**3).to("W/m^2/K")

def _is_boiling(P: Q_, h: Q_) -> bool:
    hf = WaterProps.h_f(P).to("J/kg")
    hg = WaterProps.h_g(P).to("J/kg")
    hJ = h.to("J/kg")
    return hf.magnitude <= hJ.magnitude <= hg.magnitude

def _h_water_singlephase(w: WaterStream, spec: dict, Tw: Q_) -> Q_:
    # Dittus–Boelter for water duct
    # pick hydraulic diameter
    if "cold_Di" in spec:
        Dh = spec["cold_Di"].to("m")
    elif "cold_Ai" in spec and "cold_Pi" in spec:
        Dh = (4 * spec["cold_Ai"] / spec["cold_Pi"]).to("m")
    else:
        Dh = Q_(0.02, "m")

    # use PT-based properties to avoid two-phase μ=None
    # nudge Tw slightly into liquid region if at saturation
    Tsat = WaterProps.Tsat(w.P)
    Tw_eff = (Tw - Q_(0.5, "K")) if abs((Tw - Tsat).to("K").magnitude) < 1e-6 else Tw

    V   = _cold_velocity(w, spec, Tw_eff)
    rho = WaterProps.rho_from_PT(w.P, Tw_eff)
    mu  = WaterProps.mu_from_PT(w.P, Tw_eff)
    k   = WaterProps.k_from_PT(w.P, Tw_eff)
    cpw = WaterProps.cp_from_PT(w.P, Tw_eff)

    Re = (rho * V * Dh / mu).to("dimensionless").magnitude
    Pr = (cpw * mu / k).to("dimensionless").magnitude
    Nu = 3.66 if Re < 2300 else 0.023 * (Re**0.8) * (Pr**0.3)
    return (Nu * k / Dh).to("W/m^2/K")

def _h_water_boil(w: WaterStream, spec: dict, Tw: Q_) -> Q_:
    # Placeholder hook; replace with Rohsenow/Chen/Gungor-Winterton as needed
    # Default to a high effective h to mimic nucleate boiling regime
    return Q_(5000.0, "W/m^2/K")

# ---------- overall UA per meter ----------
@trace_calls(values=True)
def ua_per_m(g: GasStream, w: WaterStream, stage: HXStage) -> Q_:
    spec = stage.spec
    # per-length areas
    Agh = _Aprime_hot(spec)
    Ac  = _Aprime_cold(spec)
    Awm = _Aprime_wall_mean(spec)

    # film temperature on cold side
    Tw = WaterProps.T_from_Ph(w.P, w.h)

    # gas-side h = convection + radiation
    h_g_conv = _h_gas_conv(g, spec)
    h_g_rad  = _h_gas_rad(g, spec)
    h_g = h_g_conv + h_g_rad

    flag = str(stage.spec.get("water_boil", "")).lower() in {"1","true","yes","nukiyama","boil"}
    if flag or _is_boiling(w.P, w.h):
        h_c = _h_water_boil(w, stage.spec, Tw)
    else:
        h_c = _h_water_singlephase(w, stage.spec, Tw)

    # fouling and wall data (optional, else zero)
    t_fg = spec.get("hot_foul_t", Q_(0.0, "m")).to("m")
    k_fg = spec.get("hot_foul_k", Q_(1e9, "W/m/K")).to("W/m/K")  # avoid div-by-zero
    t_fc = spec.get("cold_foul_t", Q_(0.0, "m")).to("m")
    k_fc = spec.get("cold_foul_k", Q_(1e9, "W/m/K")).to("W/m/K")
    t_w  = spec.get("hot_wall_t", spec.get("cold_wall_t", Q_(0.0, "m"))).to("m")
    k_w  = spec.get("hot_wall_k", spec.get("cold_wall_k", Q_(1e9, "W/m/K"))).to("W/m/K")

    # series resistances per length [K/W per m]
    Rg  = (1.0 / (h_g * Agh)).to("K*m/W")
    Rfg = (t_fg / (k_fg * Agh)).to("K*m/W")
    Rw  = (t_w  / (k_w  * Awm)).to("K*m/W")
    Rfc = (t_fc / (k_fc * Ac )).to("K*m/W")
    Rc  = (1.0 / (h_c * Ac )).to("K*m/W")

    Rtot = Rg + Rfg + Rw + Rfc + Rc
    UA_prime = (1.0 / Rtot).to("W/K/m")
    return UA_prime

# keep these wrappers working but delegate to the new UA
def heat_rate_per_length(g: GasStream, w: WaterStream, stage: HXStage) -> dict:
    UA = ua_per_m(g, w, stage)
    qprime = UA * (g.T - T_water(w))  # W/m
    return {"qprime": qprime, "UA_per_m": UA}

def rhs(g: GasStream, w: WaterStream, qprime: Q_, stage: HXStage) -> dict:
    dTgdx = - qprime / (g.mass_flow * cp_gas(g))
    dhwdx = + qprime / w.mass_flow
    dpgdx = Q_(0.0, "Pa/m")
    return {"dTgdx": dTgdx, "dhwdx": dhwdx, "dpgdx": dpgdx}

from math import pi, log10, log
from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage
from props import GasProps, WaterProps
from water_htc import water_htc

sigma = Q_(5.670374419e-8, "W/m^2/K^4")

P_CRIT_WATER = Q_(22.064, "MPa")
MW_WATER = 18.01528

_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")

def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})
def cp_water(w: WaterStream) -> Q_: return WaterProps.cp_from_Ph(w.P, w.h)
def T_water(w: WaterStream) -> Q_:  return WaterProps.T_from_Ph(w.P, w.h)

def _gas_velocity(g: GasStream, spec: dict) -> Q_:
    Di = spec["inner_diameter"].to("m")
    nt = int(spec.get("tubes_number", Q_(1.0, "dimensionless")).to("dimensionless").magnitude)
    Aflow = nt * (pi * (Di/2)**2)
    rho = _gas.rho(g.T, g.P, g.comp or {})
    return (g.mass_flow / (rho * Aflow)).to("m/s")

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

def _h_gas_rad(g: GasStream, spec: dict, Tgw) -> Q_:
    eps = spec["eps_in"].to("dimensionless").magnitude
    Tf = 0.5 * (g.T - Tgw).to("kelvin")
    return (4.0 * eps * sigma * Tf**3).to("W/m^2/K")

def gas_htc(g, spec, T_gw) -> Q_:
    return _h_gas_conv(g, spec) + _h_gas_rad(g, spec, T_gw)

def fouling_resistances(spec: dict) -> tuple[Q_, Q_]:
    di = spec["inner_diameter"].to("m")
    wall_t = spec["wall_t"].to("m")
    tfi = spec["foul_t_in"].to("m")
    kfi = spec["foul_k_in"].to("W/m/K")
    tfo = spec["foul_t_out"].to("m")
    kfo = spec["foul_k_out"].to("W/m/K")
    Rfi = (log(di / (di - tfi)) / (2 * pi * kfi)).to("K*m/W")
    Rfo = (log((di + wall_t + tfo) / (di + wall_t)) / (2 * pi * kfo)).to("K*m/W")
    return Rfi, Rfo

def wall_resistance(spec: dict) -> Q_:
    di = spec["inner_diameter"].to("m")
    wall_t = spec["wall_t"].to("m")
    k  = spec["wall_k"].to("W/m/K")
    return (log((di + wall_t)/ di) / (2 * pi * k)).to("K*m/W")

@trace_calls(values=True)
def ua_per_m(g: GasStream, w:WaterStream, stage: HXStage, T_gw: Q_, T_ww: Q_, qpp: Q_):
    spec = stage.spec
    Pg = stage.spec["inner_perimeter"]
    Pw = stage.spec["cold_wet_P"]
    Tw  = WaterProps.T_from_Ph(w.P, w.h)

    h_g = gas_htc(g, spec, T_gw)
    h_c, boiling = water_htc(w, stage, T_ww, qpp)

    Rg = (1.0 / (h_g * Pg)).to("K*m/W")
    Rfg, Rfc = fouling_resistances(spec)
    Rw = wall_resistance(spec)  # same area used previously (Awm = Agh)
    Rc = (1.0 / (h_c * Pw)).to("K*m/W")

    UA_prime = (1.0 / (Rg + Rfg + Rw + Rfc + Rc)).to("W/K/m")

    # physical heat flux
    dT = (g.T - Tw).to("K")
    qprime = (UA_prime * dT).to("W/m")        # per length
    qpp_phys = (qprime / Pw).to("W/m^2")

    return UA_prime, qpp_phys, boiling

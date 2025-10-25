from math import pi, log
from units import Q_
from models import GasStream
from props import GasProps

sigma = Q_(5.670374419e-8, "W/m^2/K^4")

P_CRIT_WATER = Q_(22.064, "MPa")
MW_WATER = 18.01528

_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")

def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})

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
    Tf = 0.5 * (g.T + Tgw).to("kelvin")
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

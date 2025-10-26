from units import Q_
from models import GasStream
from props import GasProps
import numpy as np
from typing import Dict, Any

_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")

def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})

def _mean_beam_length(spec: dict) -> Q_:
    Dh = spec["hot_Dh"].to("m")
    return (0.9 * Dh).to("m")

def _gas_partials(g: GasStream) -> tuple[float, float]:
    P = g.P.to("Pa").magnitude
    y = g.comp or {}
    yH2O = y.get("H2O", Q_(0.0, "")).to("").magnitude
    yCO2 = y.get("CO2", Q_(0.0, "")).to("").magnitude
    return yH2O * P, yCO2 * P

# Smith, Shen & Friedman (ASME J. Heat Transfer, 1982)
# 4-gray model coefficients for CO2–H2O flue gas at ~1 atm
_A = np.array([0.434, 0.313, 0.180, 0.073])        # weighting factors a_j
_K = np.array([0.0, 2.3, 11.6, 30.4])              # base absorption coeffs [1/m at 1 atm]

def emissivity(T_K: float, pH2O_Pa: float, pCO2_Pa: float, L_m: float) -> float:
    p_ratio = (pH2O_Pa + pCO2_Pa) / 101325.0
    tau = _K * p_ratio * L_m
    eps = 1.0 - np.sum(_A * np.exp(-tau))
    return float(np.clip(eps, 0.0, 1.0))

def h_rad(Tfilm_K: float, eps_g: float, F: float = 1.0) -> float:
    sigma = 5.670374419e-8
    return 4.0 * sigma * F * eps_g * Tfilm_K**3

def h_gas_rad_smith(T_K: float, pH2O_Pa: float, pCO2_Pa: float, L_m: float, Twall_K: float, F: float = 1.0) -> float:
    eps = emissivity(T_K, pH2O_Pa, pCO2_Pa, L_m)
    Tfilm = 0.5 * (T_K + Twall_K)
    return h_rad(Tfilm, eps, F)

def h_gas_rad(g: GasStream, spec: Dict[str, Any], Tgw: Q_) -> Q_:
    Tg = g.T
    pH2O, pCO2 = _gas_partials(g)
    Lb = _mean_beam_length(spec).to("m").magnitude
    F = float(1)
    h_val = h_gas_rad_smith(Tg.to("K").magnitude, pH2O, pCO2, Lb, Tgw.to("K").magnitude, F)
    return Q_(h_val, "W/m^2/K")

def _gas_velocity(g: GasStream, spec: dict) -> Q_:
    Aflow = spec["hot_flow_A"]
    rho = _gas.rho(g.T, g.P, g.comp)
    return (g.mass_flow / (rho * Aflow)).to("m/s")

def h_gas_conv(g: GasStream, spec: dict) -> Q_:
    Di = spec["inner_diameter"].to("m")
    V = _gas_velocity(g, spec)
    rho = _gas.rho(g.T, g.P, g.comp)
    mu  = _gas.mu(g.T, g.P, g.comp)
    k   = _gas.k(g.T, g.P, g.comp)
    Re = (rho * V * Di / mu).to("dimensionless").magnitude
    Pr = (cp_gas(g) * mu / k).to("dimensionless").magnitude
    Nu = 3.66 if Re < 2300 else 0.023 * (Re**0.8) * (Pr**0.4)
    return (Nu * k / Di).to("W/m^2/K")

def gas_htc(g, spec, Tgw) -> Q_:
    return h_gas_conv(g, spec) + h_gas_rad(g, spec, Tgw)
# physics.py
from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage
from props import GasProps, WaterProps

# Cantera mechanism: make sure phase name matches your YAML
_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")
@trace_calls
def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})
@trace_calls
def cp_water(w: WaterStream) -> Q_: return WaterProps.cp_from_Ph(w.P, w.h)
@trace_calls
def T_water(w: WaterStream) -> Q_:  return WaterProps.T_from_Ph(w.P, w.h)

@trace_calls()
def heat_rate_per_length(g: GasStream, w: WaterStream, stage: HXStage) -> dict:
    UA_per_m: Q_ = stage.hot.get("UA_per_m", Q_(100.0, "W/K/m"))  # still placeholder UA′
    qprime = UA_per_m * (g.T - T_water(w))                        # W/m
    return {"qprime": qprime, "UA_per_m": UA_per_m}

@trace_calls()
def rhs(g: GasStream, w: WaterStream, qprime: Q_, stage: HXStage) -> dict:
    dTgdx = - qprime / (g.mass_flow * cp_gas(g))  # K/m
    dhwdx = + qprime / w.mass_flow                # (J/kg)/m
    dpgdx = Q_(0.0, "Pa/m")                       # add friction later
    return {"dTgdx": dTgdx, "dhwdx": dhwdx, "dpgdx": dpgdx}

# physics.py
from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage
from props import GasProps, WaterProps

@trace_calls(values=True)
def ua_per_m(stage: HXStage) -> Q_:
    # stage-specific constants; fall back to default
    table = {
        "HX_1": Q_(180.0, "W/K/m"),
        "HX_2": Q_(140.0, "W/K/m"),
        "HX_3": Q_(220.0, "W/K/m"),
        "HX_4": Q_(140.0, "W/K/m"),
        "HX_5": Q_(200.0, "W/K/m"),
        "HX_6": Q_(160.0, "W/K/m"),
    }
    return stage.spec.get("UA_per_m") or table.get(stage.name, Q_(150.0, "W/K/m"))

# Cantera mechanism: make sure phase name matches your YAML
_gas = GasProps(mech_path="config/flue_cantera.yaml", phase="gas_mix")
@trace_calls(values=True)
def cp_gas(g: GasStream) -> Q_:    return _gas.cp(g.T, g.P, g.comp or {})
@trace_calls(values=True)
def cp_water(w: WaterStream) -> Q_: return WaterProps.cp_from_Ph(w.P, w.h)
@trace_calls(values=True)
def T_water(w: WaterStream) -> Q_:  return WaterProps.T_from_Ph(w.P, w.h)

@trace_calls(values=True)
def heat_rate_per_length(g: GasStream, w: WaterStream, stage: HXStage) -> dict:
    UA = ua_per_m(stage)
    qprime = UA * (g.T - T_water(w))                 # W/m
    return {"qprime": qprime, "UA_per_m": UA}

@trace_calls(values=True)
def rhs(g: GasStream, w: WaterStream, qprime: Q_, stage: HXStage) -> dict:
    dTgdx = - qprime / (g.mass_flow * cp_gas(g))  # K/m
    dhwdx = + qprime / w.mass_flow                # (J/kg)/m
    dpgdx = Q_(0.0, "Pa/m")                       # add friction later
    return {"dTgdx": dTgdx, "dhwdx": dhwdx, "dpgdx": dpgdx}

from logging_utils import trace_calls
from units import Q_, ureg
from models import GasStream, WaterStream, HXStage

CP_GAS   = Q_(1200.0, "J/kg/K")  # placeholder
CP_WATER = Q_(4200.0, "J/kg/K")  # placeholder
TREF_WATER = Q_(300.0, "K")

def water_T_from_h(h: Q_) -> Q_:
    return TREF_WATER + h / CP_WATER

@trace_calls()
def heat_rate_per_length(g: GasStream, w: WaterStream, stage: HXStage) -> dict:
    UA_per_m: Q_ = stage.hot.get("UA_per_m", Q_(100.0, "W/K/m"))
    Tw = water_T_from_h(w.h)
    qprime = UA_per_m * (g.T - Tw)       # W/m
    return {"qprime": qprime, "UA_per_m": UA_per_m}

@trace_calls()
def rhs(g: GasStream, w: WaterStream, qprime: Q_, stage: HXStage) -> dict:
    dTgdx = - qprime / (g.mass_flow * CP_GAS)   # K/m
    dhwdx = + qprime / w.mass_flow              # (J/kg)/m
    dpgdx = Q_(0.0, "Pa/m")                     # placeholder
    return {"dTgdx": dTgdx, "dhwdx": dhwdx, "dpgdx": dpgdx}

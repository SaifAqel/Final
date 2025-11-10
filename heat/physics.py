from math import pi, log
from common.units import Q_

def _nt(spec) -> float:
    Nt = spec.get("tubes_number", None)
    return (Nt.to("").magnitude if Nt is not None else 1.0)

def fouling_resistances(spec: dict) -> tuple[Q_, Q_]:
    di = spec["inner_diameter"].to("m")
    do = spec["outer_diameter"].to("m")
    tfi = spec["foul_t_in"].to("m")
    kfi = spec["foul_k_in"].to("W/m/K")
    tfo = spec["foul_t_out"].to("m")
    kfo = spec["foul_k_out"].to("W/m/K")
    Di_new = di - 2*tfi
    Rfi = log(di/Di_new)/(2*pi*kfi)
    do_new = do + 2*tfo
    Rfo = log(do_new/do)/(2*pi*kfo)

    Nt = _nt(spec)
    return (Rfi / Nt).to("K*m/W"), (Rfo / Nt).to("K*m/W")

def wall_resistance(spec: dict) -> Q_:
    di = spec["inner_diameter"].to("m")
    do = spec["outer_diameter"].to("m")
    k  = spec["wall_k"].to("W/m/K")
    R = log(do/di) / (2*pi*k)
    Nt = _nt(spec)
    return (R / Nt).to("K*m/W")
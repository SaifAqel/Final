from math import pi, log
from common.units import Q_

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
    return Rfi, Rfo

def wall_resistance(spec: dict) -> Q_:
    di = spec["inner_diameter"].to("m")
    do = spec["outer_diameter"].to("m")
    k  = spec["wall_k"].to("W/m/K")
    return (log(do/di) / (2*pi*k)).to("K*m/W")

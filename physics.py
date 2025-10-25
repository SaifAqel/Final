from math import pi, log
from units import Q_

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

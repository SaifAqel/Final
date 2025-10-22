from dataclasses import dataclass
from typing import Dict, Any
from units import Q_

@dataclass
class GasStream:
    mass_flow: Q_     # kg/s
    T: Q_             # K
    P: Q_             # Pa
    comp: Dict[str, Q_] | None = None
    stage: str = "-"

@dataclass
class WaterStream:
    mass_flow: Q_     # kg/s
    h: Q_             # J/kg
    P: Q_             # Pa
    stage: str = "-"

@dataclass
class HXStage:
    name: str
    kind: str
    spec: Dict[str, Any]

@dataclass
class Drum:
    Di: Q_          # inner diameter, m
    L: Q_     
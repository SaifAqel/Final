from dataclasses import dataclass
from typing import Dict, Iterable
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
    spec: Dict[str, Q_]        # free-form bag; must be Q_

    def __iter__(self) -> Iterable[tuple[str, Q_]]:
        yield from (("L", self.L),)

@dataclass
class Drum:
    Di: Q_          # inner diameter, m
    L: Q_     
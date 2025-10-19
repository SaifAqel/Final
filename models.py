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
    L: Q_                      # m
    hot: Dict[str, Q_]         # free-form bag; must be Q_
    cold: Dict[str, Q_]

    def __iter__(self) -> Iterable[tuple[str, Q_]]:
        yield from (("L", self.L),)

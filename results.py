from dataclasses import dataclass
from units import Q_
from typing import List
from models import GasStream, WaterStream

from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class StepResult:
    i: Optional[int]
    x: Optional[Q_]; dx: Optional[Q_]
    gas: GasStream
    water: WaterStream
    Tgw: Q_; Tww: Q_
    UA_prime: Q_
    qprime: Q_
    boiling: bool


@dataclass(frozen=True)
class StageResult:
    stage_name: str
    stage_kind: str
    steps: List[StepResult]
    Q_stage: Q_
    UA_stage: Q_

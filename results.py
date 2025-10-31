from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Types provided by the caller's codebase
from units import Q_
from models import GasStream, WaterStream


@dataclass(frozen=True)
class StepResult:
    i: Optional[int]
    x: Optional[Q_]
    dx: Optional[Q_]
    gas: GasStream
    water: WaterStream
    Tgw: Q_
    Tww: Q_
    UA_prime: Q_
    qprime: Q_
    boiling: bool
    stage_name: Optional[str] = None
    stage_index: Optional[int] = None


@dataclass(frozen=True)
class StageResult:
    stage_name: str
    stage_kind: str
    steps: List[StepResult]
    Q_stage: Q_
    UA_stage: Q_
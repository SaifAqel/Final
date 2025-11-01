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
    h_g: Q_          # <-- add
    h_c: Q_          # <-- add
    stage_name: Optional[str] = None
    stage_index: Optional[int] = None


@dataclass(frozen=True)
class StageResult:
    stage_name: str
    stage_kind: str
    steps: List[StepResult]
    Q_stage: Q_
    UA_stage: Q_

# results.py
from dataclasses import dataclass
from typing import List
from units import Q_
from models import GasStream, WaterStream

@dataclass(frozen=True)
class GlobalProfile:
    x: List[Q_]    
    dx: List[Q_]          # global gas coordinate from 0 → sum L_stage
    gas: List[GasStream]     # gas(x)
    water: List[WaterStream] # water(x) remapped from L−x
    stage_index: List[int]   # owning stage per point
    stage_name: List[str]    # owning stage name per point
    qprime: List[Q_]         # q′ at gas x (from StepResult)
    UA_prime: List[Q_]       # UA′ at gas x (from StepResult)
    h_g: List[Q_]      # <-- add
    h_c: List[Q_]      # <-- add

# results.py (or profiles.py)
from typing import List
from units import Q_

def build_global_profile(stage_results: List[StageResult]) -> GlobalProfile:
    x_glob: List[Q_] = []
    dx_glob: List[Q_] = []
    gas_glob: List[GasStream] = []
    water_glob: List[WaterStream] = []
    idx_glob: List[int] = []
    name_glob: List[str] = []
    qp_glob: List[Q_] = []          # <-- fix
    UA_glob: List[Q_] = []  
    h_g_glob: List[Q_] = []     # <-- add
    h_c_glob: List[Q_] = []     # <-- add

    x0 = Q_(0.0, "m")
    for sr in stage_results:
        steps = sr.steps
        n = len(steps)
        if n == 0: continue
        for i in range(n):
            s = steps[i]
            x_glob.append((x0 + s.x).to("m"))
            dx_glob.append(s.dx.to("m"))                # <-- new
            gas_glob.append(s.gas)
            water_glob.append(steps[n-1-i].water)
            idx_glob.append(s.stage_index)
            name_glob.append(sr.stage_name)
            qp_glob.append(s.qprime)
            UA_glob.append(s.UA_prime)
            h_g_glob.append(s.h_g)    # <-- add
            h_c_glob.append(s.h_c)    # <-- add
        L_stage = (steps[-1].x + steps[-1].dx).to("m")
        x0 = (x0 + L_stage).to("m")

    return GlobalProfile(
        x=x_glob, dx=dx_glob,
        gas=gas_glob, water=water_glob,
        stage_index=idx_glob, stage_name=name_glob,
        qprime=qp_glob, UA_prime=UA_glob,
        h_g=h_g_glob, h_c=h_c_glob  
    )

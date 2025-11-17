from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Sequence, Tuple
from common.units import Q_
from common.models import GasStream, WaterStream
from pathlib import Path
import pandas as pd
from common.results import GlobalProfile, CombustionResult
from heat.postproc import profile_to_dataframe, summary_from_profile


@dataclass(frozen=True)
class CombustionResult:
    LHV: Q_
    Q_in: Q_
    T_ad: Q_
    flue: GasStream
    fuel_LHV_mass: Q_ | None = None   # e.g. kJ/kg
    fuel_P_LHV: Q_ | None = None 

# ---------------- Streams are imported from common.models ----------------

@dataclass(frozen=True)
class StepResult:
    # marching index and geometry
    i: int
    x: Q_
    dx: Q_
    # stream snapshots at start of step
    gas: object     # GasStream
    water: object   # WaterStream
    # wall state
    Tgw: Q_
    Tww: Q_
    # local per-length metrics
    UA_prime: Q_
    qprime: Q_
    boiling: bool
    # HTCs for diagnostics
    h_g: Q_
    h_c: Q_
    # ---- appended fields (binary compatible via defaults) ----
    stage_name: str = ""
    stage_index: int = -1
    dP_fric: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))
    dP_minor: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))
    dP_total: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))

    qprime_conv: Q_ = field(default_factory=lambda: Q_(0.0, "W/m"))
    qprime_rad: Q_ = field(default_factory=lambda: Q_(0.0, "W/m"))

@dataclass(frozen=True)
class StageResult:
    stage_name: str
    stage_kind: str
    steps: Sequence[StepResult]
    Q_stage: Q_
    UA_stage: Q_
    # ---- appended fields (binary compatible via defaults) ----
    dP_stage_fric: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))
    dP_stage_minor: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))
    dP_stage_total: Q_ = field(default_factory=lambda: Q_(0.0, "Pa"))

    hot_flow_A: Q_ = field(default_factory=lambda: Q_(0.0, "m^2"))
    cold_flow_A: Q_ = field(default_factory=lambda: Q_(0.0, "m^2"))
    hot_Dh: Q_ = field(default_factory=lambda: Q_(0.0, "m"))
    cold_Dh: Q_ = field(default_factory=lambda: Q_(0.0, "m"))

@dataclass(frozen=True)
class GlobalProfile:
    # flattened along exchanger x
    x: List[Q_]
    dx: List[Q_]
    gas: List[object]
    water: List[object]
    qprime: List[Q_]
    UA_prime: List[Q_]
    h_g: List[Q_]
    h_c: List[Q_]
    stage_index: List[int]
    stage_name: List[str]
    # ΔP per-step
    dP_fric: List[Q_]
    dP_minor: List[Q_]
    dP_total: List[Q_]
    # keep full stage_results for summary
    stage_results: List[StageResult]

def build_global_profile(stage_results: Sequence[StageResult]) -> GlobalProfile:
    xs: List[Q_] = []
    dxs: List[Q_] = []
    gas: List[object] = []
    water: List[object] = []
    qprime: List[Q_] = []
    UA_prime: List[Q_] = []
    h_g: List[Q_] = []
    h_c: List[Q_] = []
    sidx: List[int] = []
    sname: List[str] = []
    dP_fric: List[Q_] = []
    dP_minor: List[Q_] = []
    dP_total: List[Q_] = []

    for k, sr in enumerate(stage_results):
        for st in sr.steps:
            xs.append(st.x)
            dxs.append(st.dx)
            gas.append(st.gas)
            water.append(st.water)
            qprime.append(st.qprime)
            UA_prime.append(st.UA_prime)
            h_g.append(st.h_g)
            h_c.append(st.h_c)
            sidx.append(k if st.stage_index < 0 else st.stage_index)
            sname.append(sr.stage_name if not st.stage_name else st.stage_name)
            dP_fric.append(st.dP_fric)
            dP_minor.append(st.dP_minor)
            dP_total.append(st.dP_total)

    return GlobalProfile(
        x=xs, dx=dxs, gas=gas, water=water,
        qprime=qprime, UA_prime=UA_prime, h_g=h_g, h_c=h_c,
        stage_index=sidx, stage_name=sname,
        dP_fric=dP_fric, dP_minor=dP_minor, dP_total=dP_total,
        stage_results=list(stage_results),
    )

def write_results_csvs(
    global_profile: GlobalProfile,
    combustion: CombustionResult | None,
    outdir: str | Path,
    run_id: str,
) -> Tuple[str, str]:
    """
    Write step-level and stage summary CSVs.

    Returns:
        (steps_csv_path, summary_csv_path) as strings.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    steps_path = outdir / f"{run_id}_steps.csv"
    summary_path = outdir / f"{run_id}_summary.csv"

    df_steps = profile_to_dataframe(global_profile, remap_water=True)
    rows, _, _ = summary_from_profile(global_profile, combustion=combustion)

    df_steps.to_csv(steps_path, index=False)

    pd.DataFrame(
        rows,
        columns=[
            "stage_index", "stage_name", "stage_kind",
            "Q_stage[MW]", "UA_stage[MW/K]",
            "gas_in_T[°C]", "gas_out_T[°C]",
            "water_in_h[kJ/kg]", "water_out_h[kJ/kg]",
            "ΔP_stage_fric[Pa]", "ΔP_stage_minor[Pa]", "ΔP_stage_total[Pa]",
            "Q_conv_stage[MW]", "Q_rad_stage[MW]",
            "η_direct[-]", "η_indirect[-]",
            "Q_total_useful[MW]", "Q_in_total[MW]",
            "P_LHV[MW]", "LHV_mass[kJ/kg]",
        ],
    ).to_csv(summary_path, index=False)

    return str(steps_path), str(summary_path)

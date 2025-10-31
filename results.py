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


# ----------------------------
# Reporting / Visualization
# ----------------------------

class HeatExchangerReport:
    """Minimal analyzer for a list[StageResult].

    Outputs:
      - per-stage summary table
      - two quick plots: axial profiles and per-stage duty

    Assumptions:
      - All temperature quantities convertible to kelvin
      - Enthalpy convertible to J/kg
      - Axial coordinate in meters; if per-stage x resets, the code will
        build a monotone sequence using cumulative offsets.
    """

    def __init__(self, stages: List[StageResult]):
        if not stages:
            raise ValueError("stages list is empty")
        self.stages = stages

        # Build long-form step dataframe on init for convenience
        self._steps_df = self._build_steps_dataframe()
        self._stage_df = self._build_stage_dataframe()

    # ---------- public API ----------
    @property
    def steps(self) -> pd.DataFrame:
        """Long-form per-step dataframe.

        Columns:
          stage_name, stage_kind, step_index, x_m, gas_T_K, water_h_MJkg,
          Tgw_K, Tww_K, UA_prime_W_per_K_per_m, qprime_W_per_m, boiling
        """
        return self._steps_df.copy()

    @property
    def stages_summary(self) -> pd.DataFrame:
        """Per-stage compact summary table.

        Columns:
          stage_name, stage_kind, x_start_m, x_end_m,
          gas_T_in_K, gas_T_out_K, water_h_in_MJkg, water_h_out_MJkg,
          Q_stage_MW, UA_stage_kW_per_K, n_steps, boiling_any
        """
        return self._stage_df.copy()

    def plot_axial_profiles(self, *, show: bool = True, savepath: Optional[str] = None) -> plt.Figure:
        """Plot axial profiles: gas T, wall Ts (left axis), water h (right axis)."""
        df = self._steps_df
        fig, ax1 = plt.subplots(figsize=(9, 5))

        # Left y: temperatures
        ax1.plot(df["x_m"], df["gas_T_K"], color="tab:red", label="Gas T [K]")
        ax1.plot(df["x_m"], df["Tgw_K"], label="Gas wall Tgw [K]", linestyle=":")
        ax1.plot(df["x_m"], df["Tww_K"], label="Water wall Tww [K]", linestyle="--")
        ax1.set_xlabel("Axial position x [m]")
        ax1.set_ylabel("Temperature [K]")

        # Right y: water enthalpy
        ax2 = ax1.twinx()
        ax2.plot(df["x_m"], df["water_h_MJkg"], color="tab:blue", label="Water h [MJ/kg]")
        ax2.set_ylabel("Water enthalpy [MJ/kg]")

        # Stage separators
        for _, row in self._stage_df.iterrows():
            ax1.axvline(row["x_start_m"], alpha=0.15)
        ax1.axvline(self._stage_df["x_end_m"].iloc[-1], alpha=0.15)

        # Legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")
        fig.tight_layout()
        if savepath:
            fig.savefig(savepath, dpi=150)
        if show:
            plt.show()
        return fig

    def plot_stage_duty(self, *, show: bool = True, savepath: Optional[str] = None) -> plt.Figure:
        """Bar chart of Q_stage by stage."""
        df = self._stage_df
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(df["stage_name"], df["Q_stage_MW"])
        ax.set_ylabel("Stage duty Q [MW]")
        ax.set_xlabel("Stage")
        ax.set_xticklabels(df["stage_name"], rotation=45, ha="right")
        fig.tight_layout()
        if savepath:
            fig.savefig(savepath, dpi=150)
        if show:
            plt.show()
        return fig

    # ---------- internals ----------
    def _build_steps_dataframe(self) -> pd.DataFrame:
        rows: List[Dict] = []

        # Detect whether x is global or resets per stage. If it resets, accumulate offsets.
        stage_offsets: Dict[str, float] = {}
        offset = 0.0
        prev_last_x = None

        for s in self.stages:
            xs = [self._to_float(step.x, "meter") for step in s.steps if step.x is not None]
            if xs:
                first_x = xs[0]
                last_x = xs[-1]
            else:
                first_x = last_x = None

            # If this stage's first x is <= previous last x, start a new offset
            if prev_last_x is not None and first_x is not None and first_x <= prev_last_x:
                offset = prev_last_x  # continue from where previous left off
            stage_offsets[s.stage_name] = offset
            if last_x is not None:
                prev_last_x = max(prev_last_x or -math.inf, last_x + offset)

            for step in s.steps:
                x_m = self._to_float(step.x, "meter")
                if x_m is not None:
                    x_m += stage_offsets[s.stage_name]
                rows.append({
                    "stage_name": s.stage_name,
                    "stage_kind": s.stage_kind,
                    "step_index": step.i,
                    "x_m": x_m,
                    "gas_T_K": self._to_float(step.gas.T, "kelvin"),
                    "water_h_MJkg": self._to_float(step.water.h, "J/kg", scale=1e-6),
                    "Tgw_K": self._to_float(step.Tgw, "kelvin"),
                    "Tww_K": self._to_float(step.Tww, "kelvin"),
                    "UA_prime_W_per_K_per_m": self._to_float(step.UA_prime, "W/(K*m)"),
                    "qprime_W_per_m": self._to_float(step.qprime, "W/m"),
                    "boiling": bool(step.boiling),
                })

        df = pd.DataFrame(rows)
        # Sort by x then by an implicit ordering when x is None
        df = df.sort_values([col for col in ["x_m", "stage_name", "step_index"] if col in df.columns]).reset_index(drop=True)
        return df

    def _build_stage_dataframe(self) -> pd.DataFrame:
        rows: List[Dict] = []
        for s in self.stages:
            df_s = self._steps_df[self._steps_df["stage_name"] == s.stage_name]
            x_start = df_s["x_m"].min()
            x_end = df_s["x_m"].max()
            gas_T_in = df_s["gas_T_K"].iloc[0]
            gas_T_out = df_s["gas_T_K"].iloc[-1]
            water_h_in = df_s["water_h_MJkg"].iloc[0]
            water_h_out = df_s["water_h_MJkg"].iloc[-1]
            rows.append({
                "stage_name": s.stage_name,
                "stage_kind": s.stage_kind,
                "x_start_m": x_start,
                "x_end_m": x_end,
                "gas_T_in_K": gas_T_in,
                "gas_T_out_K": gas_T_out,
                "water_h_in_MJkg": water_h_in,
                "water_h_out_MJkg": water_h_out,
                "Q_stage_MW": self._to_float(s.Q_stage, "W", scale=1e-6),
                "UA_stage_kW_per_K": self._to_float(s.UA_stage, "W/K", scale=1e-3),
                "n_steps": len(s.steps),
                "boiling_any": bool(df_s["boiling"].any()),
            })
        return pd.DataFrame(rows)

    # ---------- utilities ----------
    @staticmethod
    def _to_float(q: Optional[Q_], unit: str, *, scale: float = 1.0) -> Optional[float]:
        if q is None:
            return None
        try:
            return float(q.to(unit).magnitude) / scale
        except Exception:
            # If already unitless but scale requested
            try:
                return float(q.magnitude) / scale  # type: ignore[attr-defined]
            except Exception:
                return float(q) / scale  # best effort

    # ---------- quick checks ----------
    def check_counterflow_monotonic(self) -> Dict[str, bool]:
        """Simple trend checks useful during commissioning."""
        df = self._steps_df.dropna(subset=["x_m"]).sort_values("x_m")
        gas_decreasing = bool(np.all(np.diff(df["gas_T_K"]) <= 1e-9))
        water_increasing = bool(np.all(np.diff(df["water_h_MJkg"]) >= -1e-12))
        return {"gas_T_decreases": gas_decreasing, "water_h_increases": water_increasing}


# ----------------------------
# Example usage (to be adapted in caller tests)
# ----------------------------
#
# report = HeatExchangerReport(stages=stage_results)
# df_steps = report.steps
# df_stages = report.stages_summary
# report.plot_axial_profiles()
# report.plot_stage_duty()
# print(report.check_counterflow_monotonic())

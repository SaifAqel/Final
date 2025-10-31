# main.py
from __future__ import annotations

import argparse
import sys
import logging
from math import ceil
from typing import List, Optional, Tuple

from units import Q_
from io_loader import load_config
from geometry import GeometryBuilder
from solver import solve_exchanger, solve_stage
from logging_utils import setup_logging
from postproc import results_to_dataframe
from models import GasStream, WaterStream, HXStage
from props import GasProps


# ------------------------- helpers -------------------------

def q_or_none(s: Optional[str]) -> Optional[Q_]:
    """Parse a pint quantity from a string. Return None on empty."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    return Q_(s)


def _compute_single_stage_steps(
    stage: HXStage,
    target_dx: Optional[Q_],
    min_steps: int,
    max_steps: int,
) -> int:
    L = stage.spec["inner_length"].to("m")
    dx_target = (L / 50).to("m") if target_dx is None else target_dx.to("m")
    n = int(ceil((L / dx_target).to("").magnitude))
    n = max(min_steps, min(max_steps, n))
    return n


def _build_summary(
    stage_results,
) -> Tuple[List[dict], float, float]:
    """
    Build list[dict] for per-stage summary and compute totals.
    Returns (rows, Q_total_W, UA_total_W_per_K)
    """
    rows: List[dict] = []
    Q_total = 0.0
    UA_total = 0.0
    for idx, st in enumerate(stage_results):
        s0 = st.steps[0]
        sN = st.steps[-1]
        row = {
            "stage_index": idx,
            "stage_name": st.stage_name,
            "stage_kind": st.stage_kind,
            "Q_stage[W]": st.Q_stage.to("W").magnitude,
            "UA_stage[W/K]": st.UA_stage.to("W/K").magnitude,
            "gas_in_T[K]": s0.gas.T.to("K").magnitude,
            "gas_out_T[K]": sN.gas.T.to("K").magnitude,
            "water_in_h[J/kg]": s0.water.h.to("J/kg").magnitude,
            "water_out_h[J/kg]": sN.water.h.to("J/kg").magnitude,
        }
        rows.append(row)
        Q_total += row["Q_stage[W]"]
        UA_total += row["UA_stage[W/K]"]
    return rows, Q_total, UA_total


def _print_text_summary(stage_results: List, is_six_stage: bool, gas_in: GasStream, gas_out: GasStream,
                        water_in: WaterStream, water_out: WaterStream) -> None:
    # One line per stage
    for idx, st in enumerate(stage_results):
        s0 = st.steps[0]
        sN = st.steps[-1]
        line = (
            f"{idx} {st.stage_name} {st.stage_kind} "
            f"Q_stage={st.Q_stage:~P} UA={st.UA_stage:~P} "
            f"gas_in={s0.gas.T:~P}→{sN.gas.T:~P} "
            f"water_in={s0.water.h:~P}→{sN.water.h:~P}"
        )
        print(line)

    if is_six_stage:
        # End-to-end check
        gprops = GasProps()
        h_g_in = gprops.h(gas_in.T, gas_in.P, gas_in.comp)
        h_g_out = gprops.h(gas_out.T, gas_out.P, gas_out.comp)
        Q_gas = (gas_in.mass_flow * (h_g_in - h_g_out)).to("W")

        h_w_in = water_in.h
        h_w_out = water_out.h
        Q_wat = (water_in.mass_flow * (h_w_out - h_w_in)).to("W")

        Q_total = sum((sr.Q_stage for sr in stage_results), Q_(0.0, "W")).to("W")
        denom = max(abs(Q_total).to("W").magnitude, 1e-12)
        rel = abs(Q_gas.to("W").magnitude - Q_wat.to("W").magnitude) / denom
        print(f"Q_total={Q_total:~P}")
        print(f"gas: m_dot*(h_in-h_out)={Q_gas:~P}")
        print(f"water: m_dot*(h_out-h_in)={Q_wat:~P}")
        print(f"rel_mismatch={rel*100:.3f}%")

# --------------------------- CLI ---------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run HX solver.")
    p.add_argument("--stages", default="config/stages.yaml")
    p.add_argument("--streams", default="config/streams.yaml")
    p.add_argument("--drum", default="config/drum.yaml")

    p.add_argument("--target-dx", default="", help="Pint string like '0.1 m'. Empty means auto.")
    p.add_argument("--min-steps", type=int, default=20)
    p.add_argument("--max-steps", type=int, default=400)
    p.add_argument("--max-passes", type=int, default=20)
    p.add_argument("--tol-Q", default="1e-3 W")
    p.add_argument("--tol-end", default="1e-3 J/kg")

    p.add_argument("--log-level",
                   choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR"],
                   default="INFO")

    p.add_argument("--csv", default="results.csv",
                   help="Path to write per-step results CSV.")
    p.add_argument("--summary", default="summary.csv",
                   help="Path to write per-stage summary CSV.")
    p.add_argument("--print-summary", action="store_true",
                   help="Print concise text summary to stdout.")

    return p.parse_args(argv)


# -------------------------- main ---------------------------

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # Setup logging early
    setup_logging(args.log_level)

    # Units from CLI
    target_dx = q_or_none(args.target_dx)
    tol_Q = Q_(args.tol_Q)
    tol_end = Q_(args.tol_end)

    try:
        stages_raw, gas, water, drum = load_config(args.stages, args.streams, args.drum)
    except Exception as e:
        logging.exception(f"Failed to load configs: {e}")
        return 1

    # Validate required inputs
    if gas is None or water is None or drum is None:
        print("error: missing required inputs. 'gas', 'water', and 'drum' must be provided.", file=sys.stderr)
        return 2

    # Geometry enrichment
    stages: List[HXStage] = GeometryBuilder(drum).enrich(stages_raw)

    try:
        # Solver selection
        if len(stages) == 6:
            stage_results, gas_out, water_out = solve_exchanger(
                stages,
                gas,
                water,
                target_dx=target_dx,
                min_steps_per_stage=args.min_steps,
                max_steps_per_stage=args.max_steps,
                max_passes=args.max_passes,
                tol_Q=tol_Q,
                tol_end=tol_end,
                log_level=args.log_level,
            )
        elif len(stages) == 1:
            n_steps = _compute_single_stage_steps(stages[0], target_dx, args.min_steps, args.max_steps)
            g_out, w_out, st_res = solve_stage(gas, water, stages[0], n_steps, stage_index=0)
            stage_results = [st_res]
            gas_out, water_out = g_out, w_out
        else:
            print("error: solver expects exactly 6 stages, or 1 stage for the fallback single-stage march.", file=sys.stderr)
            return 3

        # Post-processing: per-step CSV
        df = results_to_dataframe(stage_results)
        df.to_csv(args.csv, index=False)

        # Build per-stage summary and write CSV
        rows, _, _ = _build_summary(stage_results)
        import pandas as pd  # local import to avoid unnecessary dependency if not saving
        df_sum = pd.DataFrame(rows, columns=[
            "stage_index", "stage_name", "stage_kind",
            "Q_stage[W]", "UA_stage[W/K]",
            "gas_in_T[K]", "gas_out_T[K]",
            "water_in_h[J/kg]", "water_out_h[J/kg]",
        ])
        df_sum.to_csv(args.summary, index=False)

        # Optional stdout summary
        if args.print_summary:
            _print_text_summary(
                stage_results,
                is_six_stage=(len(stages) == 6),
                gas_in=gas, gas_out=gas_out,
                water_in=water, water_out=water_out
            )

        return 0

    except Exception as e:
        logging.exception(f"Solver failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())


# plots.py
# Minimal plotting utilities for the HX solver outputs.
# Reads results.csv and summary.csv, and generates:
# 1) Bar chart of Q_stage by stage.
# 2) Per-stage line plots of gas_T and water_T versus x.
#
# Usage examples:
#   python plots.py --summary summary.csv --results results.csv --out figs
#   python plots.py --results results.csv --stage 0 --out figs  # single stage temp plot
#
# Notes:
# - Uses matplotlib only. No seaborn.
# - Each chart is a single figure. No subplots.
# - No explicit colors or styles are set.

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def bar_Q_by_stage(summary_csv: str, out_dir: str) -> Path:
    df = pd.read_csv(summary_csv)
    fig, ax = plt.subplots()
    ax.bar(df["stage_name"], df["Q_stage[W]"])
    ax.set_ylabel("Q_stage [W]")
    ax.set_xlabel("Stage")
    ax.set_title("Stage heat duties")
    fig.tight_layout()
    out_path = Path(out_dir) / "Q_stage_bar.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path

def temps_along_x_for_stage(results_csv: str, stage_index: int, out_dir: str) -> Path:
    df = pd.read_csv(results_csv)
    sdf = df[df["stage_index"] == stage_index].copy()
    if sdf.empty:
        raise ValueError(f"No rows for stage_index={stage_index}")
    # Sort by x for clean lines
    sdf = sdf.sort_values(by="x[m]")

    fig, ax = plt.subplots()
    ax.plot(sdf["x[m]"], sdf["gas_T[K]"], label="gas_T [K]")
    ax.plot(sdf["x[m]"], sdf["water_T[K]"], label="water_T [K]")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Temperature [K]")
    title = f"Stage {stage_index}: gas_T and water_T vs x"
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    out_path = Path(out_dir) / f"temps_stage_{stage_index}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path

def make_all(results_csv: str, summary_csv: str, out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    # Bar chart for Q by stage
    bar_Q_by_stage(summary_csv, out)
    # Temperature plots for every stage present
    df = pd.read_csv(results_csv)
    stage_ids = sorted(df["stage_index"].dropna().unique().astype(int).tolist())
    for s in stage_ids:
        temps_along_x_for_stage(results_csv, s, out)

def parse_args():
    p = argparse.ArgumentParser(description="Plot HX results")
    p.add_argument("--results", default="results.csv")
    p.add_argument("--summary", default="summary.csv")
    p.add_argument("--out", default="figs")
    p.add_argument("--stage", type=int, default=None, help="If set, plot only this stage's temperature profile")
    return p.parse_args()

def main():
    args = parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    if args.stage is None:
        make_all(args.results, args.summary, args.out)
    else:
        temps_along_x_for_stage(args.results, args.stage, args.out)
    # Always create bar chart when summary exists
    if Path(args.summary).exists():
        bar_Q_by_stage(args.summary, args.out)

if __name__ == "__main__":
    main()

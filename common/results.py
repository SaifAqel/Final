from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Sequence, Tuple
from common.units import Q_
from common.models import GasStream, WaterStream
from pathlib import Path
import pandas as pd


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
) -> Tuple[str, str, str]:
    """
    Write CSVs:
      - <run_id>_steps.csv           : step-level data (unchanged).
      - <run_id>_stages_summary.csv  : per-stage summary (no TOTAL_BOILER row).
      - <run_id>_boiler_summary.csv  : single-row boiler summary, based on
                                       TOTAL_BOILER with extra boiler-level data.

    Returns:
        (steps_csv_path, stages_summary_csv_path) as strings.
    """
    from heat.postproc import profile_to_dataframe, summary_from_profile

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Paths
    steps_path = outdir / f"{run_id}_steps.csv"
    stages_summary_path = outdir / f"{run_id}_stages_summary.csv"
    boiler_summary_path = outdir / f"{run_id}_boiler_summary.csv"

    # 2) Build step dataframe (unchanged)
    df_steps = profile_to_dataframe(global_profile, remap_water=True)
    df_steps.to_csv(steps_path, index=False)

    # 3) Get summary rows (unchanged structure: stages + TOTAL_BOILER at the end)
    rows, _, _ = summary_from_profile(global_profile, combustion=combustion)

    df_summary = pd.DataFrame(
        rows,
        columns=[
            "stage_index", "stage_name", "stage_kind",
            "Q_stage[MW]", "UA_stage[MW/K]",

            "gas_in_P[Pa]", "gas_in_T[°C]", "gas_in_h[kJ/kg]",
            "gas_out_P[Pa]", "gas_out_T[°C]", "gas_out_h[kJ/kg]",

            "water_in_P[Pa]", "water_in_T[°C]", "water_in_h[kJ/kg]",
            "water_out_P[Pa]", "water_out_T[°C]", "water_out_h[kJ/kg]",

            "ΔP_stage_fric[Pa]", "ΔP_stage_minor[Pa]", "ΔP_stage_total[Pa]",
            "Q_conv_stage[MW]", "Q_rad_stage[MW]",
            "η_direct[-]", "η_indirect[-]",
            "Q_total_useful[MW]", "Q_in_total[MW]",
            "P_LHV[MW]", "LHV_mass[kJ/kg]",
        ],
    )

    # 4) Split into:
    #    - pure per-stage table (no TOTAL_BOILER)
    #    - boiler row (TOTAL_BOILER), later turned into boiler summary
    df_stages = df_summary[df_summary["stage_name"] != "TOTAL_BOILER"].copy()
    df_boiler = df_summary[df_summary["stage_name"] == "TOTAL_BOILER"].copy()

    # 4a) Add steam capacity per stage (same for all stages, from water[0])
    try:
        m_w = global_profile.water[0].mass_flow.to("kg/s").magnitude
    except Exception:
        m_w = float("nan")
    df_stages["steam_capacity[t/h]"] = m_w * 3600.0 / 1000.0  # kg/s → t/h

    # 4b) Reindex by stage name so columns = stages
    df_stages = df_stages.set_index("stage_name")

    # 4c) Build the requested row layout (rows = quantities, columns = stages)
    table = pd.DataFrame(
        {
            # identifiers
            "name": df_stages.index,
            "kind": df_stages["stage_kind"],

            # gas side (Pa, °C, kJ/kg)
            "gas in pressure": df_stages["gas_in_P[Pa]"],
            "gas in temp": df_stages["gas_in_T[°C]"],
            "gas in enthalpy": df_stages["gas_in_h[kJ/kg]"],
            "gas out pressure": df_stages["gas_out_P[Pa]"],
            "gas out temp": df_stages["gas_out_T[°C]"],
            "gas out enthalpy": df_stages["gas_out_h[kJ/kg]"],

            # water side (water pressure taken at inlet)
            "water pressure": df_stages["water_in_P[Pa]"],
            "water in temp": df_stages["water_in_T[°C]"],
            "water in enthalpy": df_stages["water_in_h[kJ/kg]"],
            "water out temp": df_stages["water_out_T[°C]"],
            "water out enthalpy": df_stages["water_out_h[kJ/kg]"],

            # pressure drops
            "pressure drop fric": df_stages["ΔP_stage_fric[Pa]"],
            "pressure drop minor": df_stages["ΔP_stage_minor[Pa]"],
            "pressure drop total": df_stages["ΔP_stage_total[Pa]"],

            # heat duties and UA (MW / MW/K)
            "Q conv": df_stages["Q_conv_stage[MW]"],
            "Q rad": df_stages["Q_rad_stage[MW]"],
            "Q total": df_stages["Q_stage[MW]"],
            "UA": df_stages["UA_stage[MW/K]"],

            # steam capacity (t/h)
            "steam capacity": df_stages["steam_capacity[t/h]"],
        },
        index=df_stages.index,
    )

    # rows = requested headers, columns = stage names
    table.T.to_csv(stages_summary_path)




    # 4b) Build boiler summary in "column" fashion, with added boiler-level data
    if not df_boiler.empty:
        # Start from the TOTAL_BOILER row
        boiler_row = df_boiler.iloc[0].copy()

        # Add water mass flow (taken as constant from first water snapshot)
        try:
            m_w = global_profile.water[0].mass_flow.to("kg/s").magnitude
        except Exception:
            m_w = float("nan")
        boiler_row["water_mass_flow[kg/s]"] = m_w
        boiler_row["steam_capacity[t/h]"] = m_w * 3600.0 / 1000.0


        # Add combustion-related scalars if available
        if combustion is not None:
            # Adiabatic flame temperature
            try:
                boiler_row["T_ad[°C]"] = combustion.T_ad.to("degC").magnitude
            except Exception:
                boiler_row["T_ad[°C]"] = ""

            # Mass-based LHV of fuel (if provided)
            if combustion.fuel_LHV_mass is not None:
                boiler_row["fuel_LHV_mass[kJ/kg]"] = combustion.fuel_LHV_mass.to("kJ/kg").magnitude
            else:
                boiler_row["fuel_LHV_mass[kJ/kg]"] = ""

            # Firing capacity based on LHV (if provided)
            if combustion.fuel_P_LHV is not None:
                boiler_row["fuel_P_LHV[MW]"] = combustion.fuel_P_LHV.to("MW").magnitude
            else:
                boiler_row["fuel_P_LHV[MW]"] = ""
        else:
            boiler_row["T_ad[°C]"] = ""
            boiler_row["fuel_LHV_mass[kJ/kg]"] = ""
            boiler_row["fuel_P_LHV[MW]"] = ""

        # Single-row boiler summary CSV

        # Single-row boiler summary CSV

        boiler_df = pd.DataFrame([boiler_row])

        # Build reordered / reduced boiler summary in "row" fashion
        # New row labels (index) and corresponding source columns:
        summary_mapping = [
            ("LHV",                         "LHV_mass[kJ/kg]"),
            ("P-LHV",                       "P_LHV[MW]"),
            ("Q_in total",                  "Q_in_total[MW]"),
            ("Tad",                         "T_ad[°C]"),
            ("UA",                          "UA_stage[MW/K]"),
            ("Quseful",                     "Q_total_useful[MW]"),
            ("pressue drop firc total",     "ΔP_stage_fric[Pa]"),
            ("pressure drop minor total",   "ΔP_stage_minor[Pa]"),
            ("pressure drop total",         "ΔP_stage_total[Pa]"),
            ("eta direct",                  "η_direct[-]"),
            ("eta indirect",                "η_indirect[-]"),
            ("water flow",                  "water_mass_flow[kg/s]"),
            ("steam capacity",              "steam_capacity[t/h]"),
            ("stack temperature",           "gas_out_T[°C]"),
        ]

        out_row = {}
        row0 = boiler_df.iloc[0]
        for new_name, src_col in summary_mapping:
            # Keep units as they are in the existing columns; if a source
            # column is missing, leave the value blank.
            out_row[new_name] = row0.get(src_col, "")

        boiler_out_df = pd.DataFrame([out_row])
        # Transpose so the CSV has your requested labels as row headers
        boiler_out_df.T.to_csv(boiler_summary_path)



    else:
        # No TOTAL_BOILER row – write an empty boiler summary with the same columns
        empty_cols = list(df_summary.columns) + [
            "water_mass_flow[kg/s]",
            "T_ad[°C]",
            "fuel_LHV_mass[kJ/kg]",
            "fuel_P_LHV[MW]",
        ]
        pd.DataFrame(columns=empty_cols).to_csv(boiler_summary_path, index=False)

    # Keep function return type: second return value now points to stages_summary
    return str(steps_path), str(stages_summary_path), str(boiler_summary_path)   


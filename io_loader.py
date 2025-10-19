from typing import Tuple, List, Dict, Any
import yaml
from logging_utils import trace_calls
from units import Q_, ureg
from models import HXStage, GasStream, WaterStream

def _get(d: Dict[str, Any], path: str, default=None):
    cur = d
    for k in path.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def _q(node: Any, default_unit: str | None = None) -> Q_:
    # Accept {value, unit} or bare number with default_unit or bare string "3 m"
    if isinstance(node, dict) and "value" in node and "unit" in node:
        return Q_(node["value"], str(node["unit"]))
    if isinstance(node, (int, float)) and default_unit:
        return Q_(node, default_unit)
    if isinstance(node, str):
        # e.g. "3 m" or "150 W/(m*K)"
        return Q_(*node.split(maxsplit=1)) if " " in node else Q_(float(node), default_unit or "")
    raise ValueError(f"Cannot parse quantity from: {node!r}")

@trace_calls()
def load_config(stages_path: str = "config/stages.yaml",
                streams_path: str | None = "config/streams.yaml"
               ) -> Tuple[List[HXStage], GasStream, WaterStream]:
    with open(stages_path, "r", encoding="utf-8") as fh:
        sdoc = yaml.safe_load(fh)

    stages: List[HXStage] = []
    for name, node in sdoc["stages"].items():
        L = (_q(_get(node, "hot_side.inner_length")) if _get(node, "hot_side.inner_length")
             else _q(_get(node, "cold_side.inner_length")) if _get(node, "cold_side.inner_length")
             else Q_(1.0, "m"))

        # Unified stage specification; backward-compat mapping.
        spec: Dict[str, Q_] = {}

        # UA baseline (overridable later)
        spec["UA_per_m"] = Q_(150.0, "W/K/m")

        # Hot-side mapped fields
        if _get(node, "hot_side.inner_diameter"):
            spec["hot_Di"] = _q(_get(node, "hot_side.inner_diameter"))
        if _get(node, "hot_side.curvature_radius"):
            spec["hot_Rc"] = _q(_get(node, "hot_side.curvature_radius"))
        if _get(node, "hot_side.tubes_number"):
            spec["hot_ntubes"] = _q(_get(node, "hot_side.tubes_number"), "dimensionless")
        if _get(node, "hot_side.pitch"):
            spec["hot_pitch"] = _q(_get(node, "hot_side.pitch"))

        # Cold-side mapped fields
        if _get(node, "cold_side.inner_diameter"):
            spec["cold_Di"] = _q(_get(node, "cold_side.inner_diameter"))
        if _get(node, "cold_side.flow_area"):
            spec["cold_Ai"] = _q(_get(node, "cold_side.flow_area"))
        if _get(node, "cold_side.wetted_perimeter"):
            spec["cold_Pi"] = _q(_get(node, "cold_side.wetted_perimeter"))

        # Optional direct UA override at stage root
        if _get(node, "UA_per_m"):
            spec["UA_per_m"] = _q(_get(node, "UA_per_m"))

        stages.append(HXStage(name=name, kind="generic", L=L, spec=spec))

    # streams
    if streams_path:
        with open(streams_path, "r", encoding="utf-8") as fh:
            tdoc = yaml.safe_load(fh)
        g = tdoc["gas_stream"]
        w = tdoc["water_stream"]
        gas = GasStream(
            mass_flow=_q(g["mass_flow_rate"]),
            T=_q(g["temperature"]),
            P=_q(g["pressure"]),
            comp={k: _q(v, "dimensionless") for k, v in g.get("composition", {}).items()},
        )
        water = WaterStream(
            mass_flow=_q(w["mass_flow_rate"]),
            h=_q(w["enthalpy"]),
            P=_q(w["pressure"]),
        )
    else:
        gas = GasStream(Q_(9.427,"kg/s"), Q_(2176.59,"K"), Q_(115000,"Pa"))
        water = WaterStream(Q_(3.333,"kg/s"), Q_(439000,"J/kg"), Q_(1_000_000,"Pa"))

    return stages, gas, water

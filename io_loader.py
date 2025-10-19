from typing import Tuple, List, Dict, Any
import yaml
from units import Q_
from models import HXStage, GasStream, WaterStream, Drum


def _q(node: Any) -> Q_:
    if isinstance(node, dict) and "value" in node and "unit" in node:
        unit = str(node["unit"])
        if unit == "micrometer":
            unit = "um"
        if unit in ("dimensionless", "1"):
            unit = ""
        return Q_(node["value"], unit)
    raise ValueError(f"Invalid quantity format: {node!r}")


def _get(d: Dict[str, Any], key: str, default=None):
    return d.get(key, default) if isinstance(d, dict) else default


def _wall_to_spec(prefix: str, wall: Dict[str, Any], spec: Dict[str, Q_]):
    if not wall:
        return

    if _get(wall, "thickness"):
        spec[f"{prefix}_wall_t"] = _q(_get(wall, "thickness"))
    if _get(wall, "conductivity"):
        spec[f"{prefix}_wall_k"] = _q(_get(wall, "conductivity"))

    surf_in = _get(_get(wall, "surfaces"), "inner") or {}
    if _get(surf_in, "roughness"):
        spec[f"{prefix}_roughness_in"] = _q(_get(surf_in, "roughness"))
    if _get(surf_in, "emissivity"):
        spec[f"{prefix}_eps_in"] = _q(_get(surf_in, "emissivity"))
    if _get(surf_in, "fouling_thickness"):
        spec[f"{prefix}_foul_t_in"] = _q(_get(surf_in, "fouling_thickness"))
    if _get(surf_in, "fouling_conductivity"):
        spec[f"{prefix}_foul_k_in"] = _q(_get(surf_in, "fouling_conductivity"))

    surf_out = _get(_get(wall, "surfaces"), "outer") or {}
    if _get(surf_out, "roughness"):
        spec[f"{prefix}_roughness_out"] = _q(_get(surf_out, "roughness"))
    if _get(surf_out, "emissivity"):
        spec[f"{prefix}_eps_out"] = _q(_get(surf_out, "emissivity"))
    if _get(surf_out, "fouling_thickness"):
        spec[f"{prefix}_foul_t_out"] = _q(_get(surf_out, "fouling_thickness"))
    if _get(surf_out, "fouling_conductivity"):
        spec[f"{prefix}_foul_k_out"] = _q(_get(surf_out, "fouling_conductivity"))


def _map_nozzles(prefix: str, side: Dict[str, Any], spec: Dict[str, Q_]):
    nozzles = _get(side, "nozzles") or {}
    inlet = _get(nozzles, "inlet") or {}
    outlet = _get(nozzles, "outlet") or {}

    if _get(inlet, "k"):
        spec[f"{prefix}_nozzle_k_in"] = _q(_get(inlet, "k"))
    if _get(outlet, "k"):
        spec[f"{prefix}_nozzle_k_out"] = _q(_get(outlet, "k"))


def load_config(stages_path: str, streams_path: str | None = None, drum_path: str | None = None):
    with open(stages_path, "r", encoding="utf-8") as fh:
        sdoc = yaml.safe_load(fh)

    stages: List[HXStage] = []
    for name, node in sdoc["stages"].items():
        hot = _get(node, "hot_side") or {}
        spec: Dict[str, Q_] = {}

        if _get(hot, "inner_diameter"):
            spec["hot_Di"] = _q(_get(hot, "inner_diameter"))
        if _get(hot, "inner_length"):
            spec["hot_L"] = _q(_get(hot, "inner_length"))
        if _get(hot, "curvature_radius"):
            spec["hot_Rc"] = _q(_get(hot, "curvature_radius"))
        if _get(hot, "tubes_number"):
            spec["hot_ntubes"] = _q(_get(hot, "tubes_number"))
        if _get(hot, "pitch"):
            spec["hot_pitch"] = _q(_get(hot, "pitch"))

        _wall_to_spec("hot", _get(hot, "wall"), spec)
        _map_nozzles("hot", hot, spec)

        L = spec.get("hot_L", Q_(1, "m"))
        stages.append(HXStage(name=name, kind="generic", L=L, spec=spec))

    drum = None
    if drum_path:
        ddoc = yaml.safe_load(open(drum_path, "r", encoding="utf-8"))
        w = ddoc["water_drum"]
        drum = Drum(
            Di=_q(w["inner_diameter"]).to("m"),
            L=_q(w["length"]).to("m"),
        )

    if not streams_path:
        return stages, None, None

    with open(streams_path, "r", encoding="utf-8") as fh:
        tdoc = yaml.safe_load(fh)

    g = tdoc["gas_stream"]
    w = tdoc["water_stream"]

    gas = GasStream(
        mass_flow=_q(g["mass_flow_rate"]),
        T=_q(g["temperature"]),
        P=_q(g["pressure"]),
        comp={k: _q(v) for k, v in _get(g, "composition", {}).items()},
    )
    water = WaterStream(
        mass_flow=_q(w["mass_flow_rate"]),
        h=_q(w["enthalpy"]),
        P=_q(w["pressure"]),
    )

    return stages, gas, water, drum

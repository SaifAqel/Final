# geometry.py
from math import pi
from dataclasses import replace
from typing import List
from units import Q_
from models import HXStage, Drum

class GeometryBuilder:
    def __init__(self, drum: Drum):
        self.drum = drum

    def enrich(self, stages: List[HXStage]) -> List[HXStage]:
        out: List[HXStage] = []
        for stg in stages:
            spec = dict(stg.spec)  # copy bag

            Di_hot = spec.get("hot_Di")
            t_hot  = spec.get("hot_wall_t")
            nt_q   = spec.get("hot_ntubes", Q_(1.0, "dimensionless"))

            if Di_hot is None or t_hot is None:
                raise KeyError(f"{stg.name}: hot_Di and hot_wall_t required")

            nt = int(nt_q.to("dimensionless").magnitude)
            Do = (Di_hot + 2*t_hot).to("m")

            A_drum = (pi * (self.drum.Di/2)**2) * Q_(1, "dimensionless")
            A_tubes = nt * (pi * (Do/2)**2) * Q_(1, "dimensionless")
            A_flow  = (A_drum - A_tubes).to("m^2")
            if A_flow.magnitude <= 0:
                raise ValueError(f"{stg.name}: nonpositive cold flow area")

            P_wet = (nt * pi * Do).to("m")    # interface perimeter per unit length
            Dh    = (4 * A_flow / P_wet).to("m")

            spec["cold_Ai"] = A_flow
            spec["cold_Pi"] = P_wet
            spec["cold_Di"] = Dh

            # propagate water-side surface data if present
            if self.drum.roughness: spec["cold_roughness_in"] = self.drum.roughness
            if self.drum.emissivity: spec["cold_eps_in"] = self.drum.emissivity
            if self.drum.foul_t: spec["cold_foul_t_in"] = self.drum.foul_t
            if self.drum.foul_k: spec["cold_foul_k_in"] = self.drum.foul_k

            out.append(replace(stg, spec=spec))
        return out

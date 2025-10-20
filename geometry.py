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
            spec = dict(stg.spec)

            if "inner_diameter" not in spec or "wall_t" not in spec:
                raise KeyError(f"{stg.name}: 'inner_diameter' and 'wall_t' required")

            Di_hot = spec["inner_diameter"].to("m")
            t_hot  = spec["wall_t"].to("m")
            nt_q   = spec.get("tubes_number", Q_(1.0, "dimensionless"))
            nt     = int(nt_q.to("dimensionless").magnitude)

            Do = (Di_hot + 2*t_hot).to("m")

            A_drum  = (pi * (self.drum.Di/2)**2) * Q_(1, "dimensionless")
            A_tubes = nt * (pi * (Do/2)**2) * Q_(1, "dimensionless")
            A_flow  = (A_drum - A_tubes).to("m^2")
            if A_flow.magnitude <= 0:
                raise ValueError(f"{stg.name}: nonpositive cold flow area")

            P_wet = (nt * pi * Do).to("m")
            Dh    = (4 * A_flow / P_wet).to("m")

            spec["cold_Ai"] = A_flow
            spec["cold_Pi"] = P_wet
            spec["cold_Di"] = Dh

            # propagate optional water-side surface data
            out.append(replace(stg, spec=spec))
        return out

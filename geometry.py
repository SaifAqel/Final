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
            if stg.kind == "single_tube":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                drum_inner_diameter = self.drum.Di.to("m")

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (pi * inner_diameter).to("m")
                outer_perimeter = (pi * outer_diameter).to("m")
                A_drum  = (pi * (drum_inner_diameter/2)**2)
                A_tube_out = (pi * (outer_diameter/2)**2)
                A_flow_cold  = (A_drum - A_tube_out).to("m^2")
                A_flow_hot = (pi * (inner_diameter/2)**2)
                P_wet = (pi * outer_diameter).to("m")
                Dh    = (4 * A_flow_cold / P_wet).to("m")

                spec["outer_diameter"] = outer_diameter
                spec["inner_perimeter"] = inner_perimeter
                spec["outer_perimeter"] = outer_perimeter
                spec["cold_flow_A"] = A_flow_cold
                spec["cold_wet_P"] = P_wet
                spec["cold_Dh"] = Dh
                spec["hot_flow_A"] = A_flow_hot

                out.append(replace(stg, spec=spec))

            elif stg.kind == "tube_bank":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                tubes_number   = spec["tubes_number"].to("")
                ST      = spec["ST"].to("m")
                SL      = spec["SL"].to("m")
                N_rows  = spec["N_rows"].to("")
                arrangement = spec["arrangement"]
                B = spec["baffle_spacing"]

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                n_cols = (tubes_number + N_rows - 1) // N_rows
                bundle_width = n_cols * ST
                A_fr = (B * bundle_width).to("m^2")
                P_wet = (tubes_number * pi * outer_diameter).to("m")

                spec["outer_diameter"]   = outer_diameter
                spec["inner_perimeter"]  = (tubes_number * pi * inner_diameter).to("m")
                spec["hot_flow_A"]       = (tubes_number * (pi * (inner_diameter/2)**2)).to("m^2")
                spec["cold_flow_A"] = A_fr
                spec["cold_wet_P"] = P_wet

                if arrangement == "inline":
                    spec["umax_factor"] = (ST / (ST - outer_diameter)).to("dimensionless")
                else:
                    spec["umax_factor"] = (ST / (ST - 0.5*outer_diameter)).to("dimensionless")

                out.append(replace(stg, spec=spec))

            elif stg.kind == "reversal_chamber":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                drum_inner_diameter = self.drum.Di.to("m")

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (pi * inner_diameter).to("m")
                A_drum  = (pi * (drum_inner_diameter/2)**2)
                A_tube_out = (pi * (outer_diameter/2)**2)
                A_flow  = (A_drum - A_tube_out).to("m^2")
                A_hot = (pi * (inner_diameter/2)**2)
                P_wet = (pi * outer_diameter).to("m")
                Dh    = (4 * A_flow / P_wet).to("m")

                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_flow_A"] = A_flow
                spec["cold_wet_P"] = P_wet
                spec["cold_Dh"] = Dh
                spec["hot_flow_A"] = A_hot

                out.append(replace(stg, spec=spec))

            elif stg.kind == "economiser":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                tubes_number   = spec["tubes_number"]
                drum_inner_diameter = self.drum.Di.to("m")

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (tubes_number * pi * inner_diameter).to("m")
                A_drum  = (pi * (drum_inner_diameter/2)**2)
                A_tube_out = tubes_number * (pi * (outer_diameter/2)**2)
                hot_flow_A  = (A_drum - A_tube_out).to("m^2")
                cold_flow_A = tubes_number * (pi * (inner_diameter/2)**2)
                P_wet = (tubes_number * pi * outer_diameter).to("m")
                Dh    = (4 * A_flow / P_wet).to("m")

                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_flow_A"] = cold_flow_A
                spec["cold_wet_P"] = (inner_perimeter * tubes_number)
                spec["hot_wet_P"] = P_wet
                spec["hot_Dh"] = Dh
                spec["hot_flow_A"] = hot_flow_A

                out.append(replace(stg, spec=spec))
            else:
                raise ValueError("kind of stage must be given")
        return out

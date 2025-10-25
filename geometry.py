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

                spec["roughness_cold_surface"] = spec["roughness_out"]
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

                nt = int(spec["tubes_number"].to("").magnitude)
                nrows = int(spec["N_rows"].to("").magnitude)
                n_cols = (nt + nrows - 1) // nrows


                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                bundle_width = Q_(n_cols, "") * ST

                Ds = self.drum.Di
                B  = spec["baffle_spacing"].to("m")
                pt = spec["ST"].to("m")
                do = (inner_diameter + 2*wall_t).to("m")

                FAR = (1 - do/pt).to("dimensionless")
                A_cross = (Ds * B * FAR).to("m^2")

                spec["cold_flow_A"] = A_cross
                spec["umax_factor"] = ((Ds * B) / A_cross).to("dimensionless")

                P_wet = (tubes_number * pi * outer_diameter).to("m")

                spec["roughness_cold_surface"] = spec["roughness_out"]
                spec["outer_diameter"]   = outer_diameter
                spec["inner_perimeter"]  = (tubes_number * pi * inner_diameter).to("m")
                spec["hot_flow_A"]       = (tubes_number * (pi * (inner_diameter/2)**2)).to("m^2")
                spec["cold_wet_P"] = P_wet

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

                spec["roughness_cold_surface"] = spec["roughness_out"]
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
                Dh    = (4 * hot_flow_A / P_wet).to("m")

                spec["roughness_cold_surface"] = spec["roughness_in"]
                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_flow_A"] = cold_flow_A
                spec["cold_wet_P"] = inner_perimeter
                spec["hot_wet_P"] = P_wet
                spec["hot_Dh"] = Dh
                spec["hot_flow_A"] = hot_flow_A

                out.append(replace(stg, spec=spec))
            else:
                raise ValueError("kind of stage must be given")
        return out

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
                tubes_number   = spec.get("tubes_number", Q_(1.0, "dimensionless"))

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (tubes_number * pi * inner_diameter).to("m")
                A_drum  = (pi * (self.drum.Di/2)**2)
                A_tubes = tubes_number * (pi * (outer_diameter/2)**2)
                A_flow  = (A_drum - A_tubes).to("m^2")

                A_hot = tubes_number * (pi * (inner_diameter/2)**2)

                P_wet = (tubes_number * pi * outer_diameter).to("m")
                Dh    = (4 * A_flow / P_wet).to("m")
                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_Ai"] = A_flow
                spec["cold_Pi"] = P_wet
                spec["cold_Di"] = Dh
                spec["hot_A"] = A_hot

                out.append(replace(stg, spec=spec))

            elif stg.kind == "tube_bank":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                tubes_number   = spec.get("tubes_number", Q_(1.0, "dimensionless"))
                ST      = spec["ST"].to("m")
                SL      = spec["SL"].to("m")
                N_rows  = spec["N_rows"].to("")
                arrangement = spec["arrangement"]
                B = spec["baffle_spacing"]

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                n_cols = (tubes_number + N_rows - 1) // N_rows  # ceil

                # bundle size and frontal crossflow area
                bundle_width = Q_(n_cols, "dimensionless") * ST          # normal to flow
                A_fr = (B * bundle_width).to("m^2")                      # approach/frontal flow area

                # wetted perimeter (outside of tubes, for fouling res etc.)
                P_wet = (tubes_number * pi * outer_diameter).to("m")

                # store
                spec["outer_diameter"]   = outer_diameter
                spec["inner_perimeter"]  = (tubes_number * pi * inner_diameter).to("m")
                spec["hot_A"]            = (tubes_number * (pi * (inner_diameter/2)**2)).to("m^2")

                # Use frontal area for shell-side hydraulics
                spec["cold_Ai"] = A_fr
                spec["cold_Pi"] = P_wet

                # Force downstream code to use Do in Re; set "cold_Di" = Do and provide Umax factor
                if arrangement == "inline":
                    spec["umax_factor"] = (ST / (ST - outer_diameter)).to("dimensionless")
                else:  # staggered
                    spec["umax_factor"] = (ST / (ST - 0.5*outer_diameter)).to("dimensionless")

                out.append(replace(stg, spec=spec))

            elif stg.kind == "reversal_chamber":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                tubes_number   = spec.get("tubes_number", Q_(1.0, "dimensionless"))

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (tubes_number * pi * inner_diameter).to("m")
                A_drum  = (pi * (self.drum.Di/2)**2)
                A_tubes = tubes_number * (pi * (outer_diameter/2)**2)
                A_flow  = (A_drum - A_tubes).to("m^2")

                A_hot = tubes_number * (pi * (inner_diameter/2)**2)

                P_wet = (tubes_number * pi * outer_diameter).to("m")
                Dh    = (4 * A_flow / P_wet).to("m")
                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_Ai"] = A_flow
                spec["cold_Pi"] = P_wet
                spec["cold_Di"] = Dh
                spec["hot_A"] = A_hot

                out.append(replace(stg, spec=spec))
            elif stg.kind == "economiser":
                spec = dict(stg.spec)

                inner_diameter = spec["inner_diameter"].to("m")
                wall_t  = spec["wall_t"].to("m")
                tubes_number   = spec.get("tubes_number", Q_(1.0, "dimensionless"))

                outer_diameter = (inner_diameter + 2*wall_t).to("m")
                inner_perimeter = (tubes_number * pi * inner_diameter).to("m")
                A_drum  = (pi * (self.drum.Di/2)**2)
                A_tubes = tubes_number * (pi * (outer_diameter/2)**2)
                A_flow  = (A_drum - A_tubes).to("m^2")

                A_hot = tubes_number * (pi * (inner_diameter/2)**2)

                P_wet = (tubes_number * pi * outer_diameter).to("m")
                Dh    = (4 * A_flow / P_wet).to("m")
                spec["inner_perimeter"] = inner_perimeter
                spec["outer_diameter"] = outer_diameter
                spec["cold_Ai"] = A_flow
                spec["cold_Pi"] = P_wet
                spec["cold_Di"] = Dh
                spec["hot_A"] = A_hot

                out.append(replace(stg, spec=spec))
            else:
                raise "kind of stage must be given"
        return out

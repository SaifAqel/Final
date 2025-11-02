from math import pi
from dataclasses import replace
from typing import List
from common.models import HXStage, Drum

class GeometryBuilder:
    def __init__(self, drum: Drum):
        self.drum = drum

    def enrich(self, stages: List[HXStage]) -> List[HXStage]:
        out: List[HXStage] = []
        for stg in stages:
            if stg.kind == "single_tube":
                spec = dict(stg.spec)
                Di_t = spec["inner_diameter"].to("m")
                t = spec["wall_t"].to("m")
                Do_t = (Di_t + 2*t).to("m")
                Di_drum = self.drum.Di.to("m")

                A_drum = (pi * (Di_drum/2)**2).to("m^2")
                A_tube_out = (pi * (Do_t/2)**2).to("m^2")

                spec["outer_diameter"] = Do_t
                spec["roughness_cold_surface"] = spec["roughness_out"]

                hot_wet_P = (pi * Di_t).to("m")
                hot_flow_A = (pi * (Di_t/2)**2).to("m^2")
                hot_Dh = (4 * hot_flow_A / hot_wet_P).to("m")

                cold_wet_P = (pi * Do_t).to("m")
                cold_flow_A = (A_drum - A_tube_out).to("m^2")
                cold_Dh = (4 * cold_flow_A / cold_wet_P).to("m")

                spec.update({
                    "hot_wet_P": hot_wet_P, "hot_flow_A": hot_flow_A, "hot_Dh": hot_Dh,
                    "cold_wet_P": cold_wet_P, "cold_flow_A": cold_flow_A, "cold_Dh": cold_Dh
                })
                out.append(replace(stg, spec=spec))

            elif stg.kind == "tube_bank":
                spec = dict(stg.spec)
                Di_t = spec["inner_diameter"].to("m")
                t = spec["wall_t"].to("m")
                Nt = spec["tubes_number"].to("")
                Do_t = (Di_t + 2*t).to("m")

                Ds = spec["shell_inner_diameter"].to("m")
                B  = spec["baffle_spacing"].to("m")
                pt = spec["ST"].to("m")

                FAR = (1 - (Do_t/pt)).to("dimensionless")
                A_gross = (Ds * B).to("m^2")
                A_cross = (A_gross * FAR).to("m^2")
                spec["umax_factor"] = (A_gross / A_cross).to("dimensionless")  # = 1/FAR
                spec["roughness_cold_surface"] = spec["roughness_out"]
                spec["outer_diameter"] = Do_t


                A_drum = (pi * (Ds/2)**2).to("m^2")
                A_tube_out = (pi * (Do_t/2)**2).to("m^2")

                cold_wet_P = (Nt * pi * Do_t).to("m")
                cold_flow_A = A_cross
                cold_Dh = (4 * cold_flow_A / cold_wet_P).to("m")

                hot_wet_P = (Nt * pi * Di_t).to("m")
                hot_flow_A = (Nt * (pi * (Di_t/2)**2)).to("m^2")
                hot_Dh = (4 * hot_flow_A / hot_wet_P).to("m")

                spec.update({
                    "hot_wet_P": hot_wet_P, "hot_flow_A": hot_flow_A, "hot_Dh": hot_Dh,
                    "cold_wet_P": cold_wet_P, "cold_flow_A": cold_flow_A, "cold_Dh": cold_Dh
                })
                out.append(replace(stg, spec=spec))


            elif stg.kind == "reversal_chamber":
                spec = dict(stg.spec)

                Di_t = spec["inner_diameter"].to("m")
                t = spec["wall_t"].to("m")
                Do_t = (Di_t + 2*t).to("m")

                Di_drum = self.drum.Di.to("m")
                A_drum = (pi * (Di_drum/2)**2).to("m^2")
                A_tube_out = (pi * (Do_t/2)**2).to("m^2")

                spec["outer_diameter"] = Do_t
                spec["roughness_cold_surface"] = spec["roughness_out"]

                hot_wet_P = (pi * Di_t).to("m")
                hot_flow_A = (pi * (Di_t/2)**2).to("m^2")
                hot_Dh = (4 * hot_flow_A / hot_wet_P).to("m")

                cold_wet_P = (pi * Do_t).to("m")
                cold_flow_A = (A_drum - A_tube_out).to("m^2")
                cold_Dh = (4 * cold_flow_A / cold_wet_P).to("m")

                spec.update({
                    "hot_wet_P": hot_wet_P, "hot_flow_A": hot_flow_A, "hot_Dh": hot_Dh,
                    "cold_wet_P": cold_wet_P, "cold_flow_A": cold_flow_A, "cold_Dh": cold_Dh
                })
                out.append(replace(stg, spec=spec))


            elif stg.kind == "economiser":
                spec = dict(stg.spec)

                Di_t = spec["inner_diameter"].to("m")
                t = spec["wall_t"].to("m")
                Nt = spec["tubes_number"].to("")
                Do_t = (Di_t + 2*t).to("m")

                Di_drum = self.drum.Di.to("m")
                A_drum = (pi * (Di_drum/2)**2).to("m^2")
                A_tube_out = (pi * (Do_t/2)**2).to("m^2")

                spec["outer_diameter"] = Do_t
                spec["roughness_cold_surface"] = spec["roughness_in"]

                hot_wet_P = (Nt * pi * Do_t).to("m")
                hot_flow_A = (A_drum - Nt * A_tube_out).to("m^2")
                hot_Dh = (4 * hot_flow_A / hot_wet_P).to("m")

                cold_wet_P = (Nt * pi * Di_t).to("m")
                cold_flow_A = (Nt * (pi * (Di_t/2)**2)).to("m^2")
                cold_Dh = (4 * cold_flow_A / cold_wet_P).to("m")

                spec.update({
                    "hot_wet_P": hot_wet_P, "hot_flow_A": hot_flow_A, "hot_Dh": hot_Dh,
                    "cold_wet_P": cold_wet_P, "cold_flow_A": cold_flow_A, "cold_Dh": cold_Dh
                })
                out.append(replace(stg, spec=spec))
            else:
                raise ValueError("unknown stage kind")
        return out

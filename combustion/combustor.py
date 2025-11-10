from combustion.adiabatic_flame_temperature import adiabatic_flame_T
from common.results import CombustionResult
from common.models import GasStream
from common.units import Q_
from combustion.heat import total_input_heat
from combustion.flue import air_flow_rates
from combustion.adiabatic_flame_temperature import adiabatic_flame_T


class Combustor:
    def __init__(self, air: GasStream, fuel: GasStream, excess_air_ratio: Q_):
        self.air = air
        self.fuel = fuel
        self.excess_air_ratio = excess_air_ratio

    def run(self):
        air = self.air
        fuel = self.fuel
        flue = GasStream(
            mass_flow= Q_(0, "kelvin"),
            T = air.T,
            P = air.P,
            comp = {}
        )

        air.mass_flow = air_flow_rates(air, fuel, self.excess_air_ratio)

        power_LHV, Q_in = total_input_heat(fuel, air)


        # compute flue directly from air and fuel streams
        flue = adiabatic_flame_T(air, fuel)
        T_ad = flue.T


        flue.T = T_ad
        flue.P = air.P

        return CombustionResult(power_LHV, Q_in, T_ad, flue)

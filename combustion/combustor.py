from combustion.adiabatic_flame_temperature import adiabatic_flame_T
from common.results import CombustionResult
from common.models import GasStream
from common.units import Q_
from combustion.heat import total_input_heat
from combustion.flue import from_fuel_and_air, air_flow_rates
from combustion.mass_mole import to_mole, molar_flow

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

        flue.comp, flue.mass_flow = from_fuel_and_air(fuel, air)

        power_LHV, Q_in = total_input_heat(fuel, air)

        T_ad = adiabatic_flame_T(flue, Q_in)

        flue.T = T_ad
        flue.P = air.P

        return CombustionResult(power_LHV, Q_in, T_ad, flue)

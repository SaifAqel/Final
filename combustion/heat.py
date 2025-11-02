from common.units import ureg, Q_
import re
from common.props import WaterProps, GasProps
from common.models import GasStream
from combustion.mass_mole import mix_molar_mass, to_mole
from common.constants import formation_enthalpies, molar_masses, T_ref, P_ref

_dHf = formation_enthalpies
_gasprops = GasProps()

def parse_CH(s: str):
    m = re.fullmatch(r'C(\d*)H(\d+)', s)
    if not m: return None, None
    C = int(m.group(1)) if m.group(1) else 1
    H = int(m.group(2))
    return C, H

def compute_LHV_HHV(fuel: GasStream) -> Q_:
    latent_H2O = WaterProps.h_g(P_ref) - WaterProps.h_f(P_ref)
    # product enthalpies (kJ/mol)
    H2O_liq = _dHf["H2O"]                                # kJ/mol
    H2O_vap = _dHf["H2O"] + latent_H2O * molar_masses["H2O"]          # (kJ/kg)*(kg/mol) = kJ/mol

    # initialize as quantities, not floats
    react = 0 * _dHf["CO2"]                              # kJ/mol
    HHV_p = 0 * _dHf["CO2"]
    LHV_p = 0 * _dHf["CO2"]

    mol_comp = to_mole(fuel.comp)

    for comp, x in mol_comp.items():
        # use quantity zero as default to avoid dimensionless 0.0
        dh = _dHf.get(comp, 0 * _dHf["CO2"])              # kJ/mol
        react += x * dh

        C, H = parse_CH(comp)
        if C is not None:
            HHV_p += x * (C * _dHf["CO2"] + (H/2) * H2O_liq)
            LHV_p += x * (C * _dHf["CO2"] + (H/2) * H2O_vap)
        elif comp == "H2S":
            HHV_p += x * (_dHf["SO2"] + H2O_liq)
            LHV_p += x * (_dHf["SO2"] + H2O_vap)
        else:
            HHV_p += x * dh
            LHV_p += x * dh

    M_mix = mix_molar_mass(mol_comp)  # kg/mol

    HHV_mol = react - HHV_p           # kJ/mol
    LHV_mol = react - LHV_p           # kJ/mol
    HHV_kg  = HHV_mol / M_mix         # kJ/kg
    LHV_kg  = LHV_mol / M_mix         # kJ/kg

    P_HHV = (HHV_kg * fuel.mass_flow).to('kW')
    P_LHV = (LHV_kg * fuel.mass_flow).to('kW')                     # (kJ/kg)*(kg/s)=kJ/s = kW

    return HHV_kg.to('kJ/kg'), LHV_kg.to('kJ/kg'), P_HHV, P_LHV


def sensible_heat(stream: GasStream) -> Q_:
    s = stream
    cp = (_gasprops.cp(s.T, s.P, s.comp)).to('kJ/(kg*K)')
    dT = ((s.T) - (T_ref)).to('K')
    return (s.mass_flow * cp * dT).to('kW')

def total_input_heat(fuel, air) -> Q_:
    _,_,_,power_LHV = compute_LHV_HHV(fuel)
    fuel_sens = sensible_heat(fuel)
    air_sens  = sensible_heat(air)
    Q_in = ((power_LHV).to('kW') +
            (fuel_sens).to('kW') +
            (air_sens).to('kW')).to('kW')
    return power_LHV, Q_in

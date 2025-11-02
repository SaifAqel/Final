from common.models import GasStream
from common.units import Q_
from common.constants import O2_per_mol
from combustion.mass_mole import to_mole, molar_flow, mass_flow, to_mass, mix_molar_mass

def stoich_O2_required_per_mol_fuel(fuel: GasStream) -> Q_:
    fuel_x = to_mole(fuel.comp)
    total = sum(fuel_x[k] * O2_per_mol.get(k, 0.0) for k in fuel_x)
    return Q_(total, "dimensionless")

def air_flow_rates(air: GasStream, fuel: GasStream, excess: Q_) -> Q_:
    air_x = to_mole(air.comp)
    fuel_n_dot = molar_flow(fuel.comp, fuel.mass_flow)
    O2_x = air_x["O2"]
    O2_req = stoich_O2_required_per_mol_fuel(fuel)
    O2_stoich = fuel_n_dot * O2_req
    O2_actual = O2_stoich * excess
    air_n = O2_actual / O2_x
    M_air = mix_molar_mass(air_x)
    air_m = air_n * M_air
    return air_m

def from_fuel_and_air(fuel: GasStream, air: GasStream) -> GasStream:
    O2_req = stoich_O2_required_per_mol_fuel(fuel)
    fuel_x = to_mole(fuel.comp)
    fuel_n = molar_flow(fuel.comp, fuel.mass_flow)
    air_x = to_mole(air.comp)
    air_n = molar_flow(air.comp, fuel.mass_flow)

    gf=lambda k:fuel_x.get(k,0.0); ga=lambda k:air_x.get(k,0.0)
    n_CO2 = air_n*ga("CO2")+fuel_n*gf("CO2")+fuel_n*(gf("CH4")+2*gf("C2H6")+3*gf("C3H8")+4*gf("C4H10"))
    n_H2O = fuel_n*gf("H2O")+fuel_n*(2*gf("CH4")+3*gf("C2H6")+4*gf("C3H8")+5*gf("C4H10"))+fuel_n*gf("H2S")+air_n*ga("H2O")
    n_SO2 = fuel_n*gf("H2S")
    n_O2  = air_n*ga("O2") - fuel_n*O2_req
    n_N2  = air_n*ga("N2") + fuel_n*gf("N2")
    n_Ar  = air_n*ga("Ar")
    flows={"CO2":n_CO2,"H2O":n_H2O,"SO2":n_SO2,"O2":n_O2,"N2":n_N2,"Ar":n_Ar}

    n_tot=sum(flows.values())
    mol_comp={k:(v/n_tot if n_tot!=0 else 0.0) for k,v in flows.items()}
    mass_comp = to_mass(mol_comp)
    m_dot = mass_flow(mol_comp, n_tot)
    
    return mass_comp, m_dot

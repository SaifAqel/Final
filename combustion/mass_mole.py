from typing import Dict
from common.constants import molar_masses

def to_mole(mass_comp: Dict[str, float]) -> Dict[str, float]:
    n = {sp: mass_comp[sp] / molar_masses[sp] for sp in mass_comp}
    tot = sum(n.values())
    return {sp: v / tot for sp, v in n.items()}

def to_mass(mol_comp: Dict[str, float]) -> Dict[str, float]:
    m = {sp: mol_comp[sp] * molar_masses[sp] for sp in mol_comp}
    tot = sum(m.values())
    return {sp: v / tot for sp, v in m.items()}

def mix_molar_mass(mol_comp: Dict[str, float]) -> float:
    return sum(mol_comp[sp] * molar_masses[sp] for sp in mol_comp)

def molar_flow(mass_comp: Dict[str, float], m_dot: float) -> float:
    return sum((mass_comp[sp] * m_dot) / (molar_masses[sp]) for sp in mass_comp)

def mass_flow(mol_comp: Dict[str, float], n_dot: float) -> float:
    return sum(mol_comp[sp] * n_dot * molar_masses[sp] for sp in mol_comp)

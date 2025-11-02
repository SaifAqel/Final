from common.units import Q_

formation_enthalpies = {
    "CH4": Q_(-74.8, "kilojoule / mole"),
    "C2H6": Q_(-84.7, "kilojoule / mole"),
    "C3H8": Q_(-103.8, "kilojoule / mole"),
    "C4H10": Q_(-126.1, "kilojoule / mole"),
    "SO2": Q_(-296.8, "kilojoule / mole"),
    "CO2": Q_(-393.5, "kilojoule / mole"),
    "H2O": Q_(-285.5, "kilojoule / mole"),
}

T_ref = Q_(298.15, "kelvin")
P_ref = Q_(101325, "pascal")

molar_masses = {
    "CH4": Q_(0.01604, "kilogram / mole"),
    "C2H6": Q_(0.03007, "kilogram / mole"),
    "C3H8": Q_(0.04410, "kilogram / mole"),
    "C4H10": Q_(0.05812, "kilogram / mole"),
    "H2S": Q_(0.03408, "kilogram / mole"),
    "N2": Q_(0.02802, "kilogram / mole"),
    "CO2": Q_(0.04401, "kilogram / mole"),
    "H2O": Q_(0.018015, "kilogram / mole"),
    "O2": Q_(0.031998, "kilogram / mole"),
    "Ar": Q_(0.039948, "kilogram / mole"),
    "SO2": Q_(0.06407, "kilogram / mole"),
}

O2_per_mol = {
    "CH4": Q_(2.0, "dimensionless"),
    "C2H6": Q_(3.5, "dimensionless"),
    "C3H8": Q_(5.0, "dimensionless"),
    "C4H10": Q_(6.5, "dimensionless"),
    "H2S": Q_(1.0, "dimensionless"),
    "N2": Q_(0.0, "dimensionless"),
    "CO2": Q_(0.0, "dimensionless"),
    "H2O": Q_(0.0, "dimensionless"),
}

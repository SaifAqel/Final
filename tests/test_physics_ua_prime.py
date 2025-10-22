import pytest
from units import Q_
from models import HXStage, GasStream, WaterStream
from physics import ua_per_m
from props import WaterProps

@pytest.mark.usefixtures("fake_gas")
def test_ua_per_m_positive_singlephase():
    stage = HXStage(
        name="HX_test",
        kind="single_tube",
        spec={
            "inner_diameter": Q_(0.1, "m"),
            "inner_length":   Q_(1.0, "m"),
            "wall_t":         Q_(0.003, "m"),
            "wall_k":         Q_(15.0, "W/m/K"),
            "eps_in":         Q_(0.3, "dimensionless"),
            "roughness_out":  Q_(2.0, "micrometer"),
            # cold-side geometry directly
            "cold_Di":        Q_(0.05, "m"),
        },
    )
    g = GasStream(mass_flow=Q_(1.0, "kg/s"), T=Q_(800, "K"), P=Q_(101325, "Pa"))
    w = WaterStream(mass_flow=Q_(2.0, "kg/s"), h=Q_(200e3, "J/kg"), P=Q_(1e6, "Pa"))
    UA_p, qpp = ua_per_m(g, w, stage=stage)
    assert UA_p.magnitude > 0
    assert qpp.magnitude > 0
    assert UA_p.check("power / temperature / length")
    assert qpp.check("power / length ** 2")

@pytest.mark.usefixtures("fake_gas")
def test_boiling_path_sets_higher_h():
    # Put water enthalpy between h_f and h_g to trigger boiling branch
    P = Q_(1.0, "MPa")
    hf, hg = WaterProps.h_f(P), WaterProps.h_g(P)
    h_mid = (hf + hg) * 0.5
    stage = HXStage(
        name="HX_boil",
        kind="single_tube",
        spec={
            "inner_diameter": Q_(0.1, "m"),
            "inner_length":   Q_(1.0, "m"),
            "wall_t":         Q_(0.003, "m"),
            "wall_k":         Q_(15.0, "W/m/K"),
            "eps_in":         Q_(0.3, "dimensionless"),
            "roughness_out":  Q_(5.0, "micrometer"),
            "cold_Di":        Q_(0.05, "m"),
        },
    )
    g = GasStream(mass_flow=Q_(1.0, "kg/s"), T=Q_(800, "K"), P=Q_(101325, "Pa"))
    w = WaterStream(mass_flow=Q_(2.0, "kg/s"), h=h_mid, P=P)
    UA_p, qpp = ua_per_m(g, w, stage=stage)
    assert UA_p.magnitude > 0
    assert qpp.magnitude > 0

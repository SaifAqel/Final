import pytest
from units import Q_
from models import HXStage, GasStream, WaterStream
from solver import StageSolver
from physics import cp_gas
from props import WaterProps

@pytest.mark.usefixtures("fake_gas")
def test_stage_energy_balance_and_monotonicity():
    stage = HXStage(
        name="HX_energy",
        kind="single_tube",
        spec={
            "inner_diameter": Q_(0.1, "m"),
            "inner_length":   Q_(2.0, "m"),
            "wall_t":         Q_(0.003, "m"),
            "wall_k":         Q_(15.0, "W/m/K"),
            "eps_in":         Q_(0.3, "dimensionless"),
            "roughness_out":  Q_(2.0, "micrometer"),
            "cold_Di":        Q_(0.05, "m"),
        },
    )
    g0 = GasStream(mass_flow=Q_(1.0, "kg/s"), T=Q_(900, "K"), P=Q_(101325, "Pa"))
    w0 = WaterStream(mass_flow=Q_(2.0, "kg/s"), h=Q_(200e3, "J/kg"), P=Q_(1e6, "Pa"))

    g_hist, w_hist = StageSolver(stage, g0, w0).solve(N=12)

    # Monotonicity
    g_T = [g.T.to("K").magnitude for g in g_hist]
    w_h = [w.h.to("J/kg").magnitude for w in w_hist]
    assert all(g_T[i+1] <= g_T[i] for i in range(len(g_T)-1))
    assert all(w_h[i+1] >= w_h[i] for i in range(len(w_h)-1))

    # Energy balance over the whole stage
    g_in, g_out = g_hist[0], g_hist[-1]
    w_out, w_in = w_hist[0], w_hist[-1]  # counterflow history order
    cp_g = cp_gas(g_in)  # fake gas cp is constant
    q_hot = (g_in.mass_flow * cp_g * (g_in.T - g_out.T)).to("W")
    q_cold = (w_in.mass_flow * (w_out.h - w_in.h) / (Q_(1.0, "s"))).to("W")
    err = abs((q_hot.magnitude - q_cold.magnitude) / max(q_hot.magnitude, 1e-9))
    assert err < 0.05  # ≤5%

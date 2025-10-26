from math import log, sqrt, exp, log10
from units import Q_
from models import WaterStream, HXStage
from props import WaterProps
from typing import Tuple


P_CRIT_WATER = Q_(22.064, "MPa")
MW_WATER = 18.01528

def velocity(w: WaterStream, Aflow, umax_factor=None) -> Q_:
    rho = WaterProps.rho_from_Ph(w.P, w.h)
    u = (w.mass_flow / (rho * Aflow)).to("m/s")
    if umax_factor is not None:
        return (umax_factor * u).to("m/s")
    return u

def reynolds_number(w:WaterStream, Aflow, char_len, umax_factor=None):
    rho = WaterProps.rho_from_Ph(w.P, w.h)
    v = velocity(w, Aflow, umax_factor)
    mu = WaterProps.mu_from_Ph(w.P, w.h)
    return (rho * v * char_len) / mu

def prandtl_number(cp: Q_, mu: Q_, k: Q_) -> Q_:
    return cp * mu / k

def film_temp(T_bulk: Q_, T_wall: Q_) -> Q_:
    return 0.5 * (T_bulk + T_wall)

def _is_boiling(P, h, T_wall: Q_ | None = None) -> bool:
    hf = WaterProps.h_f(P)
    hg = WaterProps.h_g(P)
    if hf <= h <= hg:
        return True
    if (h < hf) and (T_wall is not None):
        return T_wall > WaterProps.Tsat(P) + Q_(3, "K")
    return False

def pr(w:WaterStream) -> Q_:
    cp = WaterProps.cp_from_Ph(w.P, w.h)
    mu = WaterProps.mu_from_Ph(w.P, w.h)
    k = WaterProps.k_from_Ph(w.P, w.h)
    return prandtl_number(cp, mu, k)

def pr_s(w: WaterStream, T_wall: Q_) -> Q_:
    cp_s = WaterProps.cp_from_PT(w.P, T_wall)
    mu_s = WaterProps.mu_from_PT(w.P, T_wall)
    k_s = WaterProps.k_from_PT(w.P, T_wall)
    return prandtl_number(cp_s, mu_s, k_s)

def nu_zukauskas_bank(Re: Q_, Pr: Q_, Pr_s: Q_, arrangement: str) -> tuple[Q_, Q_]:
    # Zukauskas correlation for crossflow over tube banks, ε = (Pr/Pr_s)^0.25
    # Coeffs by Re band and arrangement
    Re = Re.to("").magnitude
    Pr = Pr.to("").magnitude
    Pr_s = Pr_s.to("").magnitude

    bands = [
        (1e3, 2e3, {"inline": (0.90, 0.40), "staggered": (1.04, 0.40)}),
        (2e3, 4e3, {"inline": (0.52, 0.50), "staggered": (0.71, 0.50)}),
        (4e3, 1e5, {"inline": (0.27, 0.63), "staggered": (0.35, 0.60)}),
        (1e5, 2e6, {"inline": (0.021, 0.84), "staggered": (0.022, 0.84)}),
    ]
    C = None; m = None
    for Re_min, Re_max, table in bands:
        if Re_min <= Re < Re_max:
            C, m = table["staggered" if arrangement == "staggered" else "inline"]
            break
    if C is None:
        raise ValueError(f"Re={Re:.0f} outside Zukauskas bands")
    n = 0.36 if Pr <= 10.0 else 0.25
    s = 0.25
    nu = C * (Re**m) * (Pr**n) * ((Pr / max(Pr_s, 1e-12))**s)
    return Q_(nu, ""),Q_(m, "")



def nu_churchill_bernstein(Re: Q_, Pr: Q_) -> Q_:
    # External crossflow over a single cylinder
    Re = Re.to("").magnitude
    Pr = Pr.to("").magnitude
    a = 0.3
    b = (0.62 * Re**0.5 * Pr**(1/3)) / (1 + (0.4/Pr)**(2/3))**0.25
    c = (1 + (Re/282000.0)**(5/8))**(4/5)
    return Q_(a + b * c, "")



def nu_gnielinski(Re: Q_, Pr: Q_, mu_ratio: Q_, L: Q_, D: Q_) -> Q_:
    Re = Re.to("").magnitude
    Pr = Pr.to("").magnitude
    L = L.to("meter").magnitude
    D = D.to("meter").magnitude
    if Re < 2300.0:
        Gz = Re * Pr * (D / max(L, 1e-12))
        return Q_(3.66 + (0.0668 * Gz) / (1 + 0.04 * Gz**(2/3)), "")
    f = (0.79 * log(Re) - 1.64) ** -2
    num = (f/8) * (Re - 1000.0) * Pr
    den = 1 + 12.7 * (f/8)**0.5 * (Pr**(2/3) - 1)
    Nu = num / max(den, 1e-12)
    mu_ratio = mu_ratio.to("").magnitude
    return Q_(Nu * (mu_ratio ** 0.11), "")


def compute_nusselt(w:WaterStream, stage: HXStage, T_wall: Q_) -> Q_:

    if stage.kind == "single_tube":
        L = stage.spec["outer_diameter"]
        Aflow = stage.spec["cold_flow_A"]
        umax = stage.spec.get("umax_factor")
        Re = reynolds_number(w, Aflow, L, umax)
        Pr = pr(w)
        return nu_churchill_bernstein(Re, Pr)

    if stage.kind == "tube_bank":
        L = stage.spec["outer_diameter"]
        Aflow = stage.spec["cold_flow_A"]
        umax = stage.spec.get("umax_factor")
        Re = reynolds_number(w, Aflow, L, umax)
        Pr = pr(w)
        Pr_s = pr_s(w, T_wall)
        N_rows = stage.spec["N_rows"]
        ST = stage.spec["ST"]
        SL = stage.spec["SL"]
        arrangement = stage.spec["arrangement"]
        Nu, m = nu_zukauskas_bank(Re, Pr, Pr_s, arrangement)
        Nu *= bank_row_factor(N_rows)
        Nu *= spacing_factor(L, ST, SL, arrangement, m)
        return Nu

    if stage.kind == "reversal_chamber":
        L = stage.spec["outer_diameter"]
        Aflow = stage.spec["cold_flow_A"]
        umax = stage.spec.get("umax_factor")
        Re = reynolds_number(w, Aflow, L, umax)
        Pr = pr(w)
        Rc = stage.spec["curvature_radius"]
        # approximate as single cylinder with optional bump
        Nu = nu_churchill_bernstein(Re, Pr)
        return Nu * bend_factor_external(L, Rc)

    if stage.kind == "economiser":
        D = stage.spec["inner_diameter"]
        L = stage.spec["inner_length"]
        Aflow = stage.spec["cold_flow_A"]
        T_bulk = WaterProps.T_from_Ph(w.P, w.h)
        umax = stage.spec.get("umax_factor")
        Re = reynolds_number(w, Aflow, D, umax)
        Pr = pr(w)
        mu_ratio = _mu_ratio(w, T_bulk, T_wall)
        return nu_gnielinski(Re, Pr, mu_ratio, L, D)

    raise ValueError(f"unknown stage kind: {stage.kind}")



# ---------- viscosity ratio (bulk/wall) ----------
def _mu_ratio(w: WaterStream, T_bulk: Q_, T_wall: Q_) -> Q_:
    mu_b = WaterProps.mu_from_PT(w.P, T_bulk)
    mu_w = WaterProps.mu_from_PT(w.P, T_wall)
    return mu_b / mu_w

# ---------- bend factor for external U-bend / reversal chamber ----------
def bend_factor_external(D: Q_, Rc: Q_) -> Q_:
    """
    Simple curvature boost for crossflow around a bent tube.
    D = tube OD [m], Rc = bend centerline radius [m].
    Returns multiplier in [1.0, 1.25].
    """
    if Rc <= 0 or D <= 0:
        return 1.0
    phi = 1.0 + 0.10 * sqrt(D / Rc)   # ~+10% at tight bends; milder for gentle bends
    return Q_(phi, "")

def spacing_factor(D: Q_, ST: Q_, SL: Q_, arrangement: str, m_exp: Q_) -> Q_:
    # Max-velocity ratio estimates (engineering approximations)
    if arrangement == "staggered":
        denom_T = ST - D
        denom_L = SL - (0.5 * D)
        vmax_ratio = (ST / denom_T) * (SL / denom_L)
        # slightly soften for very open banks
        vmax_ratio = vmax_ratio**0.5
    else:
        # inline: minimum area between neighboring cylinders
        vmax_ratio = ST / (ST - D)
    m_exp = m_exp.to("").magnitude
    phi = vmax_ratio ** m_exp
    return phi

def bank_row_factor(N_rows: Q_) -> Q_:
    n = N_rows.to("").magnitude
    f = 1.0 - 0.30 * exp(-0.30 * n)
    return Q_(f, "")

def _h_water_singlephase(w: WaterStream, stage: HXStage, T_wall) -> Q_:
    Nu = compute_nusselt(w, stage, T_wall)
    k = WaterProps.k_from_Ph(w.P, w.h)
    if stage.kind in ("single_tube", "tube_bank", "reversal_chamber"):
        Dh = stage.spec["outer_diameter"]
    else:
        Dh = stage.spec["inner_diameter"]
    return (Nu * k / Dh).to("W/m^2/K")

def _h_water_boil_cooper(P: Q_, qpp: Q_, Rp: Q_) -> Q_:
    p_r = (P.to("MPa") / P_CRIT_WATER).magnitude
    Rp_um = Rp.to("micrometer").magnitude
    q_kWm2 = qpp.to("kW/m^2").magnitude
    h_kWm2K = 55.0 * (p_r**0.12) * ((-log10(Rp_um))**-0.55) * (MW_WATER**-0.5) * (q_kWm2**0.67)
    return Q_(h_kWm2K, "kW/m^2/K").to("W/m^2/K")

def water_htc(w: WaterStream, stage: HXStage, T_wall: Q_, qpp: Q_) -> tuple[Q_, bool]:
    boiling = _is_boiling(w.P, w.h, T_wall)
    if boiling:
        h_lo = _h_liquid_only(w, stage, T_wall)
        h_nb = _h_water_boil_cooper(w.P, qpp, stage.spec["roughness_cold_surface"])
        T_sat = WaterProps.Tsat(w.P)
        mu_l  = WaterProps.mu_from_PT(w.P, T_sat)
        A = stage.spec["cold_flow_A"]
        Dh = stage.spec["cold_Dh"]
        G = _mass_flux(w, A)
        h_lv = WaterProps.h_g(w.P) - WaterProps.h_f(w.P)
        x = WaterProps.quality_from_Ph(w.P, w.h)
        Re_lo = (G * Dh / mu_l).to("")
        S = _chen_S_factor(qpp, G, h_lv, Re_lo)
        if x is not None:
            F = _chen_F_factor(w.P, x)
        else:
            F = 1
        h_c = F * h_lo + S * h_nb
    else:
        h_c = _h_water_singlephase(w, stage, T_wall)
    return h_c, boiling

def _mass_flux(w: WaterStream, Aflow: Q_) -> Q_:
    return (w.mass_flow / Aflow).to("kg/m^2/s")

def _h_liquid_only(w: WaterStream, stage: HXStage, T_wall: Q_) -> Q_:
    D_h = stage.spec["cold_Dh"]
    L   = stage.spec["inner_length"]
    A   = stage.spec["cold_flow_A"]
    T_sat = WaterProps.Tsat(w.P)
    mu_l  = WaterProps.mu_from_PT(w.P, T_sat)
    k_l   = WaterProps.k_from_PT(w.P, T_sat)
    cp_l  = WaterProps.cp_from_PT(w.P, T_sat)
    G   = _mass_flux(w, A)
    Re_lo  = (G * D_h / mu_l).to("")
    Pr  = prandtl_number(cp_l, mu_l, k_l).to("")
    # wall/bulk viscosity ratio at Tsat
    mu_ratio = (mu_l / WaterProps.mu_from_PT(w.P, T_wall)).to("")

    if stage.kind == "economiser":
        Nu = nu_gnielinski(Re_lo, Pr, mu_ratio, L, D_h)
    elif stage.kind == "tube_bank":
        # need Pr_s at wall
        Pr_s = prandtl_number(cp_l, WaterProps.mu_from_PT(w.P, T_wall), WaterProps.k_from_PT(w.P, T_wall))
        Nu, m = nu_zukauskas_bank(Re_lo, Pr, Pr_s, stage.spec["arrangement"])
        Nu *= bank_row_factor(stage.spec["N_rows"])
        Nu *= spacing_factor(D_h, stage.spec["ST"], stage.spec["SL"], stage.spec["arrangement"], m)
    else:
        Nu = nu_churchill_bernstein(Re_lo, Pr)

    return (Nu * k_l / D_h).to("W/m^2/K")

def _martinelli_Xtt(P: Q_, x: float) -> float:
    T_sat = WaterProps.Tsat(P)
    rho_l = WaterProps.rho_from_Px(P, Q_(0.0, ""))    # saturated liquid
    rho_g = WaterProps.rho_from_Px(P, Q_(1.0, ""))    # saturated vapor
    mu_l  = WaterProps.mu_from_PT(P, T_sat)
    mu_g  = WaterProps.mu_from_PT(P, T_sat)           # IAPWS μ(T,P) handles both phases
    mu_ratio  = (mu_l / mu_g).to("").magnitude
    rho_ratio = (rho_g / rho_l).to("").magnitude
    return ((1 - x) / x) ** 0.9 * (mu_ratio ** 0.1) * (rho_ratio ** 0.5)

def _chen_S_factor(qpp: Q_, G: Q_, h_lv: Q_, Re_lo: Q_) -> Q_:
    Re = max(1.0, Re_lo.to("").magnitude)
    S = 1.0 / (1.0 + 2.53e-6 * (Re ** 1.17))
    return Q_(max(0.1, min(S, 1.0)), "")

def _chen_F_factor(P: Q_, x: float) -> Q_:
    Xtt = _martinelli_Xtt(P, x)
    F = 1.0 + 0.12 * (max(1e-6, 1.0 / Xtt) ** 0.8)
    return Q_(min(5.0, max(1.0, F)), "")

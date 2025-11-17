from common.results import StepResult
from common.units import Q_
from common.models import HXStage, GasStream, WaterStream
from heat.physics import wall_resistance, fouling_resistances
from heat.water_htc import water_htc
from common.props import WaterProps
from heat.gas_htc import gas_htc, gas_htc_parts

def solve_step(g: GasStream, w: WaterStream, stage: HXStage, Tgw_guess: Q_, Tww_guess: Q_, qprime_guess: Q_, i: int, x: Q_, dx: Q_) -> StepResult:
    spec = stage.spec
    Pg = spec["hot_wet_P"]
    Pw = spec["cold_wet_P"]
    Tg = g.T
    if stage.spec["pool_boiling"]:
        # bulk water = pool at Tsat
        Tw = WaterProps.Tsat(w.P)
    else:
        Tw = WaterProps.T_from_Ph(w.P, w.h)
    Tgw = Tgw_guess
    Tww = Tww_guess
    qprime = qprime_guess
    alpha = 0.25
    tolT = Q_(1e-3,"K"); tolq = Q_(1e-3,"W/m"); maxit = 10

    for _ in range(maxit):
        h_g = gas_htc(g, spec, Tgw)
        qpp_cold = (qprime / Pw).to("W/m^2")
        h_c, boiling = water_htc(w, stage, Tww, qpp_cold)
        Rfg, Rfc = fouling_resistances(spec)
        Rw = wall_resistance(spec)
        Rg = (1/(h_g*Pg)).to("K*m/W")
        Rc = (1/(h_c*Pw)).to("K*m/W")

        UA_prime = (1/(Rg + Rfg + Rw + Rfc + Rc)).to("W/K/m")

        qprime_new = (UA_prime * (Tg - Tw)).to("W/m")

        qpp_hot  = (qprime_new / Pg).to("W/m^2")
        qpp_cold = (qprime_new / Pw).to("W/m^2")

        Tgw_new = (Tg - qpp_hot/h_g - qpp_hot*Rfg*Pg).to("K")
        Tww_new = (Tw + qpp_cold*Rw*Pw + qpp_cold*Rfc*Pw + qpp_cold/h_c).to("K")

        dTgw = abs(Tgw_new - Tgw); dTww = abs(Tww_new - Tww); dq = abs(qprime_new - qprime)
        if dTgw < tolT and dTww < tolT and dq < tolq:
            Tgw, Tww, qprime = Tgw_new, Tww_new, qprime_new
            break

        Tgw = (alpha*Tgw_new + (1-alpha)*Tgw).to("K")
        Tww = (alpha*Tww_new + (1-alpha)*Tww).to("K")
        qprime = (alpha*qprime_new + (1-alpha)*qprime).to("W/m")
    # After convergence, decompose gas-side HTC into convective and radiative parts
    h_conv, h_rad = gas_htc_parts(g, spec, Tgw)
    h_g = (h_conv + h_rad).to("W/m^2/K")  # total gas-side HTC, keep as before for diagnostics

    # Split q′ between convection and radiation based on HTC ratio
    h_tot_mag = h_g.to("W/m^2/K").magnitude
    if h_tot_mag > 0:
        frac_conv = (h_conv / h_g).to("").magnitude
        frac_conv = max(0.0, min(1.0, frac_conv))
    else:
        frac_conv = 0.0

    qprime_conv = (qprime * frac_conv).to("W/m")
    qprime_rad  = (qprime - qprime_conv).to("W/m")

    return StepResult(
        i=i, x=x, dx=dx,
        gas=g, water=w,
        Tgw=Tgw, Tww=Tww,
        UA_prime=UA_prime,
        qprime=qprime,
        boiling=boiling,
        h_g=h_g,
        h_c=h_c,
        qprime_conv=qprime_conv,
        qprime_rad=qprime_rad,
    )
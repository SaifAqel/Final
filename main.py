from new_loader import load_all
from combustion.combustor import Combustor
from heat.runner import run_hx
from common.units import Q_
from common.props import WaterProps
from common.models import WaterStream
from common.results import CombustionResult
from common.results import write_results_csvs
import logging

logger_name = "main"
log = logging.getLogger(logger_name)

def _water_mass_from_efficiency(
    eta: Q_,
    combustion: CombustionResult,
    water_template: WaterStream,
) -> Q_:
    """
    Compute feedwater/steam mass flow from a boiler efficiency guess.

    eta           : dimensionless (e.g. Q_(0.9, ""))
    combustion    : CombustionResult (provides Q_in)
    water_template: WaterStream with inlet h,P (mass_flow ignored here)

    Uses:
      Q_target = eta * Q_in
      m_w      = Q_target / (h_steam - h_feedwater)

    Steam is taken as saturated vapor at the feedwater pressure.
    """

    # Total input heat from combustion (kW → W)
    Q_in = combustion.Q_in.to("W")  # Q_

    # Feedwater inlet enthalpy
    h_in = water_template.h.to("J/kg")

    # Target steam enthalpy: saturated vapor at the feedwater pressure
    h_steam = WaterProps.h_g(water_template.P).to("J/kg")

    delta_h = (h_steam - h_in).to("J/kg")
    if delta_h.magnitude <= 0.0:
        raise ValueError("Steam enthalpy must be greater than feedwater enthalpy.")

    Q_target = (eta * Q_in).to("W")  # dimensionless * W

    m_w = (Q_target / delta_h).to("kg/s")

    return m_w


stages, air, fuel, water, drum, operation = load_all(
    stages_path="config/stages.yaml",
    air_path="config/air.yaml",
    fuel_path="config/fuel.yaml",
    water_path="config/water.yaml",
    drum_path="config/drum.yaml",
    operation_path="config/operation.yaml",
)

svc = Combustor(air, fuel, operation["excess_air_ratio"])
combustion_results = svc.run()
log.info(f"Combustion results: {combustion_results}")


# Use water from config only as template for h and P
water_template: WaterStream = water

# Iteration settings
eta_guess = Q_(0.90, "")            # 90% initial efficiency guess
tol_m = Q_(1e-3, "kg/s")            # convergence tolerance on mass flow
max_iter = 20

prev_m = None
final_result = None
final_m = None
final_eta = None

for it in range(max_iter):
    if it % 5 == 0:
        log.info(f"Iteration {it}: eta_guess={eta_guess}, prev_m={prev_m}")
    # 1) Compute water mass flow from current efficiency guess
    m_w = _water_mass_from_efficiency(eta_guess, combustion_results, water_template)

    # 2) Build a WaterStream for this iteration
    water_iter = WaterStream(
        mass_flow=m_w,
        h=water_template.h,
        P=water_template.P,
    )

    # 3) Run the heat exchanger with this water flow (no CSV spam)
    final_result = run_hx(
        stages_raw=stages,
        water=water_iter,
        gas=combustion_results.flue,
        drum=drum,
        target_dx="0.1 m",
        combustion=combustion_results,
        write_csv=False,
    )

    # 4) Read actual efficiency from TOTAL_BOILER row (indirect efficiency)
    total_row = next(
        r for r in final_result["summary_rows"]
        if r["stage_name"] == "TOTAL_BOILER"
    )

    eta_indirect = total_row["η_indirect[-]"]
    if eta_indirect == "" or eta_indirect is None:
        # Fallback: no efficiency available; stop
        final_m = m_w
        final_eta = eta_guess
        break

    eta_new = Q_(eta_indirect, "")  # make it a Quantity, dimensionless

    # Save current values as "last successful"
    final_m = m_w
    final_eta = eta_new

    # 5) Check convergence in water mass flow
    if prev_m is not None:
        dm = (m_w - prev_m).to("kg/s")
        if abs(dm).magnitude < tol_m.magnitude:
            log.info(f"Converged in {it+1} iterations.")
            break

    # 6) Prepare next iteration
    prev_m = m_w
    eta_guess = eta_new

else:
    # Only executed if loop did not break normally
    log.warning("Did not reach mass-flow convergence within max_iter.")

# Final call: write exactly one set of CSVs with converged flow
final_result = run_hx(
    stages_raw=stages,
    water=water_iter,  # from last iteration
    gas=combustion_results.flue,
    drum=drum,
    target_dx="0.1 m",
    combustion=combustion_results,
    write_csv=True,
)

result = final_result

steps_csv, summary_csv = write_results_csvs(
    global_profile=result["global_profile"],
    combustion=result["combustion"],
    outdir=result["outdir"],
    run_id=result["run_id"],
)

log.info(f"Final water mass flow: {final_m}")
log.info(f"Final indirect efficiency guess: {final_eta}")


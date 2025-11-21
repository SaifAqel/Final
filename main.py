# main.py
import logging

from common.logging_utils import setup_logging
from common.units import Q_
from boiler_loop import run_boiler_case


def run_default_case() -> None:
    run_boiler_case()  # uses all config YAMLs and defaults


def run_excess_air_sensitivity() -> None:
    """
    Example sensitivity analysis on excess air ratio.
    Results are written to separate outdirs/run_ids by run_hx/write_results_csvs.
    """
    # Define a few excess air ratios to test
    ea_values = [1.00, 1.05, 1.10, 1.15, 1.20]

    for ea in ea_values:
        logging.getLogger(__name__).info(f"Running case with excess_air_ratio={ea}")
        run_boiler_case(
            operation_overrides={"excess_air_ratio": Q_(ea, "")},
            # You can also pass different tolerances/iteration controls if desired:
            eta_guess=Q_(0.90, ""),
            tol_m=Q_(1e-3, "kg/s"),
            max_iter=20,
            write_csv=True,
        )


def main() -> None:
    setup_logging("INFO")

    run_default_case()

    # run_excess_air_sensitivity()


if __name__ == "__main__":
    main()

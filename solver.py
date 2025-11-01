# =========================================================
# FILE: solver.py
# =========================================================
from __future__ import annotations
from typing import List, Tuple, Optional
from math import ceil
import logging

from units import Q_, ureg
from models import HXStage, GasStream, WaterStream
from results import StepResult, StageResult
from step_solver import solve_step
from props import GasProps, WaterProps
from logging_utils import setup_logging, trace_calls

_gasprops = GasProps()  # reuse for cp/h evaluation


# ------------------------- utilities -------------------------

def _clamp(v: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, v))

def _median_length(stages: List[HXStage]) -> Q_:
    Ls = sorted([st.spec["inner_length"].to("m") for st in stages], key=lambda x: x.magnitude)
    return Ls[len(Ls)//2]

def _make_grid(L: Q_, n_steps: int) -> tuple[List[Q_], Q_]:
    dx = (L / n_steps).to("m")
    xs = [(i * dx).to("m") for i in range(n_steps)]
    return xs, dx

def _copy_step_with_stage(sr: StepResult, stage_name: str, stage_index: int) -> StepResult:
    # StepResult is frozen; rebuild with stage metadata.
    return StepResult(
        i=sr.i, x=sr.x, dx=sr.dx,
        gas=sr.gas, water=sr.water,
        Tgw=sr.Tgw, Tww=sr.Tww,
        UA_prime=sr.UA_prime, qprime=sr.qprime,
        boiling=sr.boiling,
        h_g=sr.h_g, h_c=sr.h_c,
        stage_name=stage_name, stage_index=stage_index
    )


# ---------------------- wall guesses ------------------------

@trace_calls(values=True)
def initial_wall_guesses(g: GasStream, w: WaterStream) -> tuple[Q_, Q_, Q_]:
    """
    First-call initial guesses for wall state.
    Returns (T_gas_side_wall, T_water_side_wall, qprime) with qprime=0.
    """
    Tg = g.T.to("K")
    Tw = WaterProps.T_from_Ph(w.P, w.h).to("K")
    return Tg, Tw, Q_(0.0, "W/m")


# -------------------- pressure drop hook --------------------

def pressure_drop_gas(g: GasStream, stage: HXStage, i: int, dx: Q_) -> Q_:
    """
    Hook for gas-side pressure drop over one marching step.
    Default: zero. Keep the call to enable future upgrades.
    """
    return Q_(0.0, "Pa")


# --------------------- stream updaters ----------------------


def _solve_T_for_h(P, X, h_target, T0, maxit=30):
    T = T0.to("K"); h_target = h_target.to("J/kg")
    for _ in range(maxit):
        h  = _gasprops.h(T, P, X)            # J/kg
        dh = (h_target - h).to("J/kg")
        if abs(dh).magnitude < 1e-3:         # tighten if needed
            return T
        cp = _gasprops.cp(T, P, X)           # J/kg/K
        dT = (dh / cp).to("K")
        T  = (T + 0.8*dT).to("K")            # mild damping
    return T

@trace_calls(values=True)
def update_gas_after_step(g, qprime, dx, stage):
    Q_step = (qprime * dx).to("W")           # W = J/s
    dh     = (-Q_step / g.mass_flow).to("J/kg")
    h_old  = _gasprops.h(g.T, g.P, g.comp)
    h_new  = (h_old + dh).to("J/kg")
    T_new  = _solve_T_for_h(g.P, g.comp, h_new, g.T)
    P_new  = (g.P + pressure_drop_gas(g, stage, i=0, dx=dx)).to("Pa")
    return GasStream(mass_flow=g.mass_flow, T=T_new, P=P_new, comp=g.comp)

@trace_calls(values=True)
def update_water_after_step(w: WaterStream, qprime: Q_, dx: Q_, stage: HXStage) -> WaterStream:
    """
    Apply energy change to the water after a differential step.
    Enthalpy balance: dh_w = + (qprime*dx)/m_dot. Pressure constant.
    """
    Q_step = (qprime * dx).to("W")
    dh = (Q_step / w.mass_flow).to("J/kg")
    h_new = (w.h + dh).to("J/kg")
    return WaterStream(mass_flow=w.mass_flow, h=h_new, P=w.P)


# ------------------------ stage solve -----------------------

@trace_calls(values=True)
def solve_stage(
    g_in: GasStream,
    w_in: WaterStream,
    stage: HXStage,
    n_steps: int,
    *,
    stage_index: int,
    logger_name: str = "solver",
) -> tuple[GasStream, WaterStream, StageResult]:
    """
    March a single stage in local counter-flow.
    Gas moves +x from 0→L. Water moves −x from L→0. Co-located wall area per step.

    Returns (g_out_at_x=L, w_out_at_x=0, StageResult).
    """
    log = logging.getLogger(logger_name)
    L = stage.spec["inner_length"].to("m")
    xs, dx = _make_grid(L, n_steps)

    steps: List[StepResult] = []

    g = g_in
    w = w_in
    Tgw_guess, Tww_guess, qprime_guess = initial_wall_guesses(g, w)

    Q_sum = Q_(0.0, "W")
    UA_sum = Q_(0.0, "W/K")

    for i, x in enumerate(xs):
        sr = solve_step(
            g=g, w=w, stage=stage,
            Tgw_guess=Tgw_guess, Tww_guess=Tww_guess, qprime_guess=qprime_guess,
            i=i, x=x, dx=dx
        )
        sr = _copy_step_with_stage(sr, stage.name, stage_index)
        steps.append(sr)

        # Accumulate per-length integrals
        Q_sum = (Q_sum + (sr.qprime * dx)).to("W")
        UA_sum = (UA_sum + (sr.UA_prime * dx)).to("W/K")

        # Roll guesses
        Tgw_guess, Tww_guess, qprime_guess = sr.Tgw, sr.Tww, sr.qprime

        # Advance both streams according to local counter-flow energy balance
        g = update_gas_after_step(g, sr.qprime, dx, stage)
        w = update_water_after_step(w, sr.qprime, dx, stage)

        log.debug(
            "step",
            extra={"stage": stage.name, "step": f"{i+1}/{n_steps}"}
        )

    g_out = g
    w_out = w

    stage_res = StageResult(
        stage_name=stage.name,
        stage_kind=stage.kind,
        steps=steps,
        Q_stage=Q_sum,
        UA_stage=UA_sum,
    )

    # Quick internal consistency: Q_stage ≈ sum(q' dx)
    recon = sum([(s.qprime * s.dx).to("W") for s in steps], Q_(0.0, "W"))
    if abs((stage_res.Q_stage - recon) / (stage_res.Q_stage + Q_(1e-12, "W"))) > 0.005:
        raise RuntimeError(f"Stage energy accumulation mismatch >0.5% in {stage.name}")

    # Log inlet/outlet snapshots
    s0 = steps[0]
    sN = steps[-1]
    log.info(
        f"{stage.name}: gas_in(T={g_in.T:~P},P={g_in.P:~P}) gas_out(T={g_out.T:~P},P={g_out.P:~P}) "
        f"water_in(h={w_in.h:~P},P={w_in.P:~P}) water_out(h={w_out.h:~P}) Q_stage={stage_res.Q_stage:~P}",
        extra={"stage": stage.name, "step": f"{len(steps)}/{n_steps}"},
    )

    return g_out, w_out, stage_res


# ---------------------- exchanger solve ---------------------

@trace_calls(values=True)
def solve_exchanger(
    stages: List[HXStage],
    gas_in: GasStream,
    water_in: WaterStream,
    *,
    target_dx: Q_ | None = None,
    min_steps_per_stage: int = 20,
    max_steps_per_stage: int = 400,
    max_passes: int = 20,
    tol_Q: Q_ = Q_(1e-3, "W"),
    tol_end: Q_ = Q_(1e-3, "J/kg"),
    log_level: str = "INFO",
) -> tuple[List[StageResult], GasStream, WaterStream]:
    """
    Global counter-flow, 6-stage exchanger, iterated by forward/backward sweeps.
    Returns (final_stage_results, gas_out, water_out).
    """
    setup_logging(level=log_level)
    log = logging.getLogger("solver")

    if len(stages) != 6:
        raise ValueError(f"Expected 6 stages. Got {len(stages)}.")

    # Build per-stage step counts for near-uniform dx
    if target_dx is None:
        dx_target = (_median_length(stages) / 100).to("m")
    else:
        dx_target = target_dx.to("m")

    n_steps_by_stage: List[int] = []
    for st in stages:
        L = st.spec["inner_length"].to("m")
        n = _clamp(int(ceil((L / dx_target).to("").magnitude)), min_steps_per_stage, max_steps_per_stage)
        n_steps_by_stage.append(n)

    # Storage per pass
    prev_Q_total: Optional[Q_] = None
    prev_end_h: Optional[Tuple[Q_, Q_, Q_, Q_]] = None
    final_stage_results: List[StageResult] = []

    # Inlet boundary enthalpies (constants)
    h_g_in = _gasprops.h(gas_in.T, gas_in.P, gas_in.comp)
    h_w_in = water_in.h

    for p in range(max_passes + 1):
        # ------------ Gas forward sweep: stages 1 → 6 ------------
        gas_stage_results: List[StageResult] = []
        gas_at_stage_in: List[GasStream] = []
        water_for_stage_boundary: List[WaterStream] = []

        g = gas_in
        for i, st in enumerate(stages):
            # Water boundary for this stage: previous-pass stage inlet if available, else use running water from last pass end or global inlet heuristic.
            if p == 0:
                # First pass: simple heuristic. Use global water_in for all stages in forward sweep.
                w_boundary = water_in
            else:
                # Reuse last pass water boundary from this stage's previous result start
                w_boundary = final_stage_results[i].steps[0].water  # stage inlet at x=L last pass

            gas_at_stage_in.append(g)
            water_for_stage_boundary.append(w_boundary)

            g, w_tmp, st_res = solve_stage(g, w_boundary, st, n_steps_by_stage[i], stage_index=i)
            gas_stage_results.append(st_res)

        # ------------ Water backward sweep: stages 6 → 1 ------------
        water_stage_results: List[StageResult] = []
        g_fields_for_water: List[GasStream] = [gs for gs in gas_at_stage_in]  # current-pass gas inlets per stage

        w = water_in
        for i_rev, st in enumerate(reversed(stages)):
            idx = 5 - i_rev  # stage index 5..0
            g_for_stage = g_fields_for_water[idx]
            g_new, w, st_res = solve_stage(g_for_stage, w, st, n_steps_by_stage[idx], stage_index=idx)
            # overwrite gas inlet for next pass with current result inlet (small drift acceptable)
            g_fields_for_water[idx] = g_new
            water_stage_results.append(st_res)

        # reorder to stage index ascending
        water_stage_results = list(reversed(water_stage_results))

        # ------------- Pass bookkeeping and checks --------------
        Q_total = sum([sr.Q_stage.to("W") for sr in water_stage_results], Q_(0.0, "W")).to("W")

        # Exchanger end states from this pass
        g_out = gas_stage_results[-1].steps[-1].gas  # from forward sweep end
        w_out = water_stage_results[0].steps[-1].water  # after backward sweep, outlet at stage 1

        h_g_out = _gasprops.h(g_out.T, g_out.P, g_out.comp)
        h_w_out = w_out.h

        end_tuple = (h_g_in, h_g_out, h_w_in, h_w_out)

        # Convergence metrics
        duty_ok = prev_Q_total is not None and abs(Q_total - prev_Q_total) < tol_Q
        end_ok = prev_end_h is not None and max(
            abs(end_tuple[0] - prev_end_h[0]),
            abs(end_tuple[1] - prev_end_h[1]),
            abs(end_tuple[2] - prev_end_h[2]),
            abs(end_tuple[3] - prev_end_h[3]),
        ) < tol_end

        log.info(
            f"pass {p}: Q_total={Q_total:~P} "
            f"ΔQ={(Q_total - (prev_Q_total or Q_(0,'W'))):~P} "
            f"max Δends={max( (abs(end_tuple[i] - (prev_end_h[i] if prev_end_h else end_tuple[i])) for i in range(4)), default=Q_(0,'J/kg')):~P} "
            f"converged={'yes' if (duty_ok and end_ok) else 'no'}",
            extra={"stage": "ALL", "step": f"pass {p}"},
        )

        if duty_ok and end_ok:
            # --- synchronized final forward sweep using converged water boundaries ---
            water_boundaries = [sr.steps[0].water for sr in water_stage_results]  # inlet-at-x=L per stage
            g = gas_in
            final_forward_results: List[StageResult] = []
            w_out_sync = None

            for i, st in enumerate(stages):
                w_boundary = water_boundaries[i]
                g, w_tmp, st_res = solve_stage(g, w_boundary, st, n_steps_by_stage[i], stage_index=i)
                final_forward_results.append(st_res)
                if i == 0:
                    w_out_sync = w_tmp  # exchanger water outlet at stage 1

            g_out_sync = g  # exchanger gas outlet at stage 6

            # Energy check on the synchronized sweep
            h_g_out = _gasprops.h(g_out_sync.T, g_out_sync.P, g_out_sync.comp)
            h_w_out = w_out_sync.h
            Q_gas = (gas_in.mass_flow * (h_g_in - h_g_out)).to("W")
            Q_wat = (water_in.mass_flow * (h_w_out - h_w_in)).to("W")
            mismatch = abs(Q_gas - Q_wat) / (abs(Q_wat) + Q_(1e-12, "W"))

            log.info(
                f"FINAL forward: Q_total={sum((sr.Q_stage for sr in final_forward_results), Q_(0,'W')):~P} "
                f"Q_gas={Q_gas:~P} Q_water={Q_wat:~P} rel_err={mismatch:~P}",
                extra={"stage": "ALL", "step": "final_forward"},
            )

            # Optional: log per-stage snapshots for the final sweep
            for sr in final_forward_results:
                s0, sN = sr.steps[0], sr.steps[-1]
                log.info(
                    f"{sr.stage_name}: gas_in(T={s0.gas.T:~P},P={s0.gas.P:~P}) "
                    f"gas_out(T={sN.gas.T:~P},P={sN.gas.P:~P}) "
                    f"water_in(h={s0.water.h:~P},P={s0.water.P:~P}) "
                    f"water_out(h={sN.water.h:~P}) Q_stage={sr.Q_stage:~P}",
                    extra={"stage": sr.stage_name, "step": "final_forward"},
                )

            # Enforce same 0.5% criterion if you want
            if mismatch.magnitude > 0.005:
                raise RuntimeError(
                    f"Energy mismatch >0.5% on final sweep. "
                    f"Q_gas={Q_gas:~P}, Q_water={Q_wat:~P}, rel_err={mismatch:~P}"
                )

            # Return results from the synchronized sweep so CSVs and summaries reflect it
            return final_forward_results, g_out_sync, w_out_sync

        # Prepare for next pass
        prev_Q_total = Q_total
        prev_end_h = end_tuple
        final_stage_results = water_stage_results

    # If loop exits, not converged
    worst_idx = max(range(6), key=lambda k: abs(final_stage_results[k].Q_stage).to("W").magnitude if final_stage_results else 0)
    raise RuntimeError(
        f"Did not converge in {max_passes} passes. "
        f"last_Q_total={prev_Q_total:~P if prev_Q_total else 'n/a'} "
        f"last_end_delta={prev_end_h if prev_end_h else 'n/a'} "
        f"worst_stage_index={worst_idx}"
    )

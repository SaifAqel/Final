from common.units import ureg, Q_
from common.models import GasStream
from scipy.optimize import root_scalar
from common.props import GasProps

_gas = GasProps()  # reuse for all rows

def solve_brentq(func, bracket, args, xtol):
    sol = root_scalar(func, bracket=bracket, method="brentq", xtol=xtol, args=args)
    return sol.root

def _flue_sensible_enthalpy(T_in: Q_, T_out: Q_, P: Q_, X: dict, m_dot: Q_) -> Q_:
    dh = _gas.h_sensible(T_out, P, X, Tref=T_in)     # J/kg
    return (m_dot * dh).to(ureg.watt)

def _residual(T_out: Q_, T_in: Q_, P: Q_, X: dict, m_dot: Q_, Q_in: Q_) -> float:
    h_flue = _flue_sensible_enthalpy(T_in, T_out, P, X, m_dot)
    return (Q_in.to(ureg.watt) - h_flue).m

def adiabatic_flame_T(flue: GasStream, Q_in: Q_) -> Q_:
    T_in = flue.T
    P = flue.P
    X = flue.comp
    m_dot = flue.mass_flow
    f = lambda T: _residual(Q_(T, ureg.kelvin), T_in, P, X, m_dot, Q_in)
    T_val = solve_brentq(f, bracket=(1500.0, 3000.0), args=(), xtol=1e-6)
    return T_val * ureg.kelvin

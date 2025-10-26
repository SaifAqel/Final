from __future__ import annotations
from typing import Dict, Optional
from units import Q_

import cantera as ct
from iapws import IAPWS97

class GasProps:
    """Cantera-backed properties for an ideal-gas mixture."""
    def __init__(self, mech_path: str = "config/flue_cantera.yaml", phase: str = "gas_mix"):
        self._sol = ct.Solution(mech_path, phase)

    def _set(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None):
        T_K = (film_T or T).to("K").magnitude
        P_Pa = P.to("Pa").magnitude
        X_map = {k: v.to("dimensionless").magnitude for k, v in X.items()}
        self._sol.TPX = T_K, P_Pa, X_map
        return self._sol

    def cp(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None) -> Q_:
        return Q_(self._set(T,P,X,film_T).cp_mass, "J/kg/K")

    def k(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None) -> Q_:
        return Q_(self._set(T,P,X,film_T).thermal_conductivity, "W/m/K")

    def mu(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None) -> Q_:
        return Q_(self._set(T,P,X,film_T).viscosity, "Pa*s")

    def rho(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None) -> Q_:
        return Q_(self._set(T,P,X,film_T).density, "kg/m^3")

    def h(self, T: Q_, P: Q_, X: Dict[str, Q_], film_T: Optional[Q_] = None) -> Q_:
        return Q_(self._set(T,P,X,film_T).enthalpy_mass, "J/kg")


class WaterProps:
    """IAPWS-97 water/steam properties using (P,h) or (P,T)."""
    @staticmethod
    def _Ph(P: Q_, h: Q_) -> IAPWS97:
        return IAPWS97(P=P.to("megapascal").magnitude, h=h.to("kJ/kg").magnitude)

    @staticmethod
    def _PT(P: Q_, T: Q_) -> IAPWS97:
        return IAPWS97(P=P.to("megapascal").magnitude, T=T.to("K").magnitude)

    # from (P,h)
    @staticmethod
    def T_from_Ph(P: Q_, h: Q_) -> Q_: return Q_(WaterProps._Ph(P,h).T, "K")
    @staticmethod
    def rho_from_Ph(P: Q_, h: Q_) -> Q_: return Q_(WaterProps._Ph(P,h).rho, "kg/m^3")
    @staticmethod
    def mu_from_Ph(P: Q_, h: Q_) -> Q_:  return Q_(WaterProps._Ph(P,h).mu, "Pa*s")
    @staticmethod
    def k_from_Ph(P: Q_, h: Q_) -> Q_:   return Q_(WaterProps._Ph(P,h).k, "W/m/K")
    @staticmethod
    def cp_from_Ph(P: Q_, h: Q_) -> Q_:  return Q_(WaterProps._Ph(P,h).cp, "kJ/kg/K").to("J/kg/K")

    @staticmethod
    def quality_from_Ph(P: Q_, h: Q_) -> Q_ | None:
        Pcrit = Q_(22.064, "MPa")
        if P >= Pcrit:
            return None
        hf = WaterProps.h_f(P)     # J/kg
        hg = WaterProps.h_g(P)     # J/kg
        dh = hg - hf
        if abs(dh.to("J/kg").magnitude) < 1e-9:
            return None

        x = ((h - hf) / dh).to("")  # dimensionless
        xm = x.magnitude
        if xm < -1e-6 or xm > 1 + 1e-6:
            return None

        # clamp to [0,1]
        xm = min(1.0, max(0.0, xm))
        return Q_(xm, "")

    # saturation
    @staticmethod
    def Tsat(P: Q_) -> Q_: return Q_(IAPWS97(P=P.to("megapascal").magnitude, x=0.0).T, "K")
    @staticmethod
    def h_f(P: Q_) -> Q_:  return Q_(IAPWS97(P=P.to("megapascal").magnitude, x=0.0).h, "kJ/kg").to("J/kg")
    @staticmethod
    def h_g(P: Q_) -> Q_:  return Q_(IAPWS97(P=P.to("megapascal").magnitude, x=1.0).h, "kJ/kg").to("J/kg")

    @staticmethod
    def cp_from_PT(P: Q_, T: Q_) -> Q_:
        return Q_(WaterProps._PT(P,T).cp, "kJ/kg/K").to("J/kg/K")
    
    # add below existing cp_from_PT
    @staticmethod
    def mu_from_PT(P: Q_, T: Q_) -> Q_:  return Q_(WaterProps._PT(P,T).mu, "Pa*s")

    @staticmethod
    def k_from_PT(P: Q_, T: Q_) -> Q_:   return Q_(WaterProps._PT(P,T).k, "W/m/K")

    @staticmethod
    def rho_from_PT(P: Q_, T: Q_) -> Q_: return Q_(WaterProps._PT(P,T).rho, "kg/m^3")

    @staticmethod
    def rho_from_Px(P: Q_, x: Q_) -> Q_:
        P_MPa = P.to("megapascal").magnitude
        if x.magnitude <= 0:
            return Q_(IAPWS97(P=P_MPa, x=0.0).rho, "kg/m^3")
        if x.magnitude >= 1:
            return Q_(IAPWS97(P=P_MPa, x=1.0).rho, "kg/m^3")

        rho_f = IAPWS97(P=P_MPa, x=0.0).rho
        rho_g = IAPWS97(P=P_MPa, x=1.0).rho

        # mixture specific volume rule: 1/rho = (1-x)/rho_f + x/rho_g
        v_mix = (1 - x.magnitude) / rho_f + x.magnitude / rho_g
        return Q_(1 / v_mix, "kg/m^3")

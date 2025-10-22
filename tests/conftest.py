import pytest
from units import Q_
import physics as phys

class _FakeGas:
    # constant-property idealized gas to avoid Cantera in tests
    def cp(self, T, P, X, film_T=None): return Q_(1000.0, "J/kg/K")
    def k(self,  T, P, X, film_T=None): return Q_(0.06,   "W/m/K")
    def mu(self, T, P, X, film_T=None): return Q_(3.0e-5, "Pa*s")
    def rho(self, T, P, X, film_T=None):return Q_(0.7,    "kg/m^3")
    def h(self, T, P, X, film_T=None):  return Q_(0.0,    "J/kg")  # not used

@pytest.fixture
def fake_gas(monkeypatch):
    monkeypatch.setattr(phys, "_gas", _FakeGas())

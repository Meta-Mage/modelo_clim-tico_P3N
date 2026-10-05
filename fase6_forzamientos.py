# fase6_forzamientos.py -- estados iniciales y forzamientos de prueba para el nucleo seco (Fase 6.2).
#
#   - Estado de Jablonowski y Williamson (2006) ✅ (formulas cotejadas con el codigo de MPAS,
#     DISENO_FASE6_2.md §6 prueba 4): equilibrio estacionario + perturbacion opcional.
#   - Forzamiento de Held y Suarez (1994) ✅² (DISENO_FASE6_2.md §5).

import math
import numpy as np

# ---------------------------------------------------------------- Jablonowski y Williamson (2006)
JW_U0, JW_ETA0, JW_T0, JW_GAMMA, JW_DT, JW_ETAT, JW_P0 = 35.0, 0.252, 288.0, 0.005, 4.8e5, 0.2, 1.0e5


def sigma_capas(nucleo):
    """sigma en el centro de cada capa segun SB81: ln sigma_k = ln sigma_{k+1/2} - alfa_k."""
    return nucleo.sh[1:] * np.exp(-nucleo.alfa)


def _jw_corchetes(phi, a, omega, u0, cv32):
    A = (-2 * np.sin(phi) ** 6 * (np.cos(phi) ** 2 + 1.0 / 3.0) + 10.0 / 63.0) * u0 * cv32
    B = (8.0 / 5.0 * np.cos(phi) ** 3 * (np.sin(phi) ** 2 + 2.0 / 3.0) - math.pi / 4.0) * a * omega
    return A, B


def estado_jw(nucleo, perturbar=False):
    """Devuelve (ps, T, u, v, phis) del estado de JW06 en la rejilla y las capas del nucleo.
    Coordenada: eta = sigma (p_s = 1000 hPa uniforme)."""
    m = nucleo
    Rj, R, g, a, om = m.Rj, m.R, m.g, m.Rj.a, m.omega
    u0 = JW_U0
    eta = sigma_capas(m)[:, None, None]
    etav = (eta - JW_ETA0) * math.pi / 2
    LAM, PHI = Rj.malla("h")
    Tm = JW_T0 * eta ** (R * JW_GAMMA / g) + np.where(eta < JW_ETAT, JW_DT * (JW_ETAT - eta) ** 5, 0.0)
    cv = np.cos(etav)
    A, B = _jw_corchetes(PHI[None], a, om, 2 * u0, cv ** 1.5)
    T = Tm + 0.75 * eta * math.pi * u0 / R * np.sin(etav) * np.sqrt(cv) * (A + B)
    evs = (1 - JW_ETA0) * math.pi / 2
    As, Bs = _jw_corchetes(PHI, a, om, u0, math.cos(evs) ** 1.5)
    phis = u0 * math.cos(evs) ** 1.5 * (As + Bs)
    LU, PU = Rj.malla("u")
    u = u0 * cv ** 1.5 * np.sin(2 * PU[None]) ** 2
    if perturbar:
        r = np.arccos(np.clip(math.sin(math.radians(40)) * np.sin(PU) + math.cos(math.radians(40)) * np.cos(PU)
                              * np.cos(LU - math.radians(20)), -1, 1))
        u = u + 1.0 * np.exp(-(r / 0.1) ** 2)[None]
    v = np.zeros((m.N, Rj.filas - 1, Rj.columnas))
    ps = np.full((Rj.filas, Rj.columnas), JW_P0)
    return ps, T * np.ones((1, 1, Rj.columnas)), u, v, phis


# ---------------------------------------------------------------- Held y Suarez (1994)
HS = dict(T0=315.0, T_estrat=200.0, dTy=60.0, dthz=10.0, sigma_b=0.7,
          kf=1.0 / 86400, ka=1.0 / (40 * 86400), ks=1.0 / (4 * 86400), p0=1.0e5)


class HeldSuarez:
    def __init__(self, nucleo):
        m = self.m = nucleo
        Rj = m.Rj
        self.sig = sigma_capas(m)[:, None, None]
        _, PHI = Rj.malla("h")
        self.phi = PHI[None]
        w = np.maximum(0.0, (self.sig - HS["sigma_b"]) / (1 - HS["sigma_b"]))
        self.kT = HS["ka"] + (HS["ks"] - HS["ka"]) * w * np.cos(self.phi) ** 4
        self.kv = HS["kf"] * w                                                     # (N,1,1)

    def T_equilibrio(self, ps):
        p = self.sig * ps[None]
        x = p / HS["p0"]
        Teq = (HS["T0"] - HS["dTy"] * np.sin(self.phi) ** 2 - HS["dthz"] * np.log(x) * np.cos(self.phi) ** 2) * x ** self.m.kappa
        return np.maximum(HS["T_estrat"], Teq)

    def __call__(self, ps, T, u, v):
        dT = -self.kT * (T - self.T_equilibrio(ps))
        return np.zeros_like(ps), dT, -self.kv * u, -self.kv * v

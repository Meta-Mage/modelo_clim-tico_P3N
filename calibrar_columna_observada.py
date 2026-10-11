# calibrar_columna_observada.py -- 3.18.0 (DISENO_FASE_6.3.md §6.20.2): calibracion del infrarrojo con la columna
# media OBSERVADA de la Tierra.
#
# Uso:  pip install pyrtlib   (solo para la parte de A y B: trae los perfiles AFGL)
#       python calibrar_columna_observada.py
#
# Las calibraciones de la v3.0 (tau fijo: tau0 y f; DISENO_FASE_2.3.md §2) y de la v3.1 (vapor radiativo, I15: A y B;
# calibrar_tau_vapor.py) usaban una columna con la superficie a 1000 hPa como "columna media de la Tierra". La
# observada tiene 983,05 hPa de aire seco (Trenberth y Smith 2005; M3N cuenta el aire seco: nucleo de masa seca).
# Este script:
#   1. REPRODUCE las dos calibraciones originales con la columna de 1000 hPa (comprobacion del metodo);
#   2. las repite con la columna observada (p_s = P_COLUMNA_OBSERVADA, normalizacion p/P_COLUMNA_OBSERVADA).
# La columna observada es LA MISMA columna en terminos relativos: su temperatura es la de la US 1976 en p*P_REF/p_s
# (la misma en cada fraccion de la columna) y el vapor tiene la misma forma y los mismos 24,9 kg/m2. Asi solo cambia
# su masa absoluta, que es lo que se corrige: la columna de 1000 hPa hacia de "columna media de la Tierra".
# Mismo metodo y mismos objetivos: Atmosfera Estandar US 1976, suelo a 289 K, flujos sin nubes de Wild et al. (2019),
# DLR = 314 y OLR = 267 W/m2; tau0 y f con N = 160 capas (convergido, DISENO_FASE_2.3.md §2); A y B con la columna
# de calibrar_tau_vapor.py (N por defecto) y el vapor escalado a 24,9 kg/m2.
# No se ajusta nada a resultados del modelo: solo cambia la columna de referencia.

import numpy as np
from scipy.optimize import fsolve

import fase30_multicapa as M30

SB = M30.CONSTANTE_SB
G = 9.80665
MW, MA = 18.015, 28.964
AGUA_PRECIPITABLE_TIERRA = 24.9


def std1976(p):
    zz = np.arange(0, 80000, 5.0)
    T = np.where(zz < 11000, 288.15 - 0.0065 * zz, np.where(zz < 20000, 216.65,
                 np.where(zz < 32000, 216.65 + 0.001 * (zz - 20000), 228.65 + 0.0028 * (zz - 32000))))
    lnp = np.log(101325.0) - np.cumsum(G / (M30.R_AIRE * T)) * 5.0
    return np.interp(-np.log(np.maximum(p, 1e-3)), -lnp, T)


def flujos(C, T, Ts):
    D, Bh = C.infrarrojo_bajada(T, np.full(C.forma, Ts))
    U = C.infrarrojo_subida(Bh, np.full(C.forma, SB * Ts ** 4))
    return float(U[0].ravel()[0]), float(D[-1].ravel()[0])


def columna(observada, n):
    C = M30.Columna(np.zeros((1, 1)), gravedad=G, n=n, p_nivel_mar=M30.P_COLUMNA_OBSERVADA if observada else M30.P_REF)
    if observada:                     # la normalizacion de la columna observada (tau0 y f se calibran despues)
        C.p_tau = M30.P_COLUMNA_OBSERVADA
        C.actualizar_ps(C.ps)
    return C


def temperatura(C):
    """US 1976 en cada fraccion de la columna (p*P_REF/p_s): en la de 1000 hPa, la de la v3.0 y la v3.1."""
    return std1976(C.pm[:, 0, 0] * M30.P_REF / float(C.ps.ravel()[0]))[:, None, None]


def calibrar_tau_fijo(observada, n=160):
    """tau0 y f del perfil tau = tau0*(f*x + (1-f)*x^4), x = p/p_tau."""
    C = columna(observada, n)
    T = temperatura(C)

    def ec(x):
        C.tau0, C.f_lw = x
        C.actualizar_ps(C.ps)
        o, d = flujos(C, T, 289.0)
        return [d - 314.0, o - 267.0]
    return fsolve(ec, [2.45, 0.24], xtol=1e-12)


def calibrar_vapor(observada, n=M30.N_CAPAS_ATM):
    """A y B de d(tau) = (A + B*q)*dp/p_tau con el perfil de vapor US de la AFGL escalado a 24,9 kg/m2."""
    from pyrtlib.climatology import AtmosphericProfiles as atmp
    C = columna(observada, n)
    pm = C.pm[:, 0, 0]
    z, p, d, t, md = atmp.gl_atm(atmp.US_STANDARD)
    x = md[:, atmp.H2O] * 1e-6                                   # fraccion molar del vapor
    perfil = x * MW / MA if not observada else x * MW / (MA * (1.0 - x))   # obs.: razon de mezcla (por kg de aire seco)
    escala = M30.P_REF / float(C.ps.ravel()[0])                  # la misma forma en cada fraccion de la columna
    q = np.interp(-np.log(pm * escala), -np.log(p * 100), perfil)[:, None, None]
    q *= AGUA_PRECIPITABLE_TIERRA / float((q * C.dp).sum() / G)
    T = temperatura(C)

    def ec(xx):
        C.actualizar_tau_vapor(q, *xx)
        o, dd = flujos(C, T, 289.0)
        return [dd - 314.0, o - 267.0]
    return fsolve(ec, [0.8, 2000.0], xtol=1e-12)


def main():
    t0, f0 = calibrar_tau_fijo(False)
    print(f"tau fijo, columna de 1000 hPa con el codigo actual: tau0 = {t0:.4f}, f = {f0:.4f}  (los de la v3.0, 2,452 / "
          f"0,240, no se reproducen: DISENO_FASE_6.3.md §6.21.5)")
    t1, f1 = calibrar_tau_fijo(True)
    print(f"tau fijo, columna observada ({M30.P_COLUMNA_OBSERVADA / 100:.2f} hPa):  tau0 = {t1:.4f}, f = {f1:.4f}  "
          f"(en fase30_multicapa: {M30.TAU0_LW_OBS}, {M30.F_LW_OBS})")
    try:
        a0, b0 = calibrar_vapor(False)
        print(f"vapor (I15), columna de 1000 hPa (comprobacion; v3.1: 0,5656 / 773,1): A = {a0:.4f}, B = {b0:.1f}")
        a1, b1 = calibrar_vapor(True)
        print(f"vapor (I15), columna observada:  A = {a1:.4f}, B = {b1:.1f}")
    except ImportError:
        print("Sin pyrtlib: no se calibran A y B (pip install pyrtlib)")


if __name__ == "__main__":
    main()

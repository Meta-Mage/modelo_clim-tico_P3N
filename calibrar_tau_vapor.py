# calibrar_tau_vapor.py -- v3.1 (I15, prototipo): calibracion del infrarrojo con vapor radiativo.
#
# Uso:  pip install pyrtlib   (solo para este script: trae los perfiles AFGL)
#       python calibrar_tau_vapor.py
#
# Forma (Byrne y O'Gorman 2013; codigo de Isca, esquema 'byrne' ✅):
#     d(tau) = (A + B*q) * dp / P0
# A (aire seco) y B (vapor) se calibran igual que el tau fijo de la v3.0 (DISENO_V3.0.md seccion 2):
# una columna media de la Tierra sin nubes debe dar DLR = 314 y OLR = 267 W/m2 (Wild et al. 2019),
# con la temperatura de la Atmosfera Estandar 1976 (suelo 289 K) y un perfil de vapor con la forma del
# de la atmosfera estandar de EE. UU. de la AFGL (Anderson et al. 1986, via pyrtlib) escalado al agua
# precipitable media de la Tierra, 24,9 kg/m2 (Trenberth y Smith 2005: 1,27e16 kg de vapor ✅).
# Despues compara, en las columnas tipo de la AFGL, el esquema fijo de la v3.0 con este.

import numpy as np
from scipy.optimize import fsolve
from pyrtlib.climatology import AtmosphericProfiles as atmp

import fase30_multicapa as M30
import fase31_agua as A

SB = M30.CONSTANTE_SB
P0 = M30.P0
G = 9.80665
MW, MA = 18.015, 28.964
AGUA_PRECIPITABLE_TIERRA = 24.9


def perfil(op):
    z, p, d, t, md = atmp.gl_atm(op)
    return p * 100, t, md[:, atmp.H2O] * 1e-6 * MW / MA


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


def main():
    C = M30.Columna(np.zeros((1, 1)), gravedad=G)
    pm = C.pm[:, 0, 0]
    interp = lambda p, x: np.interp(-np.log(pm), -np.log(p), x)[:, None, None]
    p_us, t_us, q_us = perfil(atmp.US_STANDARD)
    q = interp(p_us, q_us)
    q *= AGUA_PRECIPITABLE_TIERRA / float((q * C.dp).sum() / G)
    T = std1976(pm)[:, None, None]

    def ec(x):
        C.actualizar_tau_vapor(q, *x)
        o, d = flujos(C, T, 289.0)
        return [d - 314.0, o - 267.0]
    a, b = fsolve(ec, [0.8, 2000.0])
    print(f"Calibrado: A = {a:.4f}, B = {b:.1f}   (en fase31_agua: {A.A_LW_SECO}, {A.B_LW_VAPOR})")
    for nombre, op in (("Tropical", atmp.TROPICAL), ("Latitud media, verano", atmp.MIDLATITUDE_SUMMER),
                       ("Estandar EE. UU.", atmp.US_STANDARD), ("Subartico, invierno", atmp.SUBARCTIC_WINTER)):
        p, t, qq = perfil(op)
        Tc, qc, Ts = interp(p, t), interp(p, qq), float(t[0])
        Cf = M30.Columna(np.zeros((1, 1)), gravedad=G)
        o1, d1 = flujos(Cf, Tc, Ts)
        C.actualizar_tau_vapor(qc, a, b)
        o2, d2 = flujos(C, Tc, Ts)
        print(f"{nombre:22s} suelo {Ts:.1f} K | tau fijo: OLR {o1:.0f} DLR {d1:.0f} | con vapor: OLR {o2:.0f} DLR {d2:.0f}")


if __name__ == "__main__":
    main()

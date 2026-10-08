# fase6_superficie.py -- Fase 6.3 (paso 5): intercambio con la superficie con VIENTO REAL (sin conectar todavia).
# Diseno: DISENO_FASE_6.3.md §6.2 (decisiones 1.3 y 1.4) y §6.7.3.
#
# 1. Rugosidad del oceano del ECMWF (IFS Cy31r1/40r1, implementacion de NEMO sbcblk_algo_ecmwf.F90 ✅):
#      z0  = 0,11 nu/u* + 0,018 u*^2/g   (Charnock; 0,018 el valor tipico del ECMWF, que varia con el oleaje)
#      z0h = 0,40 nu/u*,  z0q = 0,62 nu/u*      nu = 1,5e-5 m2/s (ECMWF TM 630 ✅)
#    NEMO acota los tres a <= 1e-3 m; aqui igual.
# 2. Coeficientes de intercambio a la altura real z_a de la capa mas baja (§6.5, hallazgo 1), con la estabilidad
#    de Louis, Tiedtke y Geleyn (1982) para el momento y para el calor (fase6_capa_limite.py ✅). Sobre el oceano
#    z0 depende de u*, que depende de z0: punto fijo, iterado hasta < 1e-12 relativo.
#    EVAPORACION: la v3.1 decidio C_E = C_H (Frierson 2007; Isca). 3.17.0 (decision de Carlos del 09/10/2026,
#    DISENO_FASE_6.3.md §6.18-§6.19): sobre el AGUA, C_E con la z0q del ECMWF (C_E/C_H ~1,04; Large y Yeager 2004
#    dan 1,058); sobre tierra, C_E = C_H (z0q = z0h). Se devuelve como "c_e".
# 3. Albedo directo del oceano de Cox y Munk (1954) con el viento LOCAL: tabla en (mu, U) con la misma funcion
#    de fase2b_atmosfera.py (que hoy usa U = 5 m/s fijo). En U = 5 m/s da exactamente el albedo de la v3.1.

import math
import numpy as np

import fase2b_atmosfera as F2
from fase6_capa_limite import KARMAN, louis_momento, louis_calor

NU_AIRE = 1.5e-5             # m2/s (ECMWF TM 630)
ALFA_M, ALFA_CH = 0.11, 0.018
ALFA_H, ALFA_Q = 0.40, 0.62
Z0_MAX_OCEANO = 1e-3         # m (como NEMO)
U_ESTRELLA_MIN = 1e-4        # m/s: solo para no dividir por cero en nu/u* con aire en calma
VIENTOS_TABLA = np.arange(0.0, 26.0, 1.0)      # m/s; por encima de 25 m/s se usa el de 25


def rugosidad_oceano(u_estrella, g):
    us = np.maximum(u_estrella, U_ESTRELLA_MIN)
    z0 = np.minimum(np.abs(ALFA_M * NU_AIRE / us + ALFA_CH * us * us / g), Z0_MAX_OCEANO)
    z0h = np.minimum(ALFA_H * NU_AIRE / us, Z0_MAX_OCEANO)
    z0q = np.minimum(ALFA_Q * NU_AIRE / us, Z0_MAX_OCEANO)
    return z0, z0h, z0q


def coeficientes(z_a, viento, theta_a, theta_s, es_agua, z0m_tierra, z0h_tierra, g, iteraciones=30):
    """C_m (momento), C_h (calor) y, desde la 3.17.0, C_e (vapor) a la altura z_a, con la estabilidad de Louis
    (la del calor tambien para el vapor). Sobre el agua, z0, z0h y z0q de la rugosidad del ECMWF con el u* que
    resulta (punto fijo); sobre tierra, z0q = z0h, asi que C_e = C_h exactamente. Devuelve un dict."""
    v = np.maximum(viento, 1e-2)                                              # minimo numerico (§6.2, 1.3)
    ri = g * z_a * (theta_a - theta_s) / (0.5 * (theta_a + theta_s) * v * v)
    z0m = np.where(es_agua, 2e-4, z0m_tierra)
    z0h = np.where(es_agua, 2e-5, z0h_tierra)
    for _ in range(iteraciones):
        lm = np.log(z_a / z0m)
        c_n = (KARMAN / lm) ** 2
        c_m = c_n * louis_momento(ri, c_n, z_a / z0m)
        u_est = np.sqrt(c_m) * v
        z0o, z0ho, _ = rugosidad_oceano(u_est, g)
        z0m_n = np.where(es_agua, z0o, z0m_tierra)
        z0h_n = np.where(es_agua, z0ho, z0h_tierra)
        cambio = max(np.max(np.abs(z0m_n / z0m - 1)), np.max(np.abs(z0h_n / z0h - 1)))
        z0m, z0h = z0m_n, z0h_n
        if cambio < 1e-12:
            break
    lm = np.log(z_a / z0m); lh = np.log(z_a / z0h)
    c_n = (KARMAN / lm) ** 2
    c_m = c_n * louis_momento(ri, c_n, z_a / z0m)
    f_h = louis_calor(ri, c_n, z_a / z0m)
    c_h = KARMAN ** 2 / (lm * lh) * f_h
    _, _, z0q_o = rugosidad_oceano(np.sqrt(c_m) * v, g)          # 3.17.0: con el u* final, como z0 y z0h
    z0q = np.where(es_agua, z0q_o, z0h)
    c_e = np.where(es_agua, KARMAN ** 2 / (lm * np.log(z_a / z0q)) * f_h, c_h)
    return {"c_m": c_m, "c_h": c_h, "c_e": c_e, "u_estrella": np.sqrt(c_m) * v, "z0m": z0m, "z0h": z0h,
            "z0q": z0q, "ri": ri, "cambio_final": cambio}


_tabla_2d = None
_RUTA_TABLA = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),
                                        "outputs", "cache", "tabla_cox_munk_2d_v1.npz")


def tabla_albedo_directo_viento():
    """Albedo directo de Cox y Munk en (mu, U): misma rejilla en mu que fase2b_atmosfera.tabla_albedo_directo.
    Construirla tarda ~40 s en el entorno de la IA: se guarda en outputs/cache/ y, al cargarla, se comprueban
    EXACTAMENTE 5 valores recalculados (si no coinciden, se reconstruye)."""
    global _tabla_2d
    if _tabla_2d is not None:
        return _tabla_2d
    import os
    mus = np.linspace(0.005, 1.0, 200)
    if os.path.exists(_RUTA_TABLA):
        d = np.load(_RUTA_TABLA)
        if np.array_equal(d["mus"], mus) and np.array_equal(d["us"], VIENTOS_TABLA):
            t = d["tabla"]
            ok = all(t[k, j] == F2._albedo_directo_cox_munk(mus[j], float(VIENTOS_TABLA[k]))
                     for k, j in ((0, 0), (5, 57), (12, 199), (20, 100), (25, 3)))
            if ok:
                _tabla_2d = (mus, VIENTOS_TABLA, t)
                return _tabla_2d
    t = np.array([[F2._albedo_directo_cox_munk(m, float(u)) for m in mus] for u in VIENTOS_TABLA])
    os.makedirs(os.path.dirname(_RUTA_TABLA), exist_ok=True)
    tmp = _RUTA_TABLA + ".tmp.npz"
    np.savez(tmp, mus=mus, us=VIENTOS_TABLA, tabla=t)
    os.replace(tmp, _RUTA_TABLA)
    _tabla_2d = (mus, VIENTOS_TABLA, t)
    return _tabla_2d


def albedo_oceano_viento(mu, transmitancia, viento):
    """Como fase2b_atmosfera.albedo_oceano, con el viento de cada celda (interpolacion lineal en U entre los
    nodos de 1 m/s; en mu, np.interp como la v3.1). Vectorizado por tramos de viento."""
    mus, us, tabla = tabla_albedo_directo_viento()
    mu = np.asarray(mu, dtype=float)
    U = np.clip(np.broadcast_to(viento, mu.shape), us[0], us[-1])
    i = np.minimum(np.floor(U).astype(int), len(us) - 2)
    w = U - us[i]
    lo = np.empty(mu.shape); hi = np.empty(mu.shape)
    for k in np.unique(i):
        sel = i == k
        lo[sel] = np.interp(mu[sel], mus, tabla[k])
        hi[sel] = np.interp(mu[sel], mus, tabla[k + 1])
    directo = lo + w * (hi - lo)
    directa = transmitancia ** (F2.TAU_DIRECTO / F2.TAU_NUEVO)
    fraccion_difusa = 1 - directa / np.maximum(transmitancia, 1e-300)
    return (1 - fraccion_difusa) * directo + fraccion_difusa * F2.ALBEDO_DIFUSO_AGUA + F2.ALBEDO_SUBSUPERFICIE

# fase6_nucleo_nb.py -- v3.1-pre8: version compilada (numba) de NucleoSeco.tendencias, SIN el filtro polar.
#
# Por que: con campos de ~52 000 valores, la version de numpy gasta casi todo el tiempo en llamar a decenas de
# operaciones pequenas (DISENO_FASE6_3.md §6.6). Aqui se hace todo en un unico recorrido compilado.
#
# Regla de esta traduccion: cada expresion reproduce EXACTAMENTE el orden de operaciones (y los parentesis
# implicitos, de izquierda a derecha) de la version de numpy de fase6_nucleo.py, para que el resultado sea
# IDENTICO BIT A BIT (lo comprueba test_fase6_2.py::test_tendencias_compiladas_identicas). Las sumas sobre
# capas (dps, cumsum) se hacen en el mismo orden secuencial que numpy usa al reducir sobre el eje 0.
# Convenios de indices como fase6_nucleo.py: [k, i, j] = [capa, fila, columna]; vecino oeste j-1 (periodico).

import numpy as np

try:
    from numba import njit
    HAY_NUMBA = True
except ImportError:            # sin numba, NucleoSeco usa la version de numpy (mismo resultado)
    HAY_NUMBA = False


def _tendencias_sin_filtro(ps, lps, T, u, v, phis, dsig, dln, alfa, sh, L_u, L_v, area, dx_u, dy, A_u, A_v,
                           L_n, area_dual, f_esq, R, cp):
    """lps = np.log(ps) calculado FUERA con numpy: el log vectorizado de numpy y el de la biblioteca
    matematica que usa numba pueden diferir en el ultimo bit (medido el 07/10), y asi el resultado es
    identico. Arrays de geometria por FILA (1D): L_v (F-1), area (F), dx_u (F), A_u (F), A_v (F-1),
    L_n (F), area_dual (F-1), f_esq (F-1). L_u y dy, escalares."""
    N, F, C = T.shape
    # ---------------- presion y flujos de masa ----------------
    ps_u = np.empty((F, C)); ps_v = np.empty((F - 1, C))
    for i in range(F):
        for j in range(C):
            ps_u[i, j] = 0.5 * (ps[i, j] + ps[i, j - 1])
    for i in range(F - 1):
        for j in range(C):
            ps_v[i, j] = 0.5 * (ps[i, j] + ps[i + 1, j])
    dp = np.empty((N, F, C)); F_u = np.empty((N, F, C)); F_v = np.empty((N, F - 1, C))
    for k in range(N):
        for i in range(F):
            for j in range(C):
                dp[k, i, j] = dsig[k] * ps[i, j]
                F_u[k, i, j] = (u[k, i, j] * (dsig[k] * ps_u[i, j])) * L_u
        for i in range(F - 1):
            for j in range(C):
                F_v[k, i, j] = (v[k, i, j] * (dsig[k] * ps_v[i, j])) * L_v[i]
    # divergencia de la masa, D = nabla.(dp v)
    D = np.empty((N, F, C))
    for k in range(N):
        for i in range(F):
            for j in range(C):
                jp = j + 1 if j + 1 < C else 0
                x = F_u[k, i, jp] - F_u[k, i, j]
                if i < F - 1:
                    x = x - F_v[k, i, j]
                if i > 0:
                    x = x + F_v[k, i - 1, j]
                D[k, i, j] = x / area[i]
    # dps = -sum_k D (suma secuencial en k, como numpy sobre el eje 0) y cumsum en k
    cum = np.empty((N, F, C))
    dps = np.empty((F, C))
    for i in range(F):
        for j in range(C):
            s = D[0, i, j]
            cum[0, i, j] = s
            for k in range(1, N):
                s = s + D[k, i, j]
                cum[k, i, j] = s
    for i in range(F):
        for j in range(C):
            s = D[0, i, j]
            for k in range(1, N):
                s = s + D[k, i, j]
            dps[i, j] = -s
    W = np.zeros((N + 1, F, C))
    for k in range(1, N):
        for i in range(F):
            for j in range(C):
                W[k, i, j] = (-cum[k - 1, i, j]) - (sh[k] * dps[i, j])
    # ---------------- temperatura ----------------
    T_u = np.empty((N, F, C)); T_v = np.empty((N, F - 1, C))
    for k in range(N):
        for i in range(F):
            for j in range(C):
                T_u[k, i, j] = 0.5 * (T[k, i, j] + T[k, i, j - 1])
        for i in range(F - 1):
            for j in range(C):
                T_v[k, i, j] = 0.5 * (T[k, i, j] + T[k, i + 1, j])
    g_u = np.empty((F, C)); g_v = np.empty((F - 1, C))
    for i in range(F):
        for j in range(C):
            g_u[i, j] = (lps[i, j] - lps[i, j - 1]) / dx_u[i]
    for i in range(F - 1):
        for j in range(C):
            g_v[i, j] = (lps[i, j] - lps[i + 1, j]) / dy
    dT = np.empty((N, F, C))
    fT_u = np.empty((F, C)); fT_v = np.empty((F - 1, C)); Y_u = np.empty((F, C)); Y_v = np.empty((F - 1, C))
    for k in range(N):
        for i in range(F):
            for j in range(C):
                fT_u[i, j] = F_u[k, i, j] * T_u[k, i, j]
                Y_u[i, j] = ((((A_u[i] * (dsig[k] * ps_u[i, j])) * u[k, i, j]) * R) * T_u[k, i, j]) * g_u[i, j]
        for i in range(F - 1):
            for j in range(C):
                fT_v[i, j] = F_v[k, i, j] * T_v[k, i, j]
                Y_v[i, j] = ((((A_v[i] * (dsig[k] * ps_v[i, j])) * v[k, i, j]) * R) * T_v[k, i, j]) * g_v[i, j]
        for i in range(F):
            for j in range(C):
                jp = j + 1 if j + 1 < C else 0
                a = fT_u[i, jp] - fT_u[i, j]
                if i < F - 1:
                    a = a - fT_v[i, j]
                if i > 0:
                    a = a + fT_v[i - 1, j]
                adv_h = a / area[i]
                # flujo vertical W * T en los seminiveles (0 arriba y abajo)
                ts_ab = 0.5 * (T[k, i, j] + T[k + 1, i, j]) if k < N - 1 else 0.0
                ts_ar = 0.5 * (T[k - 1, i, j] + T[k, i, j]) if k > 0 else 0.0
                adv_v = W[k + 1, i, j] * ts_ab - W[k, i, j] * ts_ar
                d_enc = cum[k, i, j] - D[k, i, j]
                conv1 = ((-R) * T[k, i, j]) * ((dln[k] * d_enc) + (alfa[k] * D[k, i, j]))
                x = 0.5 * (Y_u[i, j] + Y_u[i, jp])
                if i < F - 1:
                    x = x + 0.5 * Y_v[i, j]
                if i > 0:
                    x = x + 0.5 * Y_v[i - 1, j]
                conv2 = x / area[i]
                d_dpT = ((-adv_h) - adv_v) + ((conv1 + conv2) / cp)
                d_dp = dsig[k] * dps[i, j]
                dT[k, i, j] = (d_dpT - T[k, i, j] * d_dp) / dp[k, i, j]
    # ---------------- momento ----------------
    # geopotencial (SB81): Phi_k = phis + sum_{j>k} cap_j + alfa_k R T_k, con la suma acumulada desde abajo
    Phi = np.empty((N, F, C))
    for i in range(F):
        for j in range(C):
            s = 0.0
            for k in range(N - 1, -1, -1):
                cap = (R * T[k, i, j]) * dln[k]
                s_k = s + cap if k < N - 1 else cap
                Phi[k, i, j] = (phis[i, j] + (s_k - cap)) + ((alfa[k] * R) * T[k, i, j])
                s = s_k
    du = np.empty((N, F, C)); dv = np.empty((N, F - 1, C))
    q = np.empty((F - 1, C)); qV = np.empty((F - 1, C)); qU = np.empty((F - 1, C)); B = np.empty((F, C))
    for k in range(N):
        for i in range(F - 1):
            for j in range(C):
                dpe = 0.25 * (((dp[k, i, j] + dp[k, i + 1, j]) + dp[k, i, j - 1]) + dp[k, i + 1, j - 1])
                circ = ((dy * (v[k, i, j] - v[k, i, j - 1])) - (u[k, i, j] * L_n[i])) + (u[k, i + 1, j] * L_n[i + 1])
                q[i, j] = (circ / area_dual[i] + f_esq[i]) / dpe
        for i in range(F - 1):
            for j in range(C):
                qV[i, j] = (q[i, j] * 0.5) * (F_v[k, i, j] + F_v[k, i, j - 1])     # numpy: q * 0.5 * (...)
                qU[i, j] = (q[i, j] * 0.5) * (F_u[k, i, j] + F_u[k, i + 1, j])
        for i in range(F):
            for j in range(C):
                jp = j + 1 if j + 1 < C else 0
                ku = (A_u[i] * u[k, i, j]) * u[k, i, j]
                ku_e = (A_u[i] * u[k, i, jp]) * u[k, i, jp]
                kv_n = (A_v[i - 1] * v[k, i - 1, j]) * v[k, i - 1, j] if i > 0 else 0.0
                kv_s = (A_v[i] * v[k, i, j]) * v[k, i, j] if i < F - 1 else 0.0
                K = (0.5 * (ku + ku_e) + 0.5 * (kv_n + kv_s)) / (2 * area[i])
                B[i, j] = K + Phi[k, i, j]
        for i in range(F):
            for j in range(C):
                fu = 0.0
                if i > 0:
                    fu = fu + 0.5 * qV[i - 1, j]
                if i < F - 1:
                    fu = fu + 0.5 * qV[i, j]
                fu = fu / dx_u[i]
                x = (fu - ((B[i, j] - B[i, j - 1]) / dx_u[i])) - ((R * T_u[k, i, j]) * g_u[i, j])
                # adveccion vertical con W en la cara oeste
                dpf = dsig[k] * ps_u[i, j]
                t = 0.0
                if k < N - 1:
                    w_ab = 0.5 * (W[k + 1, i, j] + W[k + 1, i, j - 1])
                    t = t + w_ab * (u[k + 1, i, j] - u[k, i, j])
                if k > 0:
                    w_ar = 0.5 * (W[k, i, j] + W[k, i, j - 1])
                    t = t + w_ar * (u[k, i, j] - u[k - 1, i, j])
                du[k, i, j] = x - t / (2 * dpf)
        for i in range(F - 1):
            for j in range(C):
                jp = j + 1 if j + 1 < C else 0
                fv = ((-0.5) * (qU[i, j] + qU[i, jp])) / dy
                x = (fv - ((B[i, j] - B[i + 1, j]) / dy)) - ((R * T_v[k, i, j]) * g_v[i, j])
                dpf = dsig[k] * ps_v[i, j]
                t = 0.0
                if k < N - 1:
                    w_ab = 0.5 * (W[k + 1, i, j] + W[k + 1, i + 1, j])
                    t = t + w_ab * (v[k + 1, i, j] - v[k, i, j])
                if k > 0:
                    w_ar = 0.5 * (W[k, i, j] + W[k, i + 1, j])
                    t = t + w_ar * (v[k, i, j] - v[k - 1, i, j])
                dv[k, i, j] = x - t / (2 * dpf)
    return dps, dT, du, dv


if HAY_NUMBA:
    _tendencias_sin_filtro = njit(cache=True)(_tendencias_sin_filtro)

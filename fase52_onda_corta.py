# fase52_onda_corta.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §3.4 y §7): la LUZ (onda corta) en una
# capa que absorbe y DISPERSA (nubes, aerosoles, el propio aire), con el metodo de dos flujos delta-Eddington. Es lo
# que le falta a M3N para las nubes: hoy la luz solo se absorbe o se refleja en el suelo. No se usa todavia en el modelo.
#
# Metodo (decision tecnica T7, PROVISIONAL):
#   - Ecuaciones de dos flujos en la forma general de Meador y Weaver (1980, J. Atmos. Sci. 37, 630 ⚠️), con el haz
#     directo de la estrella como fuente (tau crece hacia abajo; S = flujo del haz medido perpendicular a el):
#         dF+/dtau =  g1 F+ - g2 F- - omega S g3 exp(-tau/mu0)
#         dF-/dtau =  g2 F+ - g1 F- + omega S g4 exp(-tau/mu0)
#   - Coeficientes de EDDINGTON (Shettle y Weinman 1970 ⚠️; tabla 1 de Meador y Weaver 1980 ⚠️):
#         g1 = (7 - omega (4 + 3 g)) / 4,  g2 = -(1 - omega (4 - 3 g)) / 4,  g3 = (2 - 3 g mu0) / 4,  g4 = 1 - g3
#   - Escalado DELTA (Joseph, Wiscombe y Weinman 1976, J. Atmos. Sci. 33, 2452 ⚠️): la fraccion f = g^2 de la luz que
#     se dispersa hacia delante se trata como no dispersada:
#         tau' = (1 - omega f) tau,  omega' = (1 - f) omega / (1 - omega f),  g' = (g - f) / (1 - f) = g / (1 + g)
#   - La solucion de la capa se obtiene EXACTA para esas ecuaciones (autovalores +-k, solucion particular del haz y
#     las dos condiciones de contorno: nada difuso entra por arriba; el suelo refleja como Lambert con albedo A),
#     resolviendo un sistema de 2x2, en vez de copiar formulas cerradas: asi lo que se prueba es el sistema.
#   - Se comprueba frente a un calculo de referencia de ordenadas discretas con 16 direcciones (PythonicDISORT 1.8,
#     Ho 2024, licencia MIT ⚠️, implementacion en Python de DISORT, Stamnes et al. 1988), con valores tabulados en la
#     prueba. Lo que se exige es el error conocido del delta-Eddington (unos pocos % en la reflexion; mayor con el sol
#     bajo: King y Harshvardhan 1986, J. Atmos. Sci. 43, 784 ⚠️), no la igualdad.
#
# Varias capas (la atmosfera real) se haran sumando capas ("adding"); este archivo es la pieza de una capa.

import math

import numpy as np

OMEGA_MAX = 1.0 - 1e-9      # con omega = 1 los dos autovalores coinciden (k = 0): se usa este limite (error < 1e-6)


def escalado_delta(tau, omega, g):
    """Escalado delta-Eddington (f = g^2): (tau', omega', g')."""
    f = g * g
    return (1.0 - omega * f) * tau, (1.0 - f) * omega / (1.0 - omega * f), g / (1.0 + g)


def gammas_eddington(omega, g, mu0):
    g1 = (7.0 - omega * (4.0 + 3.0 * g)) / 4.0
    g2 = -(1.0 - omega * (4.0 - 3.0 * g)) / 4.0
    g3 = (2.0 - 3.0 * g * mu0) / 4.0
    return g1, g2, g3, 1.0 - g3


def capa(tau, omega, g, mu0, albedo_suelo=0.0, delta=True):
    """Una capa homogenea sobre un suelo de Lambert, iluminada por un haz de flujo 1 (perpendicular al haz) con coseno
    cenital mu0. Devuelve fracciones de la luz que llega (mu0 por unidad de area horizontal):
      R      reflejada hacia arriba en el tope;
      T_dir  directa que llega al suelo (sin escalar: exp(-tau/mu0));  T_dif  el resto de la que llega (difusa);
      A_capa absorbida en la capa;  A_suelo absorbida por el suelo.
    R + A_capa + A_suelo = 1."""
    tau_real = tau
    if delta:
        tau, omega, g = escalado_delta(tau, omega, g)
    omega = min(float(omega), OMEGA_MAX)
    g1, g2, g3, g4 = gammas_eddington(omega, g, mu0)
    k = math.sqrt(max(g1 * g1 - g2 * g2, 0.0))
    if abs(k * mu0 - 1.0) < 1e-6:            # resonancia de la solucion particular: se aparta mu0 un poco
        mu0 = mu0 * (1.0 + 2e-6)
    S = 1.0
    # solucion particular c exp(-tau/mu0): (-I/mu0 - M) c = b
    M = np.array([[g1, -g2], [g2, -g1]])
    b = np.array([-omega * S * g3, omega * S * g4])
    c = np.linalg.solve(-np.eye(2) / mu0 - M, b)
    vp = np.array([g2, g1 - k])              # autovector de +k (crece hacia abajo): se escribe como exp(k (tau - tau*))
    vm = np.array([g2, g1 + k])              # autovector de -k: exp(-k tau)
    e_k = math.exp(-k * tau)
    e_0 = math.exp(-tau / mu0)
    # F-(0) = 0 ; F+(tau*) - A F-(tau*) = A mu0 S exp(-tau*/mu0)
    a_s = albedo_suelo
    sistema = np.array([[vp[1] * e_k, vm[1]],
                        [vp[0] - a_s * vp[1], (vm[0] - a_s * vm[1]) * e_k]])
    lado = np.array([-c[1], a_s * mu0 * S * e_0 - (c[0] - a_s * c[1]) * e_0])
    A_, B_ = np.linalg.solve(sistema, lado)

    def flujos(t):
        y = A_ * vp * math.exp(k * (t - tau)) + B_ * vm * math.exp(-k * t) + c * math.exp(-t / mu0)
        return y[0], y[1]
    arriba_0, _ = flujos(0.0)
    arriba_b, abajo_b = flujos(tau)
    llega = mu0 * S
    R = arriba_0 / llega
    # lo que llega al suelo: el haz escalado (que incluye la luz dispersada hacia delante, en el pico) mas el difuso.
    # Se reparte en el directo DE VERDAD (sin escalar) y el resto, difuso.
    T_total = e_0 + abajo_b / llega
    T_dir = math.exp(-tau_real / mu0)
    T_dif = T_total - T_dir
    A_suelo = (1.0 - a_s) * (T_dir + T_dif)
    A_capa = 1.0 - R - (T_dir + T_dif) + arriba_b / llega
    return {"R": R, "T_dir": T_dir, "T_dif": T_dif, "A_capa": A_capa, "A_suelo": A_suelo}


# ---------------------------------------------------------------------------------------------- varias capas
def columna(taus, omegas, gs, mu0, albedo_suelo=0.0, delta=True):
    """La columna entera (capas de arriba abajo), resuelta EXACTA para las ecuaciones de dos flujos: en cada capa la
    solucion general (dos coeficientes) mas la particular del haz que le llega, con la continuidad de los dos flujos
    difusos en cada separacion entre capas y las dos condiciones de contorno (nada difuso entra por arriba; suelo de
    Lambert). Es un sistema lineal de 2N ecuaciones (equivale al metodo de "sumar capas", pero sin aproximaciones de
    orden). Flujos normalizados por la luz que llega (mu0 por unidad de area horizontal). Devuelve:
      R, T_dir (directa de verdad en el suelo), T_dif, A_capas (absorbida en cada capa), A_suelo,
      arriba / abajo_difuso / abajo_directo: flujos en las N+1 separaciones (de arriba abajo)."""
    taus = np.atleast_1d(np.asarray(taus, dtype=float))
    omegas = np.broadcast_to(np.asarray(omegas, dtype=float), taus.shape)
    gs = np.broadcast_to(np.asarray(gs, dtype=float), taus.shape)
    n = len(taus)
    tau_real_total = float(taus.sum())
    if delta:
        esc = [escalado_delta(t, w, g) for t, w, g in zip(taus, omegas, gs)]
        taus = np.array([e[0] for e in esc]); omegas = np.array([e[1] for e in esc]); gs = np.array([e[2] for e in esc])
    omegas = np.minimum(omegas, OMEGA_MAX)
    S = 1.0
    techo = np.concatenate([[0.0], np.cumsum(taus)])          # profundidad optica (escalada) de cada separacion
    capas = []
    for i in range(n):
        g1, g2, g3, g4 = gammas_eddington(omegas[i], gs[i], mu0)
        k = math.sqrt(max(g1 * g1 - g2 * g2, 0.0))
        m = mu0 * (1.0 + 2e-6) if abs(k * mu0 - 1.0) < 1e-6 else mu0
        M = np.array([[g1, -g2], [g2, -g1]])
        b = np.array([-omegas[i] * S * g3, omegas[i] * S * g4]) * math.exp(-techo[i] / m)
        c = np.linalg.solve(-np.eye(2) / m - M, b)
        capas.append({"k": k, "vp": np.array([g2, g1 - k]), "vm": np.array([g2, g1 + k]), "c": c, "mu": m,
                      "tau": taus[i]})

    def base(cp, t):
        """Columnas de A y B y termino fijo de y(t) en la capa (t local, de 0 a tau)."""
        return (cp["vp"] * math.exp(cp["k"] * (t - cp["tau"])), cp["vm"] * math.exp(-cp["k"] * t),
                cp["c"] * math.exp(-t / cp["mu"]))
    Mx = np.zeros((2 * n, 2 * n))
    rhs = np.zeros(2 * n)
    a, bb, f = base(capas[0], 0.0)                            # arriba: F-(0) = 0
    Mx[0, 0], Mx[0, 1], rhs[0] = a[1], bb[1], -f[1]
    fila = 1
    for i in range(n - 1):                                    # continuidad de F+ y F- entre la capa i y la i+1
        a1, b1, f1 = base(capas[i], capas[i]["tau"])
        a2, b2, f2 = base(capas[i + 1], 0.0)
        for comp in (0, 1):
            Mx[fila, 2 * i], Mx[fila, 2 * i + 1] = a1[comp], b1[comp]
            Mx[fila, 2 * i + 2], Mx[fila, 2 * i + 3] = -a2[comp], -b2[comp]
            rhs[fila] = f2[comp] - f1[comp]
            fila += 1
    a, bb, f = base(capas[-1], capas[-1]["tau"])               # suelo de Lambert
    A_s = albedo_suelo
    directo_suelo = mu0 * S * math.exp(-techo[-1] / mu0)
    Mx[fila, 2 * n - 2] = a[0] - A_s * a[1]
    Mx[fila, 2 * n - 1] = bb[0] - A_s * bb[1]
    rhs[fila] = A_s * directo_suelo - (f[0] - A_s * f[1])
    coef = np.linalg.solve(Mx, rhs)
    arriba, abajo = np.zeros(n + 1), np.zeros(n + 1)
    for i, cp in enumerate(capas):
        for j, t in ((i, 0.0), (i + 1, cp["tau"])):
            a, bb, f = base(cp, t)
            y = coef[2 * i] * a + coef[2 * i + 1] * bb + f
            arriba[j], abajo[j] = y[0], y[1]
    llega = mu0 * S
    directo = mu0 * S * np.exp(-techo / mu0)                  # haz escalado (incluye el pico hacia delante)
    neto = (abajo + directo - arriba) / llega                 # flujo neto hacia abajo en cada separacion
    A_capas = neto[:-1] - neto[1:]
    T_total = (abajo[-1] + directo[-1]) / llega
    T_dir = math.exp(-tau_real_total / mu0)
    return {"R": arriba[0] / llega, "T_dir": T_dir, "T_dif": T_total - T_dir, "A_capas": A_capas,
            "A_suelo": (1.0 - A_s) * T_total, "arriba": arriba / llega, "abajo_difuso": abajo / llega,
            "abajo_directo": directo / llega}

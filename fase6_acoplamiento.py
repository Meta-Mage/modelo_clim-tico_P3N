# fase6_acoplamiento.py -- Fase 6.3 (paso 5, primeras piezas): utilidades para acoplar la fisica de columna
# de M3N (que trabaja en el CENTRO de cada celda) con el nucleo dinamico (rejilla C: u en las caras oeste,
# v en las caras sur). Diseno: DISENO_FASE_6.3.md §6.7.
#
# 1. viento_en_centros: u, v de las caras -> centros (media de las dos caras; en las filas polares la cara
#    del polo tiene v = 0, como en el nucleo).
# 2. tendencia_a_caras: el cambio del viento calculado en los centros (capa limite, rozamiento) -> caras
#    (media de los dos centros vecinos). Es la traspuesta de la anterior en el producto escalar SIN pesos
#    (prueba exacta en test_fase6_acoplamiento.py). Con los pesos de masa y area del nucleo no conserva
#    exactamente ni el momento ni la energia: la energia la cierra calor_rozamiento_exacto y el momento se
#    mide (DISENO_FASE_6.3.md §6.7).
# 3. calor_rozamiento_exacto: DIAGNOSTICO. Calentamiento que devolveria exactamente la energia cinetica perdida
#    con la definicion del nucleo (NucleoSeco.energia_cinetica), dT = -(K_1 - K_0)/cp. Medido el 07/10 tras 30
#    dias de Held y Suarez: aplicarlo celda a celda cierra la energia, pero enfria un 1,2 % de las (celda, capa)
#    (hasta 5e-3 K por paso), porque llevar el viento a las caras mezcla el momento de columnas vecinas; y
#    escalarlo por columnas esta muy mal condicionado (factores de 0,46 a > 10^5). DECISION (🔶, §6.7): el
#    calor por rozamiento es el de la capa limite, LOCAL y POSITIVO en los centros; esta funcion mide el
#    residuo global (~0,2 % del rozamiento, ~0,004 W/m2), que se declara en el balance.

import numpy as np


def viento_en_centros(u, v):
    """u (..., F, C) en las caras oeste; v (..., F-1, C) en las caras sur de las filas 0..F-2.
    Devuelve (u_c, v_c), (..., F, C)."""
    u_c = 0.5 * (u + np.roll(u, -1, axis=-1))
    forma = u.shape
    v_c = np.zeros(forma)
    v_c[..., :-1, :] += 0.5 * v          # cara sur de las filas 0..F-2
    v_c[..., 1:, :] += 0.5 * v           # cara norte de las filas 1..F-1
    return u_c, v_c


def tendencia_a_caras(du_c, dv_c):
    """Traspuesta de viento_en_centros: du (caras oeste) = media de los centros oeste y este de la cara;
    dv (caras sur de las filas 0..F-2) = media de los centros norte y sur de la cara."""
    du = 0.5 * (du_c + np.roll(du_c, 1, axis=-1))
    dv = 0.5 * (dv_c[..., :-1, :] + dv_c[..., 1:, :])
    return du, dv


def calor_rozamiento_exacto(nucleo, u0, v0, u1, v1, cp):
    """Calentamiento (K) que devuelve como calor la energia cinetica perdida, con la definicion del nucleo."""
    return -(nucleo.energia_cinetica(u1, v1) - nucleo.energia_cinetica(u0, v0)) / cp

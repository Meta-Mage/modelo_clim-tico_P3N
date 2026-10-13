# test_fase52_k.py -- pruebas de las tablas k (Fase 5.2, borrador). Casos con solucion conocida.
import math

import numpy as np
from scipy.special import i0e, i1e

import fase52_k as K
import fase52_lineas as L


def _linea_lorentz(S=1.0, alfa=0.07, ancho=None, n=400001):
    """Una linea de Lorentz (Doppler despreciable) en una banda de +-ancho/2 cm^-1: (nu, dnu, k)."""
    ancho = ancho or 2000 * alfa
    nu = np.linspace(-ancho / 2, ancho / 2, n)
    return nu, nu[1] - nu[0], S * L.voigt(nu, 0.0, 1e-8, alfa)


def test_nodos_g_suman_uno_e_integran_polinomios():
    g, w = K.nodos_g()
    assert len(g) == 16 and abs(w.sum() - 1) < 1e-14 and np.all((g > 0) & (g < 1))
    for grado in range(4):                                      # Gauss en cada tramo: exacto para polinomios bajos
        assert abs((w * g ** grado).sum() - 1 / (grado + 1)) < 1e-14


def test_una_capa_es_exacta_salvo_la_cuadratura():
    """En una capa homogenea, ordenar k no cambia la transmision: con muchos puntos g se recupera la de linea a
    linea; con los 16 de T6, el error es pequeño en todo el rango de caminos (de linea debil a saturada)."""
    nu, dnu, k = _linea_lorentz()
    muchos = tuple((i / 200, (i + 1) / 200, 4) for i in range(200))
    for u in (1e-3, 1e-1, 1.0, 10.0, 1e3):
        exacta = K.transmision_lbl(k, u)
        g, w, kg = K.tabla_k(k, tramos=muchos)
        assert abs(K.transmision_ck(w, kg, u) - exacta) < 2e-4
        g, w, kg = K.tabla_k(k)
        assert abs(K.transmision_ck(w, kg, u) - exacta) < 5e-4


def test_linea_aislada_ladenburg_reiche():
    """Anchura equivalente de una linea de Lorentz aislada (Ladenburg y Reiche 1913; Goody y Yung 1989, §4.4 ⚠️):
    W = 2 pi alfa L(x), x = S u / (2 pi alfa), L(x) = x e^-x [I0(x) + I1(x)]. Se comprueba el calculo linea a linea
    de la banda (la base de las tablas k) en los regimenes debil (W = S u), intermedio y fuerte (W = 2 sqrt(S u alfa))."""
    S, alfa = 1.0, 0.07
    nu, dnu, k = _linea_lorentz(S, alfa, ancho=4000 * alfa, n=800001)
    ancho = nu[-1] - nu[0] + dnu
    for u in (0.01, 1.0, 30.0):
        x = S * u / (2 * math.pi * alfa)
        W = 2 * math.pi * alfa * x * (i0e(x) + i1e(x))
        W_lbl = (1 - K.transmision_lbl(k, u)) * ancho
        # las alas fuera de la banda faltan en el calculo: ~ 2 S u alfa / (pi (ancho/2)) en el limite debil de las alas
        falta = 2 * S * u * alfa / (math.pi * ancho / 2)
        assert abs(W_lbl + falta - W) / W < 1e-6, (u, W_lbl, W)


def test_correlacionada_exacta_con_una_linea_a_dos_presiones():
    """Una linea en dos capas a distinta presion: k de una capa es funcion monotona de |nu - nu0|, igual que la de la
    otra, asi que el orden es el mismo y la hipotesis correlacionada es exacta (salvo la cuadratura)."""
    nu = np.linspace(-100, 100, 400001)
    k1 = L.voigt(nu, 0.0, 1e-3, 0.07)          # capa baja: ancha por presion
    k2 = L.voigt(nu, 0.0, 1e-3, 0.007)         # capa alta: diez veces menos presion
    muchos = tuple((i / 200, (i + 1) / 200, 4) for i in range(200))
    g, w, kg = K.tabla_k_capas(np.array([k1, k2]), tramos=muchos)
    for u in ((0.5, 0.05), (5.0, 0.5), (50.0, 5.0)):
        assert abs(K.transmision_ck(w, kg, u) - K.transmision_lbl(np.array([k1, k2]), u)) < 3e-4


def test_correlacionada_falla_con_lineas_que_no_se_corresponden():
    """Y su limite conocido: si la linea de una capa cae donde la otra no absorbe, el orden no se corresponde y la
    hipotesis sobreestima el solape (la transmision correlacionada sale MAYOR que la verdadera). Documenta por que
    el error se medira contra el calculo linea a linea (DISENO_FASE_5.2.md §4)."""
    nu = np.linspace(-50, 50, 200001)
    k1 = L.voigt(nu, -10.0, 1e-3, 0.1)
    k2 = L.voigt(nu, +10.0, 1e-3, 0.1)
    u = (20.0, 20.0)
    exacta = K.transmision_lbl(np.array([k1, k2]), u)
    muchos = tuple((i / 200, (i + 1) / 200, 4) for i in range(200))
    g, w, kg = K.tabla_k_capas(np.array([k1, k2]), tramos=muchos)
    assert K.transmision_ck(w, kg, u) > exacta + 0.01

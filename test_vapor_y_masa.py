# test_vapor_y_masa.py -- 3.18.0 (DISENO_FASE_6.3.md §6.20): el vapor como razon de mezcla en la fisica con I16
# (§6.20.1) y la columna observada / masa de aire del modo Tierra (§6.20.2).
#
# Uso:  python -m pytest test_vapor_y_masa.py      (unos segundos)

import numpy as np
import pytest

import fase30_multicapa as M30
import fase31_agua as A
import fase2b_atmosfera as F
from test_v31 import columna_humeda

G = 9.80665


def test_razon_de_mezcla_de_saturacion_y_su_derivada():
    T = np.array([255.0, 285.0, 300.0]); p = np.array([1.0e5, 1.0e5, 1.005e5])
    ms, dms = A.ms_y_derivada(T, p)
    e = A.es_agua(T)
    assert np.allclose(ms, A.EPSILON * e / (p - e), rtol=1e-6)          # eps*e/(p - e), con la e_s de Ambaum
    ms2, _ = A.ms_y_derivada(T + 1e-4, p)
    assert np.allclose(dms, (ms2 - ms) / 1e-4, rtol=1e-3)
    qs, _ = A.qs_y_derivada(T, p)
    assert (ms > qs).all()                                               # por kg de aire seco: siempre mayor
    assert np.array_equal(A.saturacion(T, p)[0], qs)                     # por defecto, la de la v3.1


def test_presion_total_suma_el_peso_del_vapor():
    C = M30.Columna(np.zeros((1, 1)), gravedad=G)
    q = np.full((C.n, 1, 1), 0.01)
    ph_t, pm_t = A.presion_total(C.ph, C.pm, C.dp, q)
    assert ph_t[0, 0, 0] == C.ph[0, 0, 0]                               # en el tope, nada encima
    assert ph_t[-1, 0, 0] == pytest.approx(C.ps[0, 0] * 1.01)           # en superficie, p_seca*(1 + m)
    assert (pm_t > C.pm).all() and (pm_t < ph_t[1:]).all() and (pm_t > ph_t[:-1]).all()


def test_conveccion_y_condensacion_en_razon_de_mezcla_conservan():
    C, T, q = columna_humeda()
    ph_t, pm_t = A.presion_total(C.ph, C.pm, C.dp, q)
    T2, q2, P = A.conveccion_humeda(T, q, pm_t, ph_t, 992.0, G, razon=True, ph_masa=C.ph)
    assert (P > 0).all()
    energia = (((T2 - T) * A.CP_AIRE + (q2 - q) * A.L_V) * C.dp).sum(axis=0) / G   # con la masa SECA del nucleo
    agua = ((q2 - q) * C.dp).sum(axis=0) / G + P
    assert np.allclose(energia, 0.0, atol=1e-6)
    assert np.allclose(agua, 0.0, atol=1e-12)
    C, T, q = columna_humeda(rh=1.3)
    pm_t = A.presion_total(C.ph, C.pm, C.dp, q)[1]
    T2, q2, P = A.condensacion_gran_escala(T, q, pm_t, C.dp, G, razon=True)
    assert np.allclose((((T2 - T) * A.CP_AIRE + (q2 - q) * A.L_V) * C.dp).sum(axis=0), 0.0, atol=1e-6)
    assert np.allclose(((q - q2) * C.dp).sum(axis=0) / G, P)
    ms2, _ = A.ms_y_derivada(T2, pm_t)
    assert (q2 <= ms2 * 1.001).all()


def test_conveccion_sin_razon_es_la_de_la_v31():
    C, T, q = columna_humeda()
    a = A.conveccion_humeda(T, q, C.pm, C.ph, 992.0, G)
    b = A.conveccion_humeda(T, q, C.pm, C.ph, 992.0, G, razon=False, ph_masa=C.ph)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))


def test_masa_de_aire_del_modo_tierra_y_p3n_sin_cambios():
    alt = np.zeros((36, 72)); alt[5:9, 10:20] = 3000.0; alt[30:, :] = 2500.0
    peso = np.cos(np.radians(np.linspace(87.5, -87.5, 36))).reshape(-1, 1) * np.ones((1, 72))
    pnm = M30.p_nivel_mar_para_masa(alt, peso, 98305.0, G)
    ps = M30.presion_superficie(alt, G, pnm)
    assert float((ps * peso).sum() / peso.sum()) == pytest.approx(98305.0, abs=1e-6)
    ref = M30.P_REF * np.exp(-G * alt / (M30.R_AIRE * M30.T_ALTURA_ESCALA))
    assert np.array_equal(M30.presion_superficie(alt, G), ref)           # P3N: P_NIVEL_MAR = 1 bar, como siempre
    assert M30.P_NIVEL_MAR == M30.P_REF == M30.P0 == 1.0e5


def test_columna_observada_normaliza_con_su_presion():
    a = M30.Columna(np.zeros((1, 1)), gravedad=G)
    b = M30.Columna(np.zeros((1, 1)), gravedad=G, p_nivel_mar=M30.P_COLUMNA_OBSERVADA, columna_observada=True)
    assert b.p_tau == M30.P_COLUMNA_OBSERVADA and a.p_tau == M30.P_REF
    assert (a.tau0, a.f_lw) == (M30.TAU0_LW, M30.F_LW) and (b.tau0, b.f_lw) == (M30.TAU0_LW_OBS, M30.F_LW_OBS)
    x = b.ph / M30.P_COLUMNA_OBSERVADA
    assert np.allclose(b.tau, M30.TAU0_LW_OBS * (M30.F_LW_OBS * x + (1 - M30.F_LW_OBS) * x ** M30.EXP_VAPOR_LW), rtol=1e-12)
    assert b.ps[0, 0] / b.p_tau == pytest.approx(1.0)


def test_tau_fijo_recalibrado_da_los_flujos_de_wild():
    """3.18.1: con la columna observada y N = 160 (el metodo de la v3.0), DLR 314 y OLR 267 W/m2 (Wild et al. 2019)."""
    import calibrar_columna_observada as K
    C = K.columna(True, 160)
    C.tau0, C.f_lw = M30.TAU0_LW_OBS, M30.F_LW_OBS
    C.actualizar_ps(C.ps)
    olr, dlr = K.flujos(C, K.temperatura(C), 289.0)
    assert olr == pytest.approx(267.0, abs=0.05) and dlr == pytest.approx(314.0, abs=0.05)


def test_interruptores_de_la_3_18_y_su_alcance():
    assert F.VAPOR_RAZON_MEZCLA_I16 is True and F.COLUMNA_OBSERVADA_I16 is True
    # sin I16 no se usan (el clima oficial, de 2 capas, no los ve): lo comprueba la prueba bit a bit de la 3.18.0
    assert F.PRESION_SECA_MEDIA_TIERRA is None or F.PRESION_SECA_MEDIA_TIERRA == 98305.0

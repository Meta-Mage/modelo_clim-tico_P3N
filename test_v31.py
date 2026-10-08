# test_v31.py -- pruebas rapidas de la v3.1 (Fase 5.1: ciclo del agua).
#
# Uso:  python -m pytest test_v31.py      (menos de un minuto)
#
# Lo que tiene que ser EXACTO: la conveccion humeda y la condensacion conservan
# la energia (cp*T + L*q) y el agua; una simulacion corta con todo el ciclo del
# agua cierra el balance de energia y los dos balances de agua; los
# interruptores nuevos estan apagados por defecto; las funciones auxiliares
# (saturacion, espectro, suelo) dan lo que deben.

import math

import numpy as np
import pytest

import fase30_multicapa as M30
import fase31_agua as A
import fase2b_atmosfera as F

G = 9.80665


def columna_humeda(Ts=301.0, rh=0.9):
    C = M30.Columna(np.zeros((1, 2)), gravedad=G)
    pm = C.pm[:, 0, 0]
    T = np.empty(C.n); T[-1] = Ts
    for k in range(C.n - 2, -1, -1):                   # adiabatico humedo aproximado, 1 K mas frio
        Tk = T[k + 1]
        for _ in range(2):
            qs, _ = A.qs_y_derivada(np.array([Tk]), np.array([0.5 * (pm[k] + pm[k + 1])]))
            r = qs[0] / (1 - qs[0])
            a = A.KAPPA * Tk + A.L_V / A.CP_AIRE * r
            b = A.L_V ** 2 * r / (A.CP_AIRE * A.R_V * Tk ** 2)
            Tk = T[k + 1] + a / (1 + b) * math.log(pm[k] / pm[k + 1])
        T[k] = max(Tk, 200.0)
    T = np.repeat((T - 1.0)[:, None, None], 2, axis=2)
    T[-1] += 1.5
    qs, _ = A.qs_y_derivada(T, C.pm)
    return C, T, rh * qs


def test_conveccion_humeda_conserva_energia_y_agua():
    C, T, q = columna_humeda()
    T2, q2, P = A.conveccion_humeda(T, q, C.pm, C.ph, 992.0, G)
    assert (P > 0).all()                                # conveccion profunda: llueve
    energia = (((T2 - T) * A.CP_AIRE + (q2 - q) * A.L_V) * C.dp).sum(axis=0) / G
    agua = ((q2 - q) * C.dp).sum(axis=0) / G + P
    assert np.allclose(energia, 0.0, atol=1e-6)         # J/m2 (frente a ~1e6 J/m2 movidos)
    assert np.allclose(agua, 0.0, atol=1e-12)


def test_conveccion_no_actua_en_columna_seca_y_estable():
    C = M30.Columna(np.zeros((1, 1)), gravedad=G)
    T = C.perfil_inicial(np.full((1, 1), 280.0)) + np.linspace(30, 0, C.n)[:, None, None]
    q = np.full_like(T, 1e-6)
    T2, q2, P = A.conveccion_humeda(T, q, C.pm, C.ph, 992.0, G)
    assert np.array_equal(T2, T) and np.array_equal(q2, q) and P[0, 0] == 0.0


def test_condensacion_conserva_y_deja_saturado():
    C, T, q = columna_humeda(rh=1.3)
    T2, q2, P = A.condensacion_gran_escala(T, q, C.pm, C.dp, G)
    energia = (((T2 - T) * A.CP_AIRE + (q2 - q) * A.L_V) * C.dp).sum(axis=0)
    assert np.allclose(energia, 0.0, atol=1e-6)
    assert np.allclose(((q - q2) * C.dp).sum(axis=0) / G, P)
    qs2, _ = A.qs_y_derivada(T2, C.pm)
    assert (q2 <= qs2 * 1.001).all()                   # un paso de Newton: casi exactamente saturado


def test_saturacion_valores_de_referencia():
    # Ambaum (2020) frente a valores conocidos: 611,2 Pa en el punto triple; ~4246 Pa a 30 C (agua)
    assert A.es_agua(np.array([273.16]))[0] == pytest.approx(611.2)
    assert A.es_agua(np.array([303.15]))[0] == pytest.approx(4246, rel=0.005)
    assert A.es_hielo(np.array([253.15]))[0] < A.es_agua(np.array([253.15]))[0]


def test_derivada_de_saturacion():
    T = np.array([250.0, 280.0, 305.0]); p = np.full(3, 9e4)
    qs, dqs = A.qs_y_derivada(T, p)
    qs2, _ = A.qs_y_derivada(T + 1e-4, p)
    assert np.allclose(dqs, (qs2 - qs) / 1e-4, rtol=1e-3)


def test_fraccion_visible_de_planck():
    # integral numerica de la ley de Planck: el Sol emite ~48,8 % por debajo de 0,7 um
    assert A.fraccion_visible(5772.0) == pytest.approx(0.48815, abs=2e-4)
    assert A.fraccion_visible(5420.0) < A.fraccion_visible(5772.0)
    alb = A.albedos_estrella(5420.0)
    assert alb["nieve_fria"] < A.albedos_estrella(5772.0)["nieve_fria"] and alb["factor_hielo"] < 1


def test_suelo_seco_y_saturado():
    lam, C = A.propiedades_suelo(np.array([0.0, A.W_CAMPO]))
    inercia = np.sqrt(lam * C)
    assert inercia[0] == pytest.approx(531, rel=0.02) and inercia[1] == pytest.approx(2577, rel=0.02)


def test_fraccion_nieve_rampa():
    assert np.allclose(A.fraccion_nieve(np.array([270.0, 273.15, 274.15, 275.15, 280.0])), [1, 1, 0.5, 0, 0])


def test_interruptores_v31_apagados():
    for k in ("ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral"):
        assert F.INTERRUPTORES_FASE2B[k] is False


def _simulacion_corta(con_agua, vapor_radiativo=False):
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA
    from cache_simulacion import precalcular_orbita_cacheada
    import parametros as P
    from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS
    orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    pasos_dia = round(P.ROTACION_PERIODO / 992)
    orb = orb[:2 * pasos_dia]
    lon = np.ones((FILAS, 1)) * LONGITUDES_GRADOS[None]
    lat = LATITUDES_GRADOS[:, None] * np.ones((1, COLUMNAS))
    tipo = np.where((np.abs(lon) < 40) & (np.abs(lat) < 70), TIERRA, AGUA)
    alt = np.where(tipo == TIERRA, 600.0, 0.0)
    I = dict(F.INTERRUPTORES_FASE2B)
    I["atmosfera_multicapa"] = True
    for k in ("ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral"):
        I[k] = con_agua
    I["vapor_radiativo"] = vapor_radiativo
    return F.simular_fase2b(orb, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA,
                         0.55, interruptores=I, max_anos=1, acelerar=False)


def test_simulacion_corta_cierra_energia_y_agua():
    r = _simulacion_corta(True)
    assert abs(r["energia"]["cierre_relativo"]) < 1e-9
    a = r["agua"]
    assert abs(a["cierre_agua_atmosfera_kg_m2"]) < 1e-9 and abs(a["cierre_agua_tierra_kg_m2"]) < 1e-9
    assert a["recortes_kg_m2"] == 0.0
    assert (a["precipitacion"] >= 0).all() and (a["cubo_medio"] <= A.W_CAMPO + 1e-9).all()


def test_v30_sin_agua_tambien_cierra_la_energia():
    r = _simulacion_corta(False)
    assert r["agua"] is None and abs(r["energia"]["cierre_relativo"]) < 1e-9


def test_vapor_radiativo_cierra_energia_y_agua():
    r = _simulacion_corta(True, vapor_radiativo=True)
    assert abs(r["energia"]["cierre_relativo"]) < 1e-9
    assert abs(r["agua"]["cierre_agua_atmosfera_kg_m2"]) < 1e-9


def test_tau_vapor_de_referencia():
    # con el vapor de referencia (24,9 kg/m2), el espesor total es el de la v3.0 (~2,45)
    C = M30.Columna(np.zeros((1, 1)), gravedad=G)
    q = np.full((C.n, 1, 1), 0.0)
    C.actualizar_tau_vapor(q, A.A_LW_SECO, A.B_LW_VAPOR)
    seco = float(-np.log(C.trans).sum())
    assert seco == pytest.approx(A.A_LW_SECO, rel=1e-6)


def test_punto_de_control_da_lo_mismo(tmp_path):
    """Una simulacion interrumpida y retomada desde el punto de control da EXACTAMENTE lo mismo."""
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA
    from cache_simulacion import precalcular_orbita_cacheada
    import parametros as P
    from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS
    orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    orb = orb[:round(P.ROTACION_PERIODO / 992)]                     # "años" de un dia
    lon = np.ones((FILAS, 1)) * LONGITUDES_GRADOS[None]
    lat = LATITUDES_GRADOS[:, None] * np.ones((1, COLUMNAS))
    tipo = np.where((np.abs(lon) < 40) & (np.abs(lat) < 70), TIERRA, AGUA)
    alt = np.where(tipo == TIERRA, 600.0, 0.0)
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral"):
        I[k] = True
    arg = (orb, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55)
    kw = dict(interruptores=I, acelerar=False, tolerancia_convergencia=1e-9)
    seguida = F.simular_fase2b(*arg, max_anos=3, **kw)
    archivo = str(tmp_path / "estado.pkl")
    F.simular_fase2b(*arg, max_anos=2, archivo_estado=archivo, **kw)      # se "interrumpe" tras 2 años
    retomada = F.simular_fase2b(*arg, max_anos=3, archivo_estado=archivo, **kw)
    for k in ("reg_media", "suelo_media", "hielo_espesor"):
        assert np.array_equal(seguida[k], retomada[k]), k
    assert np.array_equal(seguida["agua"]["precipitacion"], retomada["agua"]["precipitacion"])

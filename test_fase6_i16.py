# test_fase6_i16.py -- v3.1-pre11: pruebas del modelo ACOPLADO (interruptor I16, "nucleo_dinamico").
#
# Uso:  python -m pytest test_fase6_i16.py      (alrededor de un minuto)
#
# Son simulaciones MUY cortas (unos pocos pasos de la fisica por "año"): comprueban propiedades que deben
# cumplirse en CADA paso, no el clima.
#   - el agua de la atmosfera y la del suelo se conservan EXACTAMENTE (redondeo);
#   - el balance de energia cierra (los residuos declarados son mucho menores que la cota);
#   - una simulacion interrumpida y retomada desde el punto de control da EXACTAMENTE lo mismo;
#   - con I16 apagado no se crea nada del nucleo (lo de siempre).
# El clima no se puede juzgar aqui: con una orbita recortada a unos pocos pasos, el estado inicial (el
# balance de Newton con la luz MEDIA de esa orbita) tiene la superficie a mas de 400 K bajo la estrella.
# La estabilidad durante semanas se comprueba aparte, desde un estado de equilibrio (DISENO_FASE6_3.md §6.9).

import numpy as np
import pytest

import fase2b_atmosfera as F

PASOS_ANO = 4          # pasos de la fisica por "año" de prueba (muy por debajo de los ~25 del arranque frio)


def _argumentos(acua=False):
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA
    from cache_simulacion import precalcular_orbita_cacheada
    import parametros as P
    from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS
    orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)[:PASOS_ANO]
    lon = np.ones((FILAS, 1)) * LONGITUDES_GRADOS[None]
    lat = LATITUDES_GRADOS[:, None] * np.ones((1, COLUMNAS))
    tipo = np.where((np.abs(lon) < 40) & (np.abs(lat) < 70), TIERRA, AGUA)
    if acua:
        tipo = np.full_like(tipo, AGUA)
    alt = np.where(tipo == TIERRA, 600.0, 0.0)
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral",
              "nucleo_dinamico"):
        I[k] = True
    arg = (orb, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55)
    kw = dict(interruptores=I, acelerar=False, tolerancia_convergencia=1e-9)
    return arg, kw


@pytest.fixture(scope="module")
def corrida():
    arg, kw = _argumentos()
    return F.simular_fase2b(*arg, max_anos=2, **kw)


def test_i16_agua_exacta(corrida):
    a = corrida["agua"]
    assert abs(a["cierre_agua_atmosfera_kg_m2"]) < 1e-12        # kg/m2 (la columna tiene ~10-60 kg/m2)
    assert abs(a["cierre_agua_tierra_kg_m2"]) < 1e-12
    assert a["recortes_kg_m2"] == 0.0                             # el transporte no deja vapor negativo


def test_i16_energia_cierra_salvo_la_dinamica(corrida):
    """Quitando lo que la dinamica no conserva (y el calor de rozamiento que devuelve la capa limite), el
    balance de energia cierra como en la v3.1 sin dinamica (redondeo). El error de la dinamica se mide y se
    declara (DISENO_FASE6_3.md §6.9); en estas simulaciones de pocos pasos desde el reposo es grande porque
    el aire se esta poniendo en movimiento a ~1500 W/m2."""
    e = corrida["energia"]
    assert abs(e["cierre_sin_dinamica"]) < 1e-9
    assert np.isfinite(e["error_dinamica_W_m2"]) and e["calor_rozamiento_W_m2"] > 0


def test_i16_punto_de_control_da_lo_mismo(tmp_path):
    arg, kw = _argumentos()
    seguida = F.simular_fase2b(*arg, max_anos=3, **kw)
    archivo = str(tmp_path / "estado.pkl")
    F.simular_fase2b(*arg, max_anos=2, archivo_estado=archivo, **kw)
    retomada = F.simular_fase2b(*arg, max_anos=3, archivo_estado=archivo, **kw)
    for k in ("reg_media", "suelo_media", "hielo_espesor"):
        assert np.array_equal(seguida[k], retomada[k]), k
    assert np.array_equal(seguida["agua"]["precipitacion"], retomada["agua"]["precipitacion"])
    assert np.array_equal(seguida["flujos"]["T_atm"], retomada["flujos"]["T_atm"])
    assert seguida["energia"]["cierre_relativo"] == retomada["energia"]["cierre_relativo"]


def test_i16_nucleo_en_reposo_isotermo_sobre_montanas():
    """T4 (DISENO_FASE6_3.md §6.7.7): el nucleo con la configuracion de I16 (niveles de M3N, avanzar con
    hiperdifusion, su energia cinetica como calor y la correccion a superficies de presion) mantiene en
    reposo una atmosfera isoterma sobre las montañas del caso 5 de Williamson (redondeo)."""
    import fase30_multicapa as M30
    from fase6_nucleo import NucleoSeco
    from fase6_aguas_someras import caso5, A_TIERRA, OMEGA_TIERRA, G_TIERRA
    hs = caso5()[3]
    m = NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, M30.sigma_seminiveles(20), phis=G_TIERRA * hs)
    m.preparar_semiimplicito(450.0)
    m.preparar_hiperdifusion(0.5, calor_rozamiento=True, correccion_presion=True)
    T0 = 280.0
    ps = 1e5 * np.exp(-G_TIERRA * hs / (287.0 * T0))
    ant, act = m.arrancar((ps, np.full((m.N, 36, 72), T0), np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72))))
    for _ in range(100):
        ant, act = m.avanzar(ant, act)
    assert np.abs(act[2]).max() < 1e-9 and np.abs(act[3]).max() < 1e-9
    assert np.abs(act[1] - T0).max() < 1e-9 and np.abs(act[0] - ps).max() < 1e-6

# test_v30.py -- pruebas rapidas de la v3.0 (atmosfera de N capas, modo Tierra).
#
# Uso:  python -m pytest test_v30.py      (unos segundos)
#
# Comprueban, sin simular un clima, lo que tiene que ser EXACTO: que el
# infrarrojo y el ajuste convectivo conservan la energia, que el ajuste
# deja la columna estable, que el diagnostico de transporte cierra, que
# el interruptor I10 sigue apagado por defecto y que el modo Tierra pone
# los valores terrestres.

import os
import subprocess
import sys

import numpy as np
import pytest

import fase30_multicapa as M30
import fase2b_atmosfera as F
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS


def columna(altitudes=(0.0, 2300.0, 5000.0)):
    return M30.Columna(np.array(altitudes, dtype=float)[None, :])


def test_presion_superficie_al_nivel_del_mar():
    assert M30.presion_superficie(np.array([0.0]))[0] == pytest.approx(M30.P0)


def test_capas_suman_toda_la_masa():
    c = columna()
    assert np.allclose((c.dp).sum(axis=0), c.ps)
    assert np.allclose(c.peso_sw.sum(axis=0), 1.0)


def test_infrarrojo_conserva_la_energia():
    # lo que ganan las capas + lo que gana el suelo + lo que sale al espacio = 0 (no se crea ni se pierde nada)
    rng = np.random.default_rng(1)
    c = columna()
    T = 200 + 90 * rng.random((c.n,) + c.forma)
    Ts = 290 + 10 * rng.random(c.forma)
    D, Bh = c.infrarrojo_bajada(T, c.aire_superficie(T))
    U = c.infrarrojo_subida(Bh, M30.CONSTANTE_SB * Ts ** 4)
    calor = c.calentamiento_infrarrojo(U, D).sum(axis=0)
    suelo = D[-1] - U[-1]
    assert np.allclose(calor + suelo + U[0], 0.0, atol=1e-9)


def test_ajuste_convectivo_conserva_energia_y_deja_columna_estable():
    rng = np.random.default_rng(2)
    c = columna()
    T = 250 + 60 * rng.random((c.n,) + c.forma)          # columnas muy inestables
    Ta = c.ajuste_convectivo(T)
    assert np.allclose((Ta * c.cap).sum(axis=0), (T * c.cap).sum(axis=0), rtol=1e-13)
    th = Ta / c.pi_ajuste
    assert (th[:-1] >= th[1:] - 1e-9).all()               # theta* no decrece hacia arriba
    assert np.array_equal(c.ajuste_convectivo(Ta), Ta) or np.allclose(c.ajuste_convectivo(Ta), Ta, atol=1e-10)


def test_ajuste_no_toca_una_columna_estable():
    c = columna()
    T = c.perfil_inicial(np.full(c.forma, 288.0)) + 0.5 * np.arange(c.n)[::-1, None, None] * 0
    T = T + np.linspace(5, 0, c.n)[:, None, None]          # mas calida arriba: estable
    assert np.array_equal(c.ajuste_convectivo(T), T)


def test_transporte_meridional_cierra():
    # con un operador de difusion conservativo, lo que recibe todo el planeta es 0
    from fase2_difusion import construir_matriz_difusion
    rng = np.random.default_rng(3)
    L = construir_matriz_difusion(np.full((FILAS, COLUMNAS), 0.7))
    T = 250 + 50 * rng.random(FILAS * COLUMNAS)
    conv = (L @ T).reshape(FILAS, COLUMNAS)
    peso = np.cos(np.radians(LATITUDES_GRADOS))
    assert abs((conv * peso[:, None]).sum()) < 1e-9 * np.abs(conv * peso[:, None]).sum()
    bordes, tr = M30.transporte_meridional(conv, 6.371e6)
    assert len(bordes) == FILAS - 1 and bordes[0] == 85 and bordes[-1] == -85


def test_transporte_hacia_el_polo_con_ecuador_caliente():
    from fase2_difusion import construir_matriz_difusion
    L = construir_matriz_difusion(np.full((FILAS, COLUMNAS), 0.7))
    T = (300 - 40 * np.sin(np.radians(LATITUDES_GRADOS)) ** 2)[:, None] * np.ones((1, COLUMNAS))
    _, tr = M30.transporte_meridional((L @ T.ravel()).reshape(FILAS, COLUMNAS), 6.371e6)
    assert tr[5] > 0 and tr[-6] < 0                      # hacia el norte en el norte, hacia el sur en el sur


def test_i10_apagado_por_defecto():
    assert F.INTERRUPTORES_FASE2B["atmosfera_multicapa"] is False


def test_modo_tierra_pone_los_valores_terrestres():
    codigo = ("import parametros as P, temperatura as T, math;"
              "print(P.MODO_TIERRA, P.ROTACION_PERIODO, T.PASO_TIEMPO, round(P.FACTOR_ROTACION_D, 5), "
              "round(P.S3N_LUMINOSIDAD / (4 * math.pi * P.SEMIEJE_MAYOR ** 2), 3), P.P3N_GRAVEDAD, P.INCLINACION_AXIAL)")
    env = dict(os.environ, M3N_MODO="tierra")
    salida = subprocess.run([sys.executable, "-c", codigo], env=env, capture_output=True, text=True,
                            cwd=os.path.dirname(os.path.abspath(__file__)), check=True).stdout.split()
    assert salida == ["True", "86400", "900", "1.0", "1361.0", "9.80665", "23.44"]


def test_sin_modo_tierra_es_p3n():
    import parametros as P
    if os.environ.get("M3N_MODO", "").lower() == "tierra":
        pytest.skip("lanzado en modo Tierra")
    assert P.MODO_TIERRA is False and P.ROTACION_PERIODO == 71424


def test_tropopausa_omm_perfil_estandar():
    # perfil tipo atmosfera estandar: 6,5 K/km hasta 11 km y luego isotermo -> tropopausa a 11 km
    import validar_v30 as V
    zz = np.arange(0, 30001, 500.0)
    Tz = np.where(zz <= 11000, 288.15 - 0.0065 * zz, 288.15 - 0.0065 * 11000)
    pz = 1e5 * np.exp(-zz / 7500)
    zt, Tt = V.tropopausa_omm(Tz, zz, pz)
    assert zt == 11000 and Tt == pytest.approx(216.65)


def test_alturas_recuperan_el_gradiente_critico():
    import validar_v30 as V
    c = columna((0.0,))
    T = c.perfil_inicial(np.full(c.forma, 300.0))[:, 0, 0]
    z = V.alturas(T, c.pm[:, 0, 0], c.ps[0, 0], c.g)
    caliente = T > 200.5                                   # por encima del suelo de 200 K del perfil inicial
    g = -np.diff(T[caliente][::-1]) / np.diff(z[caliente][::-1]) * 1000
    assert np.allclose(g, M30.GRADIENTE_CRITICO * 1000, atol=0.01)

# test_fase0.py -- Fase 0: suite de pytest.
#
# Reune, como pruebas automatizadas de verdad, las comprobaciones que se
# fueron haciendo "a mano" (con print()) al construir cada pieza de la
# Fase 0. Ejecutar con:
#   pytest test_fase0.py -v              (todo)
#   pytest test_fase0.py -v -m "not slow"  (todo menos la simulacion completa)

import math
import pytest

from geometria import declinacion_solar, angulo_horario, angulo_cenital
from atmosfera import masa_aire, trans
from radiacion import i_toa, i_inst, i_atm, i_abs
from temperatura import precalcular_orbita, simular
from orbita import distancia_S3N
from parametros import (S3N_LUMINOSIDAD, ALBEDO, PROFUNDIDAD_OPTICA, SEMIEJE_MAYOR,
                         EMISIVIDAD, INERCIA_TERMICA, INCLINACION_AXIAL_RAD)

from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS, FILAS, COLUMNAS
from fase0_geometria import angulo_cenital_rejilla
from fase0_atmosfera import masa_aire_rejilla
from fase0_radiacion import irradiancia_absorbida_rejilla
from fase0_temperatura import simular_rejilla

MOMENTO_PRUEBA = 0.7
INCLINACION_PRUEBA = math.radians(1.7)
HORA_PRUEBA = 14.5


@pytest.fixture(scope="module")
def declinacion_prueba():
    return declinacion_solar(MOMENTO_PRUEBA, INCLINACION_PRUEBA)


def test_rejilla_tiene_las_dimensiones_esperadas():
    assert LATITUDES_GRADOS.shape == (FILAS,)
    assert LONGITUDES_GRADOS.shape == (COLUMNAS,)


def test_ecuador_cae_en_la_frontera_entre_dos_filas():
    assert LATITUDES_GRADOS[FILAS // 2 - 1] == pytest.approx(2.5)
    assert LATITUDES_GRADOS[FILAS // 2] == pytest.approx(-2.5)


def test_angulo_cenital_rejilla_coincide_con_punto(declinacion_prueba):
    angulos_horarios = angulo_horario(HORA_PRUEBA, LONGITUDES_GRADOS)
    resultado = angulo_cenital_rejilla(LATITUDES_GRADOS, declinacion_prueba, angulos_horarios)

    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h = angulo_horario(HORA_PRUEBA, lon)
            esperado = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h)
            assert resultado[i, j] == pytest.approx(esperado, abs=1e-9)


def test_masa_aire_rejilla_cubre_dia_y_noche(declinacion_prueba):
    angulos_horarios = angulo_horario(HORA_PRUEBA, LONGITUDES_GRADOS)
    cenital = angulo_cenital_rejilla(LATITUDES_GRADOS, declinacion_prueba, angulos_horarios)
    masa = masa_aire_rejilla(cenital)

    hay_dia, hay_noche = False, False
    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h = angulo_horario(HORA_PRUEBA, lon)
            cenital_punto = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h)
            masa_punto = masa_aire(cenital_punto)
            if masa_punto is None:
                hay_noche = True
                assert math.isnan(masa[i, j])
            else:
                hay_dia = True
                assert masa[i, j] == pytest.approx(masa_punto, abs=1e-9)

    assert hay_dia and hay_noche


def test_irradiancia_absorbida_rejilla_coincide_con_punto(declinacion_prueba):
    distancia = distancia_S3N(MOMENTO_PRUEBA, SEMIEJE_MAYOR)
    resultado = irradiancia_absorbida_rejilla(
        distancia, S3N_LUMINOSIDAD, declinacion_prueba, HORA_PRUEBA, ALBEDO, PROFUNDIDAD_OPTICA
    )

    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h = angulo_horario(HORA_PRUEBA, lon)
            cenital_punto = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h)
            toa_punto = i_toa(distancia, S3N_LUMINOSIDAD)
            inst_punto = i_inst(toa_punto, cenital_punto)
            masa_punto = masa_aire(cenital_punto)
            if masa_punto is None:
                esperado = 0.0
            else:
                tra_punto = trans(masa_punto, PROFUNDIDAD_OPTICA)
                atm_punto = i_atm(inst_punto, tra_punto)
                esperado = i_abs(atm_punto, ALBEDO)
            assert resultado[i, j] == pytest.approx(esperado, abs=1e-6)


@pytest.mark.slow
def test_simular_rejilla_coincide_con_punto_en_celdas_representativas():
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    T_grid, _ = simular_rejilla(datos_orbita, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)

    celdas_de_prueba = [(0, 0), (9, 18), (17, 36), (18, 54), (35, 71)]
    for i, j in celdas_de_prueba:
        lat_rad = math.radians(LATITUDES_GRADOS[i])
        lon_rad = math.radians(LONGITUDES_GRADOS[j])
        datos_desplazados = [(toa, decl, ang_h + lon_rad) for (toa, decl, ang_h) in datos_orbita]
        T_punto, _ = simular(lat_rad, datos_desplazados, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        assert T_grid[i, j] == pytest.approx(T_punto, abs=0.02)

# test_fase52_lineas.py -- pruebas del calculo linea a linea (Fase 5.2, borrador). Casos con solucion conocida.
import math

import numpy as np

import fase52_lineas as L


def test_voigt_area_uno_y_limites():
    """Area 1; con Doppler despreciable es la de Lorentz y con presion despreciable la de Gauss (valor en el pico)."""
    nu = np.linspace(-50, 50, 400001)
    dnu = nu[1] - nu[0]
    for aD, gL in ((0.01, 0.05), (0.05, 0.01), (0.03, 0.03)):
        assert abs(L.voigt(nu, 0.0, aD, gL).sum() * dnu - 1) < 2e-3      # colas fuera de +-50 (Lorentz) < 0,1 %
    gL = 0.07                                                              # Lorentz: pico 1/(pi gamma)
    assert abs(L.voigt(np.array([0.0]), 0.0, 1e-7, gL)[0] * math.pi * gL - 1) < 1e-6
    aD = 0.02                                                              # Gauss: pico sqrt(ln2/pi)/alfa
    assert abs(L.voigt(np.array([0.0]), 0.0, aD, 1e-9)[0] / (math.sqrt(math.log(2) / math.pi) / aD) - 1) < 1e-6
    # semianchura de Gauss: a una distancia alfa_D del centro, la mitad del pico
    assert abs(L.voigt(np.array([aD]), 0.0, aD, 1e-9)[0] / L.voigt(np.array([0.0]), 0.0, aD, 1e-9)[0] - 0.5) < 1e-6


def test_voigt_anchura_olivero_longbothum():
    """La anchura total del Voigt numerico coincide con la aproximacion de Olivero y Longbothum (1977, precision
    ~0,02 %): f_V = 0,5346 f_L + sqrt(0,2166 f_L^2 + f_G^2)."""
    aD, gL = 0.03, 0.02
    nu = np.linspace(0, 0.2, 200001)
    v = L.voigt(nu, 0.0, aD, gL)
    mitad = nu[np.argmin(np.abs(v - v[0] / 2))]
    fL, fG = 2 * gL, 2 * aD
    fV = 0.5346 * fL + math.sqrt(0.2166 * fL ** 2 + fG ** 2)
    assert abs(2 * mitad / fV - 1) < 5e-4


def test_intensidad_igual_en_la_referencia_y_boltzmann():
    """A T_ref no cambia; el cociente entre dos temperaturas sigue la ley de Boltzmann con la energia del nivel bajo."""
    S = L.intensidad(1e-20, 1500.0, 200.0, L.T_REF, 174.58, 174.58)
    assert abs(S / 1e-20 - 1) < 1e-14
    T = 250.0
    r = L.intensidad(1.0, 1500.0, 1000.0, T, 1.0, 1.0) / L.intensidad(1.0, 1500.0, 0.0, T, 1.0, 1.0)
    assert abs(r / (math.exp(-L.C2 * 1000 / T) / math.exp(-L.C2 * 1000 / L.T_REF)) - 1) < 1e-12


def test_doppler_del_co2():
    """HWHM Doppler del CO2 (44 u) a 667 cm^-1 y 296 K: ~5,9e-4 cm^-1 (calculo directo de la formula)."""
    a = L.hwhm_doppler(667.0, 43.9898, 296.0)
    esperado = 667.0 / 2.99792458e8 * math.sqrt(2 * 1.380649e-23 * 296 * math.log(2) / (43.9898 * 1.66053906660e-27))
    assert abs(a / esperado - 1) < 1e-12 and 5.5e-4 < a < 6.3e-4


def test_coeficiente_conserva_la_intensidad():
    """Integrado en numero de onda, k devuelve la suma de las intensidades (salvo las alas cortadas)."""
    lineas = {"nu0": np.array([1000.0, 1001.3]), "S": np.array([1e-20, 3e-21]), "E": np.array([0.0, 100.0]),
              "g_air": np.array([0.07, 0.08]), "g_self": np.array([0.3, 0.35]), "n_air": np.array([0.7, 0.75]),
              "d_air": np.array([0.0, -0.002])}
    nu = np.linspace(970, 1031, 610001)
    k = L.coeficiente_absorcion(nu, lineas, L.T_REF, 1.0, 0.01, 1.0, 1.0, 18.0)
    total = k.sum() * (nu[1] - nu[0])
    # con corte a 25 cm^-1 y anchura ~0,07, lo que queda fuera es ~2 gamma/(pi 25) ~ 0,2 %
    assert abs(total / 1.3e-20 - 1) < 4e-3

# test_fase6_superficie.py -- superficie con viento real (Fase 6.3, paso 5; sin conectar). pytest.
import numpy as np
import pytest

import fase2b_atmosfera as F2
import fase6_superficie as S


def test_rugosidad_oceano_formula_del_ecmwf():
    us = np.array([0.1, 0.3, 0.6])
    z0, z0h, z0q = S.rugosidad_oceano(us, 9.81)
    assert np.allclose(z0, np.minimum(0.11 * 1.5e-5 / us + 0.018 * us ** 2 / 9.81, 1e-3), rtol=1e-13, atol=0)
    assert np.allclose(z0h, 0.40 * 1.5e-5 / us, rtol=1e-13, atol=0) and np.allclose(z0q, 0.62 * 1.5e-5 / us, rtol=1e-13, atol=0)


def test_coeficientes_punto_fijo_y_valores_de_referencia():
    U = np.array([3.0, 5.0, 10.0, 15.0])
    d = S.coeficientes(10.0, U, np.full(4, 290.0), np.full(4, 290.0), np.ones(4, bool), 0.01, 0.001, 9.81)
    assert d["cambio_final"] < 1e-11
    # z0 es solucion de su ecuacion con el u* resultante
    z0, z0h, _ = S.rugosidad_oceano(d["u_estrella"], 9.81)
    assert np.allclose(z0, d["z0m"], rtol=1e-10) and np.allclose(z0h, d["z0h"], rtol=1e-10)
    # C_H neutro a 10 m: ~1,1e-3 como el valor fijo de M3N (Large y Pond 1982); medido 1,06-1,21e-3
    assert np.all((d["c_h"] > 1.0e-3) & (d["c_h"] < 1.25e-3))
    # C_D neutro a 10 m crece con el viento (Charnock)
    assert np.all(np.diff(d["c_m"]) > 0)


def test_coeficientes_sobre_tierra_y_estabilidad():
    d = S.coeficientes(185.0, np.array([5.0, 5.0, 5.0]), np.array([290.0, 290.0, 290.0]),
                       np.array([287.0, 290.0, 293.0]), np.zeros(3, bool), 0.01, 0.001, 9.81)
    assert np.all(d["z0m"] == 0.01)                                   # en tierra no cambia
    assert d["c_h"][0] < d["c_h"][1] < d["c_h"][2]                     # estable < neutro < inestable


@pytest.mark.slow
def test_albedo_con_5_m_s_identico_a_la_v31_y_dependencia_del_viento():
    mu = np.linspace(0.01, 1.0, 50); tr = np.full(50, 0.8)
    a = F2.albedo_oceano(mu, tr)
    b = S.albedo_oceano_viento(mu, tr, np.full(50, 5.0))
    assert np.array_equal(a, b)
    bajo = [float(S.albedo_oceano_viento(np.array([0.1]), np.array([0.8]), np.array([U]))[0]) for U in (1.0, 5.0, 15.0)]
    assert bajo[0] > bajo[1] > bajo[2]                                 # con el sol bajo, el mar en calma refleja mas

# test_fase52_ozono.py -- el ozono prescrito (Fase 5.2, propuesta del §3.3).
import numpy as np
import pytest

import fase52_ozono as O


def test_la_columna_de_la_atmosfera_estandar():
    """Integrando en presion el perfil de la atmosfera estandar con g = 9,80665 sale su columna, ~345,7 DU (la misma
    integrada en altura con la densidad de la tabla; < 1 % de diferencia por la discretizacion)."""
    ph = np.concatenate([[1e-5], np.geomspace(1e-4, 1013.25, 4000)])
    pm = np.sqrt(ph[1:] * ph[:-1])
    assert O.columna_DU(ph, O.forma(pm), 9.80665) == pytest.approx(345.7, rel=0.01)


def test_perfil_escalado_a_la_columna_pedida_y_con_la_forma():
    ph = np.concatenate([[0.0], np.geomspace(0.5, 1000.0, 40)])
    v = O.perfil(ph, 9.05, 180.0)
    assert O.columna_DU(ph, v, 9.05) == pytest.approx(180.0, rel=1e-12)
    pm = np.sqrt(ph[1:] * np.maximum(ph[:-1], 1e-6))
    r = v / O.forma(pm)
    assert np.allclose(r, r[0], rtol=1e-12)                 # misma forma: solo un factor
    assert pm[np.argmax(v)] == pytest.approx(5.7, rel=0.3)   # maximo hacia 5-6 hPa


def test_cociente_de_rugheimer():
    assert O.COCIENTE_RUGHEIMER == pytest.approx(0.604, abs=1e-3)

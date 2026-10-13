# test_fase52_nubes_solape.py -- solapamiento de nubes y subcolumnas McICA (Fase 5.2, T23).
import numpy as np
import pytest

import fase52_nubes_solape as S


def test_limites_aleatorio_y_maximo():
    f = np.array([0.0, 0.3, 0.5, 0.0, 0.4, 0.2])
    # aleatorio (alfa = 0): 1 - prod(1 - f)
    assert S.cobertura_total(f, np.zeros(5)) == pytest.approx(1 - np.prod(1 - f), rel=1e-12)
    # maximo-aleatorio (alfa = 1): maximo dentro de cada bloque seguido, aleatorio entre bloques
    esperado = 1 - (1 - 0.5) * (1 - 0.4)
    assert S.cobertura_total(f, np.ones(5)) == pytest.approx(esperado, rel=1e-12)
    # una capa cubierta del todo: cobertura 1
    assert S.cobertura_total(np.array([0.2, 1.0, 0.3]), np.full(2, 0.5)) == pytest.approx(1.0)


def test_cobertura_entre_los_limites_y_creciente():
    f = np.array([0.1, 0.4, 0.4, 0.2, 0.0, 0.3])
    a = [S.cobertura_total(f, np.full(5, x)) for x in (0.0, 0.3, 0.7, 1.0)]
    assert all(b <= c + 1e-15 for b, c in zip(a[1:], a[:-1]))            # mas solapamiento, menos cobertura
    acum, _ = S.cobertura_acumulada(f, np.full(5, 0.5))
    assert np.all(np.diff(acum) >= -1e-15) and acum[-1] <= 1.0


def test_longitud_de_decorrelacion_del_ifs():
    assert S.longitud_decorrelacion_km(0.0) == pytest.approx(2.899)
    assert S.longitud_decorrelacion_km(90.0) == pytest.approx(0.75)
    a = S.alfa_capas(np.array([10.0, 9.0, 7.0]), 2.0)
    assert np.allclose(a, np.exp([-0.5, -1.0]))


@pytest.mark.parametrize("alfa", [0.0, 0.6, 0.95, 1.0])
def test_subcolumnas_reproducen_las_fracciones_y_las_coberturas(alfa):
    """Con muchas subcolumnas: la frecuencia de nube en cada capa es f, la de cada par de capas seguidas es la
    cobertura del par y la de la columna es la cobertura total (dentro de 4 desviaciones del muestreo)."""
    f = np.array([0.0, 0.15, 0.5, 0.35, 0.0, 0.25, 0.6])
    a = np.full(len(f) - 1, alfa)
    rng = np.random.default_rng(42)
    n = 40000
    sc = S.subcolumnas(f, a, n, rng)
    tol = lambda p: 4 * np.sqrt(max(p * (1 - p), 1e-6) / n)
    for k in range(len(f)):
        assert abs(sc[:, k].mean() - f[k]) < tol(f[k])
    _, par = S.cobertura_acumulada(f, a)
    for k in range(len(f) - 1):
        frec = (sc[:, k] | sc[:, k + 1]).mean()
        assert abs(frec - par[k]) < tol(par[k])
    C = S.cobertura_total(f, a)
    assert abs(sc.any(axis=1).mean() - C) < tol(C)


def test_sin_nubes():
    sc = S.subcolumnas(np.zeros(5), np.full(4, 0.5), 100, np.random.default_rng(1))
    assert not sc.any()

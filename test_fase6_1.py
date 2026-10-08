# test_fase6_1.py -- pruebas automaticas de la etapa 6.1 (aguas someras). Ejecutar con pytest.
import math
import warnings
import numpy as np
import pytest

from fase6_aguas_someras import (AguasSomeras, A_TIERRA, OMEGA_TIERRA, G_TIERRA, DIA, caso1, caso2, caso5,
                                 caso6, errores_normalizados, RAW_NU, RAW_ALFA)

warnings.filterwarnings("ignore", category=RuntimeWarning)


def _m(**kw):
    return AguasSomeras(A_TIERRA, OMEGA_TIERRA, G_TIERRA, **kw)


def _estado_ruidoso():
    rng = np.random.default_rng(1)
    h, u, v = caso6()
    return (h + 50 * rng.standard_normal(h.shape), u + 5 * rng.standard_normal(u.shape),
            v + 5 * rng.standard_normal(v.shape))


def test_masa_se_conserva_exactamente():
    h, u, v = _estado_ruidoso()
    _, _, _, hs = caso5()
    m = _m(hs=hs)
    dh, _, _ = m.tendencias(h, u, v)
    assert abs((dh * m.R.area).sum()) / abs(dh * m.R.area).sum() < 1e-14


def test_energia_se_conserva_en_el_sistema_semidiscreto():
    """Sin filtro polar, la tendencia de la energia total es cero salvo redondeo (DISENO_FASE_6.1.md §3)."""
    h, u, v = _estado_ruidoso()
    _, _, _, hs = caso5()
    m = _m(hs=hs, lat_filtro=90.0)
    dh, du, dv = m.tendencias(h, u, v)
    R = m.R
    B = m.g * (h + hs) + m.energia_cinetica(u, v)
    h_u = 0.5 * (h + np.roll(h, 1, axis=1)); h_v = 0.5 * (h[:-1] + h[1:])
    t = [(B * dh * R.area).sum(), (h_u * u * du * R.A_u).sum(), (h_v * v * dv * R.A_v).sum()]
    assert abs(sum(t)) / max(map(abs, t)) < 1e-12


def test_reposo_se_mantiene():
    h = np.full((36, 72), 5000.0); u = np.zeros((36, 72)); v = np.zeros((35, 72))
    dh, du, dv = _m().tendencias(h, u, v)
    assert abs(dh).max() == 0 and abs(du).max() < 1e-15 and abs(dv).max() < 1e-15


def test_caso2_un_dia_y_eje_girado():
    for alfa in (0.0, math.pi / 4, math.pi / 2):
        h, u, v = caso2(alfa)
        m = _m(alfa_coriolis=alfa)
        fin = m.integrar(h, u, v, 150.0, int(DIA / 150))
        assert errores_normalizados(fin[0], h, m.R.area)[1] < 2e-3


def test_caso2_converge_con_orden_2():
    e = []
    for filas in (36, 72):
        h, u, v = caso2(math.pi / 4, filas=filas)
        m = _m(alfa_coriolis=math.pi / 4, filas=filas)
        dt = 150.0 * 36 / filas
        fin = m.integrar(h, u, v, dt, int(round(DIA / dt)))
        e.append(errores_normalizados(fin[0], h, m.R.area)[1])
    assert math.log2(e[0] / e[1]) > 1.7


def test_altura_rossby_haurwitz_esta_en_equilibrio():
    """Si la formula de h del caso 6 esta bien transcrita, el desequilibrio (divergencia de la tendencia
    del momento) tiende a cero con la resolucion a orden 2."""
    r = []
    for filas in (36, 72):
        h, u, v = caso6(filas=filas)
        m = _m(lat_filtro=90.0, filas=filas)
        R = m.R
        _, du, dv = m.tendencias(h, u, v)
        Fu, Fv = du * R.L_u, dv * R.L_v
        div = np.roll(Fu, -1, axis=1) - Fu
        div[:-1] -= Fv; div[1:] += Fv; div /= R.area
        rot = m.vorticidad(du, dv)
        k = filas // 9
        r.append(math.sqrt((div[k:-k] ** 2).mean() / (rot[k:-k] ** 2).mean()))
    assert r[0] < 0.25 and r[0] / r[1] > 3.5


def test_filtro_polar_no_toca_la_media_zonal():
    from fase6_aguas_someras import _filtro_fourier
    rng = np.random.default_rng(2)
    x = rng.standard_normal((36, 72))
    R = _m().R
    y = _filtro_fourier(x, R.cos_c, math.cos(math.radians(60)), R.dlam)
    assert np.allclose(x.mean(1), y.mean(1), atol=1e-14)
    assert np.array_equal(x[12:24], y[12:24])          # fuera de las latitudes altas no cambia nada


def _amplificacion_raw(wdt, nu=RAW_NU, al=RAW_ALFA):
    M = np.zeros((2, 2), complex)
    lam = 1j * wdt
    for j, (xa, x) in enumerate(((1, 0), (0, 1))):
        xn = xa + 2 * lam * x
        d = nu / 2 * (xa - 2 * x + xn)
        M[:, j] = [x + al * d, xn + (al - 1) * d]
    return abs(np.linalg.eigvals(M)).max()


def test_limite_de_estabilidad_de_raw():
    """Propiedad del esquema (DISENO_FASE_6.1.md §5): leapfrog + RAW (0,53; 0,2) es estable hasta
    omega*dt ~ 0,45 y amplifica por encima (0,7 -> 1,0056 por paso)."""
    assert _amplificacion_raw(0.4) <= 1.0
    assert _amplificacion_raw(0.7) > 1.005


def test_caso6_estable_14_dias():
    h, u, v = caso6()
    m = _m()
    I0 = m.integrales(h, u, v)
    fin = m.integrar(h, u, v, 150.0, int(14 * DIA / 150))
    I = m.integrales(*fin)
    assert np.isfinite(fin[0]).all()
    assert abs(I["masa"] / I0["masa"] - 1) < 1e-13
    assert abs(I["energia"] / I0["energia"] - 1) < 1e-5

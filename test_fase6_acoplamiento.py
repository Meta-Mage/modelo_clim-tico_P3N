# test_fase6_acoplamiento.py -- primeras piezas del acoplamiento fisica-nucleo (Fase 6.3, paso 5). pytest.
import numpy as np

import fase30_multicapa as Mm
import fase31_agua as AG
from fase6_acoplamiento import viento_en_centros, tendencia_a_caras, calor_rozamiento_exacto
from fase6_nucleo import NucleoSeco
from fase6_forzamientos import sigma_capas
from fase6_aguas_someras import A_TIERRA, OMEGA_TIERRA, G_TIERRA

F, C, N = 36, 72, 20


def _alt(semilla=3):
    rng = np.random.default_rng(semilla)
    return np.maximum(0.0, 3000 * rng.random((F, C)) - 1000)


def test_columna_actualizar_ps_reproduce_la_construccion():
    col = Mm.Columna(_alt())
    antes = {k: np.array(v, copy=True) for k, v in vars(col).items() if isinstance(v, np.ndarray)}
    col.actualizar_ps(col.ps.copy())
    for k, v in antes.items():
        assert np.array_equal(v.view(np.uint8), np.asarray(getattr(col, k)).view(np.uint8)), k


def test_columna_ps_variable_cambia_lo_que_depende_de_ps():
    col = Mm.Columna(_alt())
    ps = col.ps * (1 + 0.01 * np.sin(np.arange(F * C).reshape(F, C)))
    col.actualizar_ps(ps)
    nueva = Mm.Columna(_alt())
    nueva.ps = None
    nueva.actualizar_ps(ps.copy())
    for k in ("ph", "pm", "dp", "cap", "tau", "trans", "peso_sw", "factor_superficie", "cap_transporte"):
        assert np.array_equal(getattr(col, k), getattr(nueva, k)), k
    assert np.allclose(col.dp.sum(axis=0), ps)                                   # la masa de la columna es p_s


def test_presion_de_capa_sb81_igual_que_el_nucleo():
    col = Mm.Columna(_alt(), presion_capa="sb81")
    m = NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, Mm.sigma_seminiveles())
    assert np.array_equal(col.pm, sigma_capas(m)[:, None, None] * col.ps[None])
    media = Mm.Columna(_alt())
    assert np.array_equal(col.pm[0], media.pm[0])                                 # capa de arriba: iguales
    assert np.max(np.abs(col.pm / media.pm - 1)) < 0.007                          # resto: < 0,7 % (DISENO §6.1)


def test_cache_de_geometria_de_la_conveccion_detecta_cambios_interiores():
    col = Mm.Columna(_alt())
    rng = np.random.default_rng(5)
    T = col.perfil_inicial(np.full((F, C), 295.0)) + rng.standard_normal((N, F, C))
    qs, _ = AG.qs_y_derivada(T, col.pm)
    q = 0.9 * qs
    AG.conveccion_humeda(T, q, col.pm, col.ph, 992.0, 9.8)                        # llena la cache
    ps2 = col.ps.copy(); ps2[10, 20] *= 0.97                                       # cambia SOLO una celda interior
    col.actualizar_ps(ps2)
    a = AG.conveccion_humeda(T, q, col.pm, col.ph, 992.0, 9.8)
    AG._CACHE_GEO.clear()
    b = AG.conveccion_humeda(T, q, col.pm, col.ph, 992.0, 9.8)                    # sin cache
    for x, y in zip(a, b):
        assert np.array_equal(x, y)


def test_viento_en_centros_y_su_traspuesta():
    rng = np.random.default_rng(1)
    u = rng.standard_normal((N, F, C)); v = rng.standard_normal((N, F - 1, C))
    uc, vc = viento_en_centros(u, v)
    a = rng.standard_normal((N, F, C)); b = rng.standard_normal((N, F, C))
    du, dv = tendencia_a_caras(a, b)
    # <P x, y> = <x, P^T y> en el producto escalar sin pesos
    assert abs((uc * a).sum() + (vc * b).sum() - (u * du).sum() - (v * dv).sum()) < 1e-10
    # un viento zonal uniforme se conserva al ir y volver
    uu = np.full((N, F, C), 7.0); vv = np.zeros((N, F - 1, C))
    uc, vc = viento_en_centros(uu, vv)
    assert np.all(uc == 7.0) and np.all(vc == 0.0)
    du, dv = tendencia_a_caras(uc, vc)
    assert np.all(du == 7.0) and np.all(dv == 0.0)


def test_calor_por_rozamiento_cierra_la_energia_del_nucleo():
    m = NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, Mm.sigma_seminiveles())
    rng = np.random.default_rng(2)
    u0 = 10 * rng.standard_normal((N, F, C)); v0 = 10 * rng.standard_normal((N, F - 1, C))
    u1 = 0.9 * u0; v1 = 0.9 * v0
    dT = calor_rozamiento_exacto(m, u0, v0, u1, v1, 1004.0)
    resto = 1004.0 * dT + m.energia_cinetica(u1, v1) - m.energia_cinetica(u0, v0)
    assert np.max(np.abs(resto)) < 1e-10
    assert np.all(dT >= 0)                                                         # frenar calienta

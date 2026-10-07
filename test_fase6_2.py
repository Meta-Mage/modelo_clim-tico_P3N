# test_fase6_2.py -- pruebas automaticas del nucleo seco (Fase 6.2). Ejecutar con pytest.
import math
import numpy as np

import fase30_multicapa as Mm
from fase6_nucleo import NucleoSeco
from fase6_forzamientos import estado_jw, HeldSuarez
from fase6_aguas_someras import caso5, A_TIERRA, OMEGA_TIERRA, G_TIERRA

R_GAS, CP = 287.0, 1004.0
SH = Mm.sigma_seminiveles()


def _nucleo(**kw):
    return NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, R_GAS, CP, kw.pop("sh", SH), **kw)


def _estado_aleatorio(m, hs):
    rng = np.random.default_rng(3)
    N = m.N
    ps = 1e5 * np.exp(-hs / 8000) * (1 + 0.02 * rng.standard_normal((36, 72)))
    T = 250 + 20 * rng.standard_normal((N, 36, 72))
    return ps, T, 10 * rng.standard_normal((N, 36, 72)), 10 * rng.standard_normal((N, 35, 72))


def test_energia_y_masa_exactas_en_el_sistema_semidiscreto():
    _, _, _, hs = caso5()
    m = _nucleo(phis=G_TIERRA * hs, lat_filtro=90.0)
    est = _estado_aleatorio(m, hs)
    tend = m.tendencias(*est)
    dE, partes = m.derivada_energia(*est, tend)
    assert abs(dE) / np.abs(partes).max() < 1e-12
    a = tend[0] * m.Rj.area
    assert abs(a.sum()) / np.abs(a).sum() < 1e-14


def test_reposo_isotermo_sobre_montanas_exacto():
    _, _, _, hs = caso5()
    m = _nucleo(phis=G_TIERRA * hs)
    T0 = 280.0
    ps = 1e5 * np.exp(-G_TIERRA * hs / (R_GAS * T0))
    est = (ps, np.full((m.N, 36, 72), T0), np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)))
    for d in m.tendencias(*est):
        assert np.abs(d).max() < 1e-13


def test_onda_de_lamb():
    """El modo externo del sistema lineal va a ~sqrt(gamma R T) (onda de Lamb), a menos de un 1 %."""
    m = _nucleo()
    _, _, M = m.matrices_verticales(300.0)
    c = math.sqrt(np.linalg.eigvals(M).real.max())
    assert abs(c / math.sqrt(CP / (CP - R_GAS) * R_GAS * 300.0) - 1) < 0.01


def test_equilibrio_jw_converge_con_orden_2_en_horizontal():
    sh = Mm.sigma_seminiveles(160)
    r = []
    for f in (36, 72):
        m0 = _nucleo(sh=sh, filas=f)
        ps, T, u, v, phis = estado_jw(m0)
        m = _nucleo(sh=sh, filas=f, phis=phis)
        d = m.tendencias(ps, T, u, v)
        k = f // 9
        r.append(math.sqrt((d[3][:, k:-k] ** 2).mean()))
    assert r[0] / r[1] > 3.5


def test_semiimplicito_estable_con_paso_largo():
    """Con dt = 900 s (6 veces el limite explicito) el reposo perturbado se mantiene acotado."""
    m = _nucleo()
    rng = np.random.default_rng(0)
    ps = 1e5 * (1 + 1e-3 * rng.standard_normal((36, 72)))
    T = 280.0 + 0.5 * rng.standard_normal((m.N, 36, 72))
    m.preparar_semiimplicito(900.0)
    I0 = m.integrales(ps, T, np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)))
    fin = m.integrar_si(ps, T, np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)), 96)
    I = m.integrales(*fin)
    assert np.isfinite(fin[2]).all() and np.abs(fin[2]).max() < 20
    assert abs(I["masa"] / I0["masa"] - 1) < 1e-14
    assert abs(I["energia"] / I0["energia"] - 1) < 1e-5


def test_held_suarez_valores():
    m = _nucleo()
    hs = HeldSuarez(m)
    ps = np.full((36, 72), 1e5)
    Teq = hs.T_equilibrio(ps)
    assert Teq.min() >= 200.0 and abs(Teq[-1].max() - 315 * 0.99 ** m.kappa) < 3.0
    assert hs.kv[-1, 0, 0] > 0 and hs.kv[0, 0, 0] == 0


def test_hiperdifusion_matrices_exactas_y_conservadoras():
    m = _nucleo()
    m.preparar_hiperdifusion(0.5)
    rng = np.random.default_rng(0)
    T = rng.standard_normal((2, 36, 72)); u = rng.standard_normal((2, 36, 72)); v = rng.standard_normal((2, 35, 72))
    Rj = m.Rj
    K = np.arange(37)
    fu = np.exp(1j * K * Rj.lam_u[0]); fh = np.exp(1j * K * Rj.lam_c[0])
    w = np.concatenate([np.fft.rfft(u, axis=-1) / fu, np.fft.rfft(v, axis=-1) / fh], axis=-2)
    wl = np.einsum("kfg,ngk->nfk", m._Lvec, w)
    Lu, Lv = m._laplaciano_vector(u, v)
    assert abs(np.fft.irfft(wl[:, :36] * fu, n=72, axis=-1) - Lu).max() < 1e-13 * abs(Lu).max()
    assert abs(np.fft.irfft(wl[:, 36:] * fh, n=72, axis=-1) - Lv).max() < 1e-13 * abs(Lv).max()
    x, _, _ = m.aplicar_hiperdifusion(np.full((1, 36, 72), 250.0), np.zeros((1, 36, 72)), np.zeros((1, 35, 72)), 900.0)
    assert abs(x - 250).max() < 1e-10
    # la media global (ponderada por areas) de T no cambia
    y, _, _ = m.aplicar_hiperdifusion(T, u, v, 900.0)
    assert abs((y * Rj.area).sum() - (T * Rj.area).sum()) < 1e-9 * abs(T * Rj.area).sum()


def test_punto_de_control_identico():
    """Cortar y reanudar da exactamente lo mismo que una corrida sin cortes."""
    m = _nucleo(sh=np.linspace(0, 1, 21))
    hs = HeldSuarez(m)
    m.preparar_semiimplicito(450.0); m.preparar_hiperdifusion(0.5)
    ps = np.full((36, 72), 1e5); T = hs.T_equilibrio(ps) + 0.1 * np.random.default_rng(1).standard_normal((m.N, 36, 72))
    z = (np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)))
    seguido = m.integrar_si(ps, T, *z, 12, forzamiento=hs)
    guardado = {}
    m.integrar_si(ps, T, *z, 6, forzamiento=hs, guardar=lambda n, a, b: guardado.update(n=n, a=a, b=b), cada_guardar=6)
    reanudado = m.integrar_si(None, None, None, None, 12, forzamiento=hs, reanudar=(guardado["n"], guardado["a"], guardado["b"]))
    for x, y in zip(seguido, reanudado):
        assert np.array_equal(x, y)


# ---------------- v3.1-pre8: version compilada (numba) y filtros en lote: IDENTICAS bit a bit ----------------

def _bits(x):
    return np.ascontiguousarray(x).view(np.int64)


def test_tendencias_compiladas_identicas_bit_a_bit():
    from fase6_nucleo_nb import HAY_NUMBA
    if not HAY_NUMBA:
        return                                                   # sin numba se usa la version de numpy
    _, _, _, hs = caso5()
    m = _nucleo(phis=G_TIERRA * hs)
    est = _estado_aleatorio(m, hs)
    m.usar_numba = True
    a = m.tendencias(*est)
    b = m._tendencias_numpy(*est)
    for x, y in zip(a, b):
        assert np.array_equal(_bits(x), _bits(y))


def test_filtro_en_lote_identico_bit_a_bit():
    _, _, _, hs = caso5()
    m = _nucleo(phis=G_TIERRA * hs)
    ps, T, u, v = _estado_aleatorio(m, hs)
    a = m._filtro_varios((ps, T, u), m.Rj.cos_c)
    for x, y in zip(a, (m._filtro3(ps, m.Rj.cos_c), m._filtro3(T, m.Rj.cos_c), m._filtro3(u, m.Rj.cos_c))):
        assert np.array_equal(_bits(x), _bits(y))


def test_held_suarez_20_pasos_numba_igual_que_numpy():
    from fase6_nucleo_nb import HAY_NUMBA
    if not HAY_NUMBA:
        return
    def corre(nb):
        m = _nucleo()
        m.usar_numba = nb
        hsf = HeldSuarez(m)
        m.preparar_semiimplicito(450.0)
        m.preparar_hiperdifusion(0.5)
        rng = np.random.default_rng(1)
        ps = np.full((36, 72), 1e5)
        T = hsf.T_equilibrio(ps) + 0.1 * rng.standard_normal((m.N, 36, 72))
        return m.integrar_si(ps, T, np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)), 20, forzamiento=hsf)
    for x, y in zip(corre(False), corre(True)):
        assert np.array_equal(_bits(x), _bits(y))


# ---------------- v3.1-pre10: opciones de la hiperdifusion (apagadas por defecto) ----------------

def test_hiperdifusion_calor_cierra_la_energia():
    m = _nucleo()
    rng = np.random.default_rng(4)
    ps = 1e5 * (1 + 0.01 * rng.standard_normal((36, 72)))
    T = 250 + 10 * rng.standard_normal((m.N, 36, 72))
    u = 10 * rng.standard_normal((m.N, 36, 72)); v = 10 * rng.standard_normal((m.N, 35, 72))
    resultados = {}
    for calor in (False, True):
        m.preparar_hiperdifusion(0.5, calor_rozamiento=calor)
        e0 = m.integrales(ps, T, u, v)["energia"]
        Tn, un, vn = m.aplicar_hiperdifusion(T, u, v, 900.0, ps=ps)
        resultados[calor] = abs(m.integrales(ps, Tn, un, vn)["energia"] - e0) / e0
    # sin calor se pierde la energia cinetica quitada; con calor el resto es >= 100 veces menor (DISENO §6.7)
    assert resultados[True] < resultados[False] / 100


def test_hiperdifusion_correccion_presion_reduce_el_calentamiento_falso_sobre_montanas():
    _, _, _, hs = caso5()
    m = _nucleo(phis=G_TIERRA * hs)
    ps = 1e5 * np.exp(-hs / 8000.0)
    from fase6_forzamientos import sigma_capas
    T = 288.0 * (sigma_capas(m)[:, None, None] * ps[None] / 1e5) ** 0.19       # T solo depende de p
    z = (np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)))
    dT = {}
    for corr in (False, True):
        m.preparar_hiperdifusion(0.5, correccion_presion=corr)
        Tn, _, _ = m.aplicar_hiperdifusion(T, *z, 900.0, ps=ps)
        dT[corr] = np.abs(Tn - T).max()
    # medido el 07/10: 3,4e-2 K -> 8,3e-3 K por paso de 900 s con la forma de CAM (x4,1); con la interpolacion
    # vertical de la v3.1-pre14 (DISENO_FASE6_3.md §6.12), 5,0e-3 K (x6,8)
    assert dT[True] < dT[False] / 3.5

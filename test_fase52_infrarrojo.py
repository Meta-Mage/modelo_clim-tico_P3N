# test_fase52_infrarrojo.py -- el infrarrojo por bandas (Fase 5.2, borrador): Planck por bandas (T20), el dos flujos
# (T5, T19) y el solapamiento de gases (T9).
import numpy as np
import pytest

import fase52_hitran as H
import fase52_infrarrojo as F
import fase52_k as K


@pytest.mark.parametrize("T", [150.0, 200.0, 255.0, 300.0, 330.0])
def test_planck_bandas_suman_sigma_T4(T):
    """Con las colas (T20), las 16 bandas suman sigma T^4 (la diferencia, 1,1e-9, es la de redondeo entre las
    constantes C1, C2 y sigma de CODATA); fuera de 10-3250 cm^-1 queda < 0,05 % (medido: 0,039 % a 330 K)."""
    b = F.planck_banda(H.BANDAS_IR, T)
    assert b.sum() == pytest.approx(F.SIGMA * T ** 4, rel=2e-9)
    fuera = F.planck_integral(1e-6, 10.0, T) + F.planck_cola_alta(3250.0, T)
    assert fuera / (F.SIGMA * T ** 4) < 5e-4
    sin = F.planck_banda(H.BANDAS_IR, T, colas=False)
    assert sin.sum() + fuera == pytest.approx(F.SIGMA * T ** 4, rel=2e-9)


def test_planck_bandas_tienen_que_ir_seguidas():
    with pytest.raises(ValueError):
        F.planck_banda(((10, 350), (400, 500)), 250.0)


def test_g1_serie_y_formula_coinciden():
    for x in (0.0099, 0.01, 0.0101):
        serie = x / 2 - x ** 2 / 3 + x ** 3 / 8 - x ** 4 / 30 + x ** 5 / 144
        directa = (1.0 - np.exp(-x) * (1.0 + x)) / x
        assert serie == pytest.approx(directa, rel=1e-9)
    assert F._g1(np.array([1e-12]))[0] == pytest.approx(5e-13, rel=1e-6)


def _perfil(n=12, semilla=1):
    r = np.random.default_rng(semilla)
    dtau = r.uniform(0.01, 2.0, n)
    Bh = np.linspace(80.0, 390.0, n + 1) + r.uniform(-5, 5, n + 1)
    return dtau, Bh


def test_isotermo_exacto():
    """Atmosfera isoterma sobre un suelo negro a la misma temperatura: U = B en todas partes y
    D = B (1 - e^(-1,66 tau desde arriba))."""
    dtau = np.array([0.1, 0.5, 2.0, 0.03])
    B = 300.0
    U, D = F.dos_flujos(dtau, [1.0], np.full(5, B), B)
    assert np.allclose(U, B, rtol=1e-13)
    tau_arriba = np.concatenate([[0.0], np.cumsum(dtau)])
    assert np.allclose(D, B * (1 - np.exp(-F.DIFUSIVIDAD * tau_arriba)), rtol=1e-13)


def test_transparente_y_suelo_gris():
    dtau, Bh = _perfil()
    U, D = F.dos_flujos(np.zeros_like(dtau), [1.0], Bh, 400.0, emisividad=0.9)
    assert np.allclose(D, 0.0) and np.allclose(U, 0.9 * 400.0)
    U, D = F.dos_flujos(dtau, [1.0], Bh, 400.0, emisividad=0.9)
    assert U[-1] == pytest.approx(0.9 * 400.0 + 0.1 * D[-1], rel=1e-14)


def test_fuente_lineal_en_tau_exacta_al_partir_capas():
    """La fuente de Planck lineal en tau es exacta si B es lineal en tau dentro de cada capa: partir cada capa en 7
    subcapas con B interpolado linealmente en tau no cambia los flujos (10^-12)."""
    dtau, Bh = _perfil()
    U, D = F.dos_flujos(dtau, [1.0], Bh, 395.0)
    sub, Bs = [], [Bh[0]]
    for k, t in enumerate(dtau):
        f = np.arange(1, 8) / 7.0
        sub += [t / 7.0] * 7
        Bs += list(Bh[k] + f * (Bh[k + 1] - Bh[k]))
    U2, D2 = F.dos_flujos(np.array(sub), [1.0], np.array(Bs), 395.0)
    assert np.allclose(U, U2[::7], rtol=1e-12)
    assert np.allclose(D, D2[::7], rtol=1e-12)


def test_igual_que_el_infrarrojo_gris_de_m3n():
    """Con un solo termino, son las mismas recurrencias que fase30_multicapa (infrarrojo_bajada/subida), escritas aqui
    con el espesor optico ya multiplicado por la difusividad (en la 2.3 venia calibrado asi)."""
    dtau, Bh = _perfil()
    dt = F.DIFUSIVIDAD * dtau
    tr, g1 = np.exp(-dt), (1 - np.exp(-dt) * (1 + dt)) / dt
    D = np.zeros(len(Bh))
    for k in range(len(dtau)):
        D[k + 1] = D[k] * tr[k] + Bh[k + 1] * (1 - tr[k]) - (Bh[k + 1] - Bh[k]) * g1[k]
    U = np.empty(len(Bh)); U[-1] = 390.0
    for k in range(len(dtau) - 1, -1, -1):
        U[k] = U[k + 1] * tr[k] + Bh[k] * (1 - tr[k]) + (Bh[k + 1] - Bh[k]) * g1[k]
    U2, D2 = F.dos_flujos(dtau, [1.0], Bh, 390.0)
    assert np.allclose(U, U2, rtol=1e-13) and np.allclose(D, D2, rtol=1e-13)


def test_terminos_pesados_y_calentamiento_que_cierra():
    dtau, Bh = _perfil()
    t = np.vstack([dtau, 0.1 * dtau, 10 * dtau])
    w = np.array([0.5, 0.3, 0.2])
    U, D = F.dos_flujos(t, w, Bh, 390.0)
    partes = [F.dos_flujos(t[i], [1.0], Bh, 390.0) for i in range(3)]
    assert np.allclose(U, sum(wi * p[0] for wi, p in zip(w, partes)), rtol=1e-13)
    h = F.calentamiento(U, D)
    assert h.sum() == pytest.approx((U[-1] - U[0]) + (D[0] - D[-1]), abs=1e-10)


# ------------------------------------------------------------------------------------------------ solapamiento
def test_rorr_con_un_gas_nulo_o_gris():
    g, w = K.nodos_g()
    r = np.random.default_rng(3)
    a = np.sort(r.lognormal(0, 2, (5, len(w))), axis=1)
    assert np.allclose(F.combinar_rorr([a, np.zeros_like(a)], w), a, rtol=1e-12)
    gris = np.repeat(r.uniform(0.1, 1, (5, 1)), len(w), axis=1)
    assert np.allclose(F.combinar_rorr([a, gris], w), a + gris, rtol=1e-12)


def test_rorr_cerca_del_solapamiento_aleatorio_en_una_capa():
    """En una capa, el solapamiento aleatorio exacto da la transmision PRODUCTO; con el reagrupamiento a 16 puntos queda
    cerca (medido: < 0,01 en estos casos)."""
    g, w = K.nodos_g()
    r = np.random.default_rng(5)
    for _ in range(5):
        a = np.sort(r.lognormal(-1, 2.5, len(w)))
        b = np.sort(r.lognormal(-1, 2.5, len(w)))
        prod = (w * np.exp(-a)).sum() * (w * np.exp(-b)).sum()
        c = F.combinar_rorr([a[None], b[None]], w)[0]
        assert abs((w * np.exp(-c)).sum() - prod) < 0.01


def test_mapa_agrupa_por_profundidad_y_la_media_aritmetica_suma():
    """El mapa ordena de mas transparente (g = 0) a mas opaco, reparte por anchura, y con la media aritmetica la
    suma de los gases en cada grupo es la media del total (es lineal)."""
    r = np.random.default_rng(7)
    n, m = 6, 4000
    k1 = r.lognormal(0, 2, (n, m)); k2 = r.lognormal(-1, 1.5, (n, m))
    dnu = r.uniform(0.5, 1.5, m)
    tau = 0.1 * (k1 + k2)
    acum = np.vstack([np.zeros(m), np.cumsum(tau, axis=0)])
    g, w = K.nodos_g()
    grupo = F.tabla_mapa(acum, dnu, w)
    # mas opaco en los grupos altos: el nivel donde el espesor desde arriba llega a 1 sube (indice menor) con el grupo
    llega = acum >= 1.0
    nivel = np.where(llega.any(axis=0), llega.argmax(axis=0), n + 1)
    maximo = [nivel[grupo == j].max() for j in range(len(w)) if (grupo == j).any()]
    minimo = [nivel[grupo == j].min() for j in range(len(w)) if (grupo == j).any()]
    assert all(mx_sig <= mn for mn, mx_sig in zip(minimo[:-1], maximo[1:]))
    frac = np.array([dnu[grupo == j].sum() for j in range(len(w))]) / dnu.sum()
    assert np.allclose(frac, w, atol=2e-3)
    s = F.media_por_grupo(0.1 * k1, dnu, grupo, len(w)) + F.media_por_grupo(0.1 * k2, dnu, grupo, len(w))
    assert np.allclose(s, F.media_por_grupo(tau, dnu, grupo, len(w)), rtol=1e-12)


def test_media_de_transmision_limites():
    """Con espesores muy pequeños, la media de la transmision es la media aritmetica pesada (sin perder precision); con
    un grupo homogeneo, el propio espesor; nunca fuera del rango del grupo."""
    r = np.random.default_rng(9)
    t = r.lognormal(0, 1, (3, 500))
    p = r.uniform(0.5, 1.5, 500)
    grupo = r.integers(0, 4, 500)
    for esc in (1e-12, 1e-8):
        m = F.media_transmision_por_grupo(esc * t, p, grupo, 4)
        ar = np.array([[np.average(esc * t[c, grupo == j], weights=p[grupo == j]) for j in range(4)] for c in range(3)])
        assert np.allclose(m, ar, rtol=1e-6)
    hom = np.ones((2, 500)) * 0.7
    assert np.allclose(F.media_transmision_por_grupo(hom, p, grupo, 4), 0.7, rtol=1e-12)
    m = F.media_transmision_por_grupo(5 * t, p, grupo, 4)
    for j in range(4):
        assert np.all(m[:, j] <= (5 * t[:, grupo == j]).max(axis=1) + 1e-12)
        assert np.all(m[:, j] >= (5 * t[:, grupo == j]).min(axis=1) - 1e-12)


def test_media_de_transmision_grupo_opaco_finito():
    t = np.full((1, 10), 1e4)
    m = F.media_transmision_por_grupo(t, np.ones(10), np.zeros(10, dtype=int), 1)
    assert np.isfinite(m).all() and m[0, 0] == pytest.approx(-np.log(1e-300) / F.DIFUSIVIDAD)

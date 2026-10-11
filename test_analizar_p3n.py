# test_analizar_p3n.py -- 4.0.0 (Fase 6.4): las piezas de analizar_p3n.py con casos construidos a mano.
import math

import numpy as np

import analizar_p3n as A


def test_un_cero_exacto_no_es_un_tramo():
    lats = np.array([20.0, 10.0, 0.0, -10.0, -20.0])
    t = A.tramos_de_signo(np.array([2.0, 1.0, 0.0, -1.0, -2.0]), lats)
    assert [x[0] for x in t] == [1, -1] and t[0][4] == 0.0 and t[1][3] == 0.0


def test_tramos_de_signo_con_bordes_interpolados():
    lats = np.array([30.0, 20.0, 10.0, 0.0, -10.0])
    perfil = np.array([1.0, 3.0, -1.0, -2.0, 2.0])
    t = A.tramos_de_signo(perfil, lats)
    assert [x[0] for x in t] == [1, -1, 1]
    assert t[0][3] == 30.0 and math.isclose(t[0][4], 12.5)       # 3 -> -1 entre 20 y 10: cruce en 12,5
    assert math.isclose(t[1][3], 12.5) and math.isclose(t[1][4], -5.0)
    assert t[2][4] == -10.0


def test_celulas_encuentra_todas_sin_ventanas():
    # seis celulas alternas de norte a sur, como una Tierra (Ferrel, Hadley... ) o lo que salga
    lat_c = np.linspace(85, -85, 35)
    perfil = np.sin(np.radians(lat_c) * 6)                        # 6 tramos de signo
    nivel = np.linspace(0, 1, 21)
    psi = np.outer(np.sin(np.pi * nivel), perfil)                  # maximo a media altura, cero en tope y suelo
    k5 = int(np.argmin(np.abs(nivel - 0.5)))
    cs = A.celulas(psi, lat_c, nivel, k5)
    assert len(cs) == 6
    assert [c["signo"] for c in cs] == [1, -1, 1, -1, 1, -1]
    assert all(not c["debil"] for c in cs)
    assert all(c["norte"] >= c["lat"] >= c["sur"] for c in cs)


def test_celula_debil_se_marca_y_no_se_oculta():
    lat_c = np.linspace(85, -85, 35)
    perfil = np.where(lat_c > 0, 10.0, -0.1)
    nivel = np.linspace(0, 1, 21)
    psi = np.outer(np.sin(np.pi * nivel), perfil)
    cs = A.celulas(psi, lat_c, nivel, 10)
    assert len(cs) == 2 and cs[1]["debil"] and not cs[0]["debil"]


def test_nucleos_de_chorro_varios_y_del_este():
    lat = np.linspace(87.5, -87.5, 36)
    nivel = np.linspace(0.025, 0.975, 20)
    forma_v = np.exp(-((nivel - 0.25) / 0.08) ** 2)[:, None]       # maximo en sigma 0,25
    lat_u = 30 * np.exp(-((lat - 52.5) / 6) ** 2) + 20 * np.exp(-((lat - 27.5) / 5) ** 2) \
        - 15 * np.exp(-((lat - 2.5) / 6) ** 2) + 25 * np.exp(-((lat + 47.5) / 6) ** 2)
    u = forma_v * lat_u[None, :]
    nc = A.nucleos_chorro(u, lat, nivel)
    lats = [round(x[1], 1) for x in nc]
    assert any(abs(la - 52.5) <= 2.5 for la in lats) and any(abs(la - 27.5) <= 2.5 for la in lats)
    assert any(abs(la + 47.5) <= 2.5 for la in lats)
    assert len(nc) == 4                                            # tres del oeste y uno del este, sin repetir
    este = [x for x in nc if x[0] < 0]
    assert len(este) == 1 and abs(este[0][1]) <= 2.5
    assert all(0.1 <= x[2] <= 0.5 for x in nc)


def test_bandas_de_viento_y_extremos_locales():
    lat = np.linspace(87.5, -87.5, 36)
    us = -5 * np.cos(np.radians(lat) * 3)                         # este en el ecuador, oeste hacia los polos
    b = A.bandas_viento(us, lat)
    assert [x["signo"] for x in b] == [1, -1, 1]                  # oeste > 30 grados, este entre 30 S y 30 N
    assert math.isclose(b[1]["norte"], 30.0, abs_tol=0.5) and math.isclose(b[1]["sur"], -30.0, abs_tol=0.5)
    z = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    la = np.array([20.0, 10.0, 0.0, -10.0, -20.0])
    assert A.extremos_locales(z, la) == [(3.0, 10.0), (5.0, -10.0)]
    assert A.extremos_locales(z, la, maximos=False) == [(1.0, 20.0), (2.0, 0.0), (4.0, -20.0)]   # en los bordes, con su unica vecina


def test_calendario_con_la_orbita_de_p3n():
    import parametros as P
    from cache_simulacion import precalcular_orbita_cacheada
    from temperatura import PASO_TIEMPO
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    pasos_dia = round(P.ROTACION_PERIODO / PASO_TIEMPO)
    n = len(orbita) // pasos_dia
    c = A.calendario(orbita, pasos_dia, n)
    # periastro y apoastro separados medio año; luz en el tope en la proporcion ((1+e)/(1-e))^2
    sep = (c["d_apo"] - c["d_peri"]) % (len(orbita) / pasos_dia)
    assert abs(sep - len(orbita) / pasos_dia / 2) < 2
    e = P.ORBITA_EXCENTRICIDAD
    assert abs(c["toa_peri"] / c["toa_apo"] - ((1 + e) / (1 - e)) ** 2) < 0.01
    assert c["temporada"][c["d_peri"]] == "periastro" and c["temporada"][c["d_apo"]] == "apoastro"
    frac = (c["temporada"] == "periastro").mean()
    assert abs(frac - 0.25) < 0.01
    assert set(np.unique(c["doceavo"])) == set(range(12))


def test_factor_held_hou_de_p3n():
    # §3 de DISENO_FASE_6.4.md: x1,08 con g, Omega y a de P3N
    assert abs(A.factor_held_hou() - 1.08) < 0.005

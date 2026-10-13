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


# ---------------------------------------------------------------- 3.22.0
def _z_sintetico(dias=8, capas=20, filas=36, seed=1):
    rng = np.random.default_rng(seed)
    return {"u": rng.standard_normal((dias, capas, filas)), "vps": rng.standard_normal((dias, capas, filas - 1)),
            "T": 250 + rng.standard_normal((dias, capas, filas)), "ps": 1e5 + rng.standard_normal((dias, filas)),
            "subpasos_por_dia": np.array([4.0] * (dias - 1) + [2.0]),          # el ultimo dia, incompleto
            "sigma_semi": np.linspace(0, 1, capas + 1), "lat_caras": np.linspace(85, -85, filas - 1),
            "cos_caras": np.cos(np.radians(np.linspace(85, -85, filas - 1))), "anos": 30}


def test_vientos_de_climatologia_es_la_media_de_todos_los_subpasos():
    z = _z_sintetico()
    temporada = np.array(["periastro", "periastro", "", "", "apoastro", "apoastro", "", "periastro"])
    v = A.vientos_de_climatologia(z, {"temporada": temporada})
    n = z["subpasos_por_dia"]
    # la media de todos los subpasos = suma de (media del dia x subpasos del dia) / subpasos
    esperado = (z["u"] * n[:, None, None]).sum(axis=0) / n.sum()
    assert np.allclose(v["anual"]["u"], esperado)
    sel = temporada == "periastro"
    assert np.allclose(v["periastro"]["vps"], (z["vps"][sel] * n[sel, None, None]).sum(axis=0) / n[sel].sum())
    assert v["anos"] == 30 and v["sigma"].shape == (21,) and "apoastro" in v
    # con la misma forma que los del año extra: la funcion de corriente sale igual
    from validar_i16 import funcion_corriente
    psi = funcion_corriente(v["anual"]["vps"], v["sigma"], v["cos_caras"], 4.6767e6, 9.0)
    assert psi.shape == (21, 35)


def test_seccion_vientos_dice_de_donde_salen():
    z = _z_sintetico()
    v = A.vientos_de_climatologia(z, {"temporada": np.array([""] * 8)})
    L, campos = [], {}
    A.seccion_vientos(v, L, campos)
    texto = "\n".join(L)
    assert "media de los 30 años de la climatologia" in texto and "psi_anual" in campos
    L = []
    A.seccion_vientos(None, L, {})
    assert "sin datos" in L[0]


def test_error_diferencia_con_la_del_resumen():
    from fase6_clima import error_media
    rng = np.random.default_rng(5)
    b1, b2 = 20 + rng.standard_normal((30, 36)), 20 + rng.standard_normal((12, 36))
    err, err_g = A.error_diferencia({"climatologia": {"bandas_aire2m_C": b1}}, {"climatologia": {"bandas_aire2m_C": b2}})
    assert np.allclose(err, np.sqrt(error_media(b1) ** 2 + error_media(b2) ** 2))
    assert err.shape == (36,) and 0 < err_g < err.max()        # la media global varia menos que una fila


def test_seccion_d_con_hielo_amplitud_y_error(tmp_path):
    import pickle
    rng = np.random.default_rng(2)
    dias = 24
    base = 10 + np.zeros((dias, 36, 72))
    hielo = np.zeros((dias, 36, 72)); hielo[:, 0] = 1.0
    rc = {"reg_media": base, "agua": {"precipitacion": np.full((36, 72), 3.0)}, "hielo_espesor": hielo,
          "mapa": {"huella": "x"}, "climatologia": {"bandas_aire2m_C": 10 + 0.1 * rng.standard_normal((30, 36))}}
    otra = dict(rc)
    ciclo = np.sin(2 * np.pi * np.arange(dias) / dias)[:, None, None]
    otra["reg_media"] = base + 1.0 + 2.0 * ciclo                     # 1 K mas y 4 K de amplitud estacional
    otra["hielo_espesor"] = np.zeros_like(hielo)
    otra["climatologia"] = {"bandas_aire2m_C": 11 + 0.1 * rng.standard_normal((12, 36)), "anos_promediados": 12}
    otra["sensibilidad"] = {"oceano_difusion": 2.0, "oceano_capa_mezcla": 1.0}
    with open(tmp_path / "resultado.pkl", "wb") as fh:
        pickle.dump(otra, fh)
    cal = {"doceavo": np.repeat(np.arange(12), 2)}
    from fase1_geografia import TIERRA
    tipo = np.zeros((36, 72), dtype=int); tipo[10:20, 10:20] = TIERRA
    L = []
    A.seccion_d(L, rc, [str(tmp_path)], tipo, cal)
    texto = "\n".join(L)
    assert "global: dT +1.00" in texto and ")*" in texto            # 1 K con error de ~0,03 K: distinguible
    assert "N max " in texto and "->  0.00" in texto                 # el hielo del norte desaparece
    assert "amplitud estacional" in texto and "cambio medio +3." in texto     # doceavos de una sinusoide de 4 K
    # con el mismo oceano lo dice
    otra["sensibilidad"] = {"oceano_difusion": 1.0, "oceano_capa_mezcla": 1.0}
    with open(tmp_path / "resultado.pkl", "wb") as fh:
        pickle.dump(otra, fh)
    L = []
    A.seccion_d(L, rc, [str(tmp_path)], tipo, cal)
    assert "el mismo oceano" in "\n".join(L)


def test_seccion_hielo_de_la_climatologia():
    hielo = np.zeros((10, 36, 72)); hielo[:5, 0] = 1.0; hielo[:, -1] = 0.5
    L = []
    h = A.seccion_hielo({"hielo_espesor": hielo, "climatologia": {"anos_promediados": 30}}, L)
    assert h["N_min"] == 0 and h["N_max"] > 0 and h["S_min"] == h["S_max"] > 0 and "30 años" in L[0]


def test_tierra_mar_separa_posicion_y_misma_latitud():
    """Un planeta donde la temperatura solo depende de la latitud y la tierra esta toda en el ecuador: la tierra sale
    mas calida por POSICION, y a la misma latitud no hay diferencia. Si ademas la tierra esta 2 K mas caliente que el
    mar de su fila, aparece en "misma latitud"."""
    lat = A.LAT
    T = np.repeat((30 - 0.5 * np.abs(lat))[:, None], 72, axis=1)
    tierra = np.zeros((36, 72), dtype=bool); tierra[16:20, :20] = True
    r = A.tierra_mar(T, tierra)
    assert r["total"] > 5 and abs(r["misma_latitud"]) < 1e-12 and abs(r["posicion"] - r["total"]) < 1e-12
    r = A.tierra_mar(T + 2.0 * tierra, tierra)
    assert abs(r["misma_latitud"] - 2.0) < 1e-12 and len(r["filas"]) == 4 and abs(r["filas"][0]["pct"] - 100 * 20 / 72) < 1e-9
    assert A.tierra_mar(T, np.zeros((36, 72), dtype=bool)) is None


# ---------------------------------------------------------------- 3.22.1
def test_bloqueantes_de_una_climatologia_continuada():
    """Una climatologia continuada (--continuar-desde) no tiene fase de equilibrio propia: A3 y A5 se leen como
    heredados de la de partida, en lugar de dar 'resumen incompleto' (fallo de la 3.22.0, visto el 13/10/2026)."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "herramientas"))
    import criterio_cierre_6_3 as CC
    filas = "\n".join(f"  {i}: cierre energia 1.0e-08 | error dinamica {-0.047 + 0.001 * ((-1) ** i):.4f} W/m2 | resto 0.0 "
                      f"| agua atm 1.0e-13 tierra 1.0e-14" for i in range(1, 31))
    res = ("RESUMEN\n\nEQUILIBRIO: el de /x/p3n_3.19.0 (15 años de equilibrio y 30 de climatologia antes de esta; año "
           "simulado de partida 45). Esta climatologia continua desde su estado final.\n\n"
           "CLIMATOLOGIA: 30 años registrados de 30 | normal\n\nPROCEDENCIA\n  tramo 1: M3N 3.21.0, codigo abc123def456\n\n"
           f"POR AÑO\n{filas}\n\nDERIVA: sin deriva significativa\n  aire 2 m global: -0.0002 +- 0.0011 /año (p 0.869) K | x\n")
    t = CC.Tabla()
    CC.seccion_a(t, res)
    d = {c: (e, txt) for c, _, e, txt in t.filas}
    assert d["A3"][0] == "✅" and "heredado" in d["A3"][1] and "p3n_3.19.0" in d["A3"][1]
    assert d["A5"][0] == "✅" and "de la de partida" in d["A5"][1]
    assert d["A1"][0] == "✅" and d["A2"][0] == "✅" and d["A4"][0] == "✅" and d["A6"][0] == "✅"

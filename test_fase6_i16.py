# test_fase6_i16.py -- v3.1-pre11: pruebas del modelo ACOPLADO (interruptor I16, "nucleo_dinamico").
#
# Uso:  python -m pytest test_fase6_i16.py      (alrededor de un minuto)
#
# Son simulaciones MUY cortas (unos pocos pasos de la fisica por "año"): comprueban propiedades que deben
# cumplirse en CADA paso, no el clima.
#   - el agua de la atmosfera y la del suelo se conservan EXACTAMENTE (redondeo);
#   - el balance de energia cierra (los residuos declarados son mucho menores que la cota);
#   - una simulacion interrumpida y retomada desde el punto de control da EXACTAMENTE lo mismo;
#   - con I16 apagado no se crea nada del nucleo (lo de siempre).
# El clima no se puede juzgar aqui: con una orbita recortada a unos pocos pasos, el estado inicial (el
# balance de Newton con la luz MEDIA de esa orbita) tiene la superficie a mas de 400 K bajo la estrella.
# La estabilidad durante semanas se comprueba aparte, desde un estado de equilibrio (DISENO_FASE6_3.md §6.9).

import numpy as np
import pytest

import fase2b_atmosfera as F

PASOS_ANO = 4          # pasos de la fisica por "año" de prueba (muy por debajo de los ~25 del arranque frio)


def _argumentos(acua=False):
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA
    from cache_simulacion import precalcular_orbita_cacheada
    import parametros as P
    from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS
    orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)[:PASOS_ANO]
    lon = np.ones((FILAS, 1)) * LONGITUDES_GRADOS[None]
    lat = LATITUDES_GRADOS[:, None] * np.ones((1, COLUMNAS))
    tipo = np.where((np.abs(lon) < 40) & (np.abs(lat) < 70), TIERRA, AGUA)
    if acua:
        tipo = np.full_like(tipo, AGUA)
    alt = np.where(tipo == TIERRA, 600.0, 0.0)
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral",
              "nucleo_dinamico"):
        I[k] = True
    arg = (orb, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55)
    kw = dict(interruptores=I, acelerar=False, tolerancia_convergencia=1e-9)
    return arg, kw


@pytest.fixture(scope="module")
def corrida():
    arg, kw = _argumentos()
    return F.simular_fase2b(*arg, max_anos=2, **kw)


def test_i16_agua_exacta(corrida):
    a = corrida["agua"]
    assert abs(a["cierre_agua_atmosfera_kg_m2"]) < 1e-12        # kg/m2 (la columna tiene ~10-60 kg/m2)
    assert abs(a["cierre_agua_tierra_kg_m2"]) < 1e-12
    assert a["recortes_kg_m2"] == 0.0                             # el transporte no deja vapor negativo


def test_i16_energia_cierra_salvo_la_dinamica(corrida):
    """Quitando lo que la dinamica no conserva (y el calor de rozamiento que devuelve la capa limite), el
    balance de energia cierra como en la v3.1 sin dinamica (redondeo). El error de la dinamica se mide y se
    declara (DISENO_FASE6_3.md §6.9); en estas simulaciones de pocos pasos desde el reposo es grande porque
    el aire se esta poniendo en movimiento a ~1500 W/m2."""
    e = corrida["energia"]
    assert abs(e["cierre_sin_dinamica"]) < 1e-9
    assert np.isfinite(e["error_dinamica_W_m2"]) and e["calor_rozamiento_W_m2"] > 0


def test_i16_punto_de_control_da_lo_mismo(tmp_path):
    arg, kw = _argumentos()
    seguida = F.simular_fase2b(*arg, max_anos=3, **kw)
    archivo = str(tmp_path / "estado.pkl")
    F.simular_fase2b(*arg, max_anos=2, archivo_estado=archivo, **kw)
    retomada = F.simular_fase2b(*arg, max_anos=3, archivo_estado=archivo, **kw)
    for k in ("reg_media", "suelo_media", "hielo_espesor"):
        assert np.array_equal(seguida[k], retomada[k]), k
    assert np.array_equal(seguida["agua"]["precipitacion"], retomada["agua"]["precipitacion"])
    assert np.array_equal(seguida["flujos"]["T_atm"], retomada["flujos"]["T_atm"])
    assert seguida["energia"]["cierre_relativo"] == retomada["energia"]["cierre_relativo"]


def test_i16_nucleo_en_reposo_isotermo_sobre_montanas():
    """T4 (DISENO_FASE6_3.md §6.7.7): el nucleo con la configuracion de I16 (niveles de M3N, avanzar con
    hiperdifusion, su energia cinetica como calor y la correccion a superficies de presion) mantiene en
    reposo una atmosfera isoterma sobre las montañas del caso 5 de Williamson (redondeo)."""
    import fase30_multicapa as M30
    from fase6_nucleo import NucleoSeco
    from fase6_aguas_someras import caso5, A_TIERRA, OMEGA_TIERRA, G_TIERRA
    hs = caso5()[3]
    m = NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, M30.sigma_seminiveles(20), phis=G_TIERRA * hs)
    m.preparar_semiimplicito(450.0)
    m.preparar_hiperdifusion(0.5, calor_rozamiento=True, correccion_presion=True)
    T0 = 280.0
    ps = 1e5 * np.exp(-G_TIERRA * hs / (287.0 * T0))
    ant, act = m.arrancar((ps, np.full((m.N, 36, 72), T0), np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72))))
    for _ in range(100):
        ant, act = m.avanzar(ant, act)
    assert np.abs(act[2]).max() < 1e-9 and np.abs(act[3]).max() < 1e-9
    assert np.abs(act[1] - T0).max() < 1e-9 and np.abs(act[0] - ps).max() < 1e-6


# ---------------- v3.1-pre12: equilibrio y arranque caliente ----------------

def test_equilibrio_pendiente_y_criterio():
    import fase6_equilibrio as EQ
    assert EQ.pendiente([1.0, 2.0, 3.0, 4.0, 5.0]) == pytest.approx(1.0, abs=1e-15)
    assert EQ.pendiente([3.0, 3.0, 3.0]) == 0.0
    base = {"N": [0.05] * 5, "T2m": [288.0, 288.01, 287.99, 288.0, 288.01], "hielo": [0.1] * 5,
            "agua_suelo": [0.5] * 5}
    ok, v = EQ.evaluar(base)
    assert ok and set(v) == {"N", "T2m", "hielo", "agua_suelo"}
    assert not EQ.evaluar({k: x[:4] for k, x in base.items()})[0]            # menos de 5 años
    for k, malo in (("N", [0.3] * 5), ("T2m", [288.0, 288.03, 288.06, 288.09, 288.12]),
                    ("hielo", [0.1, 0.102, 0.104, 0.106, 0.108]), ("agua_suelo", [0.5, 0.52, 0.54, 0.56, 0.58])):
        ok, v = EQ.evaluar(dict(base, **{k: malo}))
        assert not ok and not v[k][2], k
    ok, v = EQ.evaluar(dict(base, hielo=[None] * 5))                        # sin oceano: no cuenta
    assert ok and "hielo" not in v
    # solo cuenta la ventana: un comienzo con deriva no impide el equilibrio
    ok, _ = EQ.evaluar({k: ([x[0] + 5.0] * 3 if k != "N" else [2.0] * 3) + x for k, x in base.items()})
    assert ok


def test_arranque_caliente_desde_el_modelo_de_2_capas(tmp_path):
    """El estado final (oceano, suelo, hielo) de una simulacion de la v2.4.3 sirve de estado inicial con I16;
    el aire parte en reposo. Con el mismo estado y un punto de control, retomar da lo mismo bit a bit."""
    arg, kw = _argumentos()
    dos_capas = F.simular_fase2b(*arg, max_anos=2, acelerar=False, tolerancia_convergencia=1e-9)
    ini = dos_capas["estado_final"]
    assert ini["T_col"].ndim == 3 and set(ini["HIELO"]) == {"h", "Ts"}
    seguida = F.simular_fase2b(*arg, max_anos=3, estado_inicial=ini, **kw)
    archivo = str(tmp_path / "estado.pkl")
    F.simular_fase2b(*arg, max_anos=2, estado_inicial=ini, archivo_estado=archivo, **kw)
    retomada = F.simular_fase2b(*arg, max_anos=3, estado_inicial=ini, archivo_estado=archivo, **kw)
    assert np.array_equal(seguida["flujos"]["T_atm"], retomada["flujos"]["T_atm"])
    assert seguida["equilibrio"]["serie"] == retomada["equilibrio"]["serie"]
    assert len(seguida["equilibrio"]["serie"]["N"]) == 3 and not seguida["convergido"]
    assert abs(seguida["agua"]["cierre_agua_atmosfera_kg_m2"]) < 1e-12
    assert abs(seguida["energia"]["cierre_sin_dinamica"]) < 1e-9
    # otro estado inicial = otra simulacion: el punto de control no se mezcla
    otro = {"T_col": ini["T_col"] + 1.0, "HIELO": ini["HIELO"]}
    with pytest.raises(ValueError):
        F.simular_fase2b(*arg, max_anos=3, estado_inicial=otro, archivo_estado=archivo, **kw)
    with pytest.raises(ValueError):
        F.simular_fase2b(*arg, max_anos=1, estado_inicial={"T_col": ini["T_col"][..., :2], "HIELO": ini["HIELO"]}, **kw)


def test_climatologia_error_de_la_media():
    import fase6_clima as C
    rng = np.random.default_rng(0)
    y = rng.standard_normal((2000, 3))
    assert np.allclose(C.error_media(y), 1 / np.sqrt(2000), rtol=0.05)          # ruido blanco: sigma/sqrt(N)
    ar = np.zeros((4000, 1))
    for t in range(1, 4000):
        ar[t] = 0.6 * ar[t - 1] + rng.standard_normal()
    n_eff = 4000 * 0.4 / 1.6
    assert C.error_media(ar)[0] == pytest.approx(ar.std(ddof=1) / np.sqrt(n_eff), rel=0.05)
    cumple, info = C.criterio_climatologia(np.full((5, 36), 15.0), np.full((5, 36), 2.0))
    assert cumple and info["bandas_T_sin_cumplir"] == 0
    assert not C.criterio_climatologia(np.full((4, 36), 15.0), np.full((4, 36), 2.0))[0]      # minimo 5 años
    T = np.full((6, 36), 15.0); T[::2, 3] += 1.0                                              # una banda ruidosa
    cumple, info = C.criterio_climatologia(T, np.full((6, 36), 2.0))
    assert not cumple and info["bandas_T_sin_cumplir"] == 1
    P = np.full((6, 36), 0.01); P[::2] += 0.02                                                # casi seco
    assert C.criterio_climatologia(np.full((6, 36), 15.0), P)[0]                              # < 0,05 mm/dia
    # desviacion entre años sin perdida de precision (valores grandes, diferencias pequeñas)
    xs = [300.0 + 0.01 * rng.standard_normal((3, 4)) for _ in range(7)]
    A = C.Acumulador()
    for x in xs:
        A._cuadrado("reg_media", x)
    var = np.maximum(A.suma2["reg_media"] - A.suma2["reg_media_d"] ** 2 / 7, 0.0) / 6
    assert np.allclose(np.sqrt(var), np.std(xs, axis=0, ddof=1), rtol=1e-9, atol=0)


def test_climatologia_cortada_y_retomada_da_lo_mismo(tmp_path, monkeypatch):
    """La cadena completa (2 capas -> I16 -> años registrados) con un corte justo despues de que el simulador
    guarde un año registrado y ANTES de que se guarde la climatologia: al retomar, el resultado es identico."""
    import fase6_clima as C
    import parametros as P
    arg, kw = _argumentos()
    from cache_simulacion import precalcular_orbita_cacheada
    orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    arg = (orb[:round(P.ROTACION_PERIODO / 992)],) + arg[1:]              # "años" de 1 dia: hace falta un dia entero
    I = kw["interruptores"]
    # solo se prueba la cadena: arranque desde una superficie templada y uniforme (sin hielo), no el equilibrio
    from rejilla import FILAS, COLUMNAS
    n_suelo = F.simular_fase2b(*arg, max_anos=0)["estado_final"]["T_col"].shape[-1]
    ini = {"T_col": np.full((FILAS, COLUMNAS, n_suelo), 285.0),
           "HIELO": {"h": np.zeros((FILAS, COLUMNAS)), "Ts": np.full((FILAS, COLUMNAS), 271.35)}}
    DOS = dict(estado_inicial=ini)
    seguida = C.simular_clima(arg, I, str(tmp_path / "a"), min_anos=2, max_anos=2, max_anos_equilibrio=1,
                              informar=lambda *x: None, **DOS)
    assert seguida["climatologia"]["anos_promediados"] == 2
    original = F.simular_fase2b
    llamadas = {"n": 0}

    def con_corte(*a, **k):
        r = original(*a, **k)
        llamadas["n"] += 1
        if llamadas["n"] == 2:
            raise KeyboardInterrupt
        return r
    monkeypatch.setattr(F, "simular_fase2b", con_corte)
    with pytest.raises(KeyboardInterrupt):
        C.simular_clima(arg, I, str(tmp_path / "b"), min_anos=2, max_anos=2, max_anos_equilibrio=1,
                        informar=lambda *x: None, **DOS)
    monkeypatch.setattr(F, "simular_fase2b", original)
    retomada = C.simular_clima(arg, I, str(tmp_path / "b"), min_anos=2, max_anos=2, max_anos_equilibrio=1,
                               informar=lambda *x: None, **DOS)
    for k in ("reg_media", "reg_min", "horario_aire2m"):
        assert np.array_equal(seguida[k], retomada[k]), k
    assert np.array_equal(seguida["agua"]["precipitacion"], retomada["agua"]["precipitacion"])
    assert np.array_equal(seguida["extremos"]["reg_max"], retomada["extremos"]["reg_max"])
    assert np.array_equal(seguida["climatologia"]["bandas_aire2m_C"], retomada["climatologia"]["bandas_aire2m_C"])
    assert seguida["climatologia"]["equilibrio"]["anos"] == retomada["climatologia"]["equilibrio"]["anos"] == 1


def test_correccion_presion_limitada_en_escalones_de_relieve():
    """v3.1-pre13 (§6.11): junto a un escalon de ~2800 m (Antartida a 5 grados) la correccion de la difusion
    de T a superficies de presion equivaldria a desplazar el perfil ~3 capas (fuera de su validez, e inestable);
    se limita a media capa. Con relieve suave (montaña del caso 5 de Williamson) no cambia nada."""
    import fase30_multicapa as M30
    import fase6_nucleo as N6
    from fase6_aguas_someras import caso5, A_TIERRA, OMEGA_TIERRA, G_TIERRA
    sh = M30.sigma_seminiveles(20)
    sig_m = 0.5 * (np.asarray(sh[1:]) + np.asarray(sh[:-1]))
    T = 290.0 - 40.0 * (1 - sig_m)[:, None, None] * np.ones((1, 36, 72))     # T que cambia con la altura
    for alto, cambia in ((np.where((np.arange(36)[:, None] >= 31) & (np.abs(np.arange(72)[None] - 36) < 6), 2800.0, 0.0), True),
                         (caso5()[3], False)):
        m = N6.NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, sh, phis=G_TIERRA * alto)
        m.preparar_semiimplicito(450.0)
        m.preparar_hiperdifusion(0.5, calor_rozamiento=True, correccion_presion=True)
        ps = 1e5 * np.exp(-G_TIERRA * alto / (287.0 * 260.0))
        corr = m._correccion_presion(T, ps, 900.0)
        lim = N6.LIMITE_CORRECCION_PRESION
        N6.LIMITE_CORRECCION_PRESION = np.inf
        try:
            libre = m._correccion_presion(T, ps, 900.0)
        finally:
            N6.LIMITE_CORRECCION_PRESION = lim
        assert np.isfinite(corr).all()
        if cambia:
            assert np.abs(corr).max() < np.abs(libre).max()
            dT = np.zeros_like(T)
            dT[:-1] += np.asarray(sh)[1:-1, None, None] * (T[1:] - T[:-1])
            dT[1:] += np.asarray(sh)[1:-1, None, None] * (T[1:] - T[:-1])
            delps = corr / (0.5 / (m._d3 * ps[None]) * np.where(dT == 0, 1.0, dT))       # delps efectivo
            c = (np.abs(delps) * sig_m[:, None, None] / (m._d3 * ps[None])).max(axis=0)
            assert c.max() <= lim * (1 + 1e-12)
        else:
            assert np.array_equal(corr, libre)

# test_fase6_i16.py -- v3.9.0: pruebas del modelo ACOPLADO (interruptor I16, "nucleo_dinamico").
#
# Uso:  python -m pytest test_fase6_i16.py      (alrededor de un minuto)
#
# Son simulaciones MUY cortas (unos pocos pasos de la fisica por "año"): comprueban propiedades que deben
# cumplirse en CADA paso, no el clima.
#   - el agua de la atmosfera y la del suelo se conservan EXACTAMENTE (redondeo);
#   - el balance de energia cierra (los residuos declarados son mucho menores que la cota) y, desde la v3.11.0,
#     el corrector global de energia devuelve lo que la dinamica no conserva;
#   - una simulacion interrumpida y retomada desde el punto de control da EXACTAMENTE lo mismo;
#   - con I16 apagado no se crea nada del nucleo (lo de siempre).
# El clima no se puede juzgar aqui: con una orbita recortada a unos pocos pasos, el estado inicial (el
# balance de Newton con la luz MEDIA de esa orbita) tiene la superficie a mas de 400 K bajo la estrella.
# La estabilidad durante semanas se comprueba aparte, desde un estado de equilibrio (DISENO_FASE_6.3.md §6.9).

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
    declara (DISENO_FASE_6.3.md §6.9); en estas simulaciones de pocos pasos desde el reposo es grande porque
    el aire se esta poniendo en movimiento a ~1500 W/m2."""
    e = corrida["energia"]
    assert abs(e["cierre_sin_dinamica"]) < 1e-9
    assert np.isfinite(e["error_dinamica_W_m2"]) and e["calor_rozamiento_W_m2"] > 0


def test_i16_corrector_de_energia(corrida, monkeypatch):
    """v3.11.0 (DISENO_FASE_6.3.md §6.13): el corrector global devuelve a la atmosfera lo que la dinamica no
    conserva, con un incremento uniforme de T. Con el, el balance TOTAL cierra (salvo el desfase de un paso del
    calor de rozamiento, que en el equilibrio se compensa); sin el, el error de la dinamica queda en el balance."""
    e = corrida["energia"]
    assert e["correccion_energia_W_m2"] != 0.0
    assert abs(e["correccion_energia_W_m2"] + e["error_dinamica_W_m2"]) < 0.05 * abs(e["error_dinamica_W_m2"])
    monkeypatch.setattr(F, "CORRECTOR_ENERGIA_I16", False)
    arg, kw = _argumentos()
    sin = F.simular_fase2b(*arg, max_anos=2, **kw)["energia"]
    assert sin["correccion_energia_W_m2"] == 0.0 and abs(sin["cierre_sin_dinamica"]) < 1e-9
    assert abs(e["cierre_relativo"]) < 0.05 * abs(sin["cierre_relativo"])


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
    """T4 (DISENO_FASE_6.3.md §6.7.7): el nucleo con la configuracion de I16 (niveles de M3N, avanzar con
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


# ---------------- v3.10.0: equilibrio y arranque caliente ----------------

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


def test_correccion_presion_interpola_sin_extrapolar():
    """v3.10.2 (§6.12): la correccion de la difusion de T a superficies de presion es una interpolacion lineal
    entre cada capa y su vecina en la direccion del desplazamiento, con peso <= 0,5:
      - con un perfil lineal en sigma y desplazamientos pequeños, coincide con la forma de CAM (Taylor centrado)
        en las capas interiores;
      - nunca extrapola con el perfil del modelo: en la capa de abajo, con una inversion fortisima (meseta
        antartica en la noche polar), no la extrapola (gradiente acotado entre isotermo y adiabatico seco);
      - en ningun caso crea extremos nuevos."""
    import fase30_multicapa as M30
    import fase6_nucleo as N6
    from fase6_aguas_someras import caso5, A_TIERRA, OMEGA_TIERRA, G_TIERRA
    sh = np.asarray(M30.sigma_seminiveles(20))
    sm = 0.5 * (sh[1:] + sh[:-1])
    lim = N6.LIMITE_CORRECCION_PRESION

    def nucleo(alto):
        m = N6.NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, sh, phis=G_TIERRA * alto)
        m.preparar_semiimplicito(450.0)
        m.preparar_hiperdifusion(0.5, calor_rozamiento=True, correccion_presion=True)
        return m, 1e5 * np.exp(-G_TIERRA * alto / (287.0 * 260.0))

    # 1) relieve suave (montaña del caso 5 a la decima parte) y T lineal en sigma, T = a + b sigma: la
    #    interpolacion es EXACTA, corr = b dsigma con dsigma = sigma delps/p_s (lo que la forma de CAM aproxima;
    #    con capas desiguales, la de CAM se aparta ~8 %); arriba, hacia el tope, 0; abajo, hacia el suelo, con el
    #    gradiente de las dos capas de abajo en ln(sigma), acotado entre isotermo y adiabatico seco
    m, ps = nucleo(0.1 * caso5()[3])
    b_ = 70.0
    T = (220.0 + b_ * sm)[:, None, None] * np.ones((1, 36, 72))
    corr = m._correccion_presion(T, ps, 900.0)
    Is, _ = m._inversas_difusion(900.0)
    pk = np.ascontiguousarray(np.fft.rfft(ps, axis=-1).T)
    delps = ps - np.fft.irfft((np.matmul(Is, pk[..., None].real)[..., 0]
                               + 1j * np.matmul(Is, pk[..., None].imag)[..., 0]).T, n=72, axis=-1)
    exacta = b_ * sm[:, None, None] * (delps / ps)[None]
    gap = np.diff(sm)
    d_ab = np.minimum(np.maximum(sm[-1] * delps / ps, 0.0), lim * gap[-1])               # abajo, hacia el suelo:
    pend = np.clip((T[-1] - T[-2]) / np.log(sm[-1] / sm[-2]), 0.0, (287.0 / 1004.0) * T[-1])  # gradiente acotado
    abajo = pend * np.log1p(d_ab / sm[-1])
    exacta[-1] = np.where(delps > 0, abajo, exacta[-1])
    exacta[0] = np.where(delps < 0, 0.0, exacta[0])
    assert np.abs(exacta).max() > 1e-3 and np.allclose(corr, exacta, rtol=1e-10, atol=1e-14)

    # 2) meseta de 2800 m y una inversion fortisima en la capa de abajo
    alto = np.where((np.arange(36)[:, None] >= 31) & (np.abs(np.arange(72)[None] - 36) < 6), 2800.0, 0.0)
    m, ps = nucleo(alto)
    T = (220.0 + 40.0 * sm)[:, None, None] * np.ones((1, 36, 72))
    T[-1] = 70.0
    corr = m._correccion_presion(T, ps, 900.0)
    assert np.isfinite(corr).all() and np.abs(corr).max() > 0
    # la capa de abajo solo se acerca a la de arriba, como mucho la mitad del camino
    assert (corr[-1] * (T[-2] - T[-1]) >= 0).all() and (np.abs(corr[-1]) <= lim * np.abs(T[-2] - T[-1]) + 1e-12).all()
    # ningun extremo nuevo: T + corr queda entre los valores de la capa y sus vecinas
    Tn = T + corr
    vec_min = np.minimum(T, np.minimum(np.concatenate([T[:1], T[:-1]]), np.concatenate([T[1:], T[-1:]])))
    vec_max = np.maximum(T, np.maximum(np.concatenate([T[:1], T[:-1]]), np.concatenate([T[1:], T[-1:]])))
    assert (Tn >= vec_min - 1e-9).all() and (Tn <= vec_max + 1e-9).all()


def test_procedencia_registra_y_avisa_si_cambia_el_codigo(tmp_path):
    """3.14.0 (DISENO_FASE_6.3.md §6.16): cada vez que una simulacion larga empieza o continua queda un tramo en
    procedencia.json; si el codigo no es el del tramo anterior, avisa (sin parar); una simulacion empezada antes
    de la 3.14.0 (sin procedencia y con años hechos) tambien avisa."""
    import fase6_clima as C
    avisos = []
    t = C.registrar_procedencia(str(tmp_path), "9.9.9", "aaa", 0, 0, avisos.append)
    assert len(t) == 1 and not avisos and t[0]["version"] == "9.9.9" and t[0]["codigo"] == "aaa"
    t = C.registrar_procedencia(str(tmp_path), "9.9.9", "aaa", 7, 0, avisos.append)
    assert len(t) == 2 and not avisos
    t = C.registrar_procedencia(str(tmp_path), "9.9.10", "bbb", 12, 3, avisos.append)
    assert len(t) == 3 and len(avisos) == 1 and "ha cambiado" in avisos[0]
    assert C.leer_procedencia(str(tmp_path)) == t
    otra = tmp_path / "antigua"
    otra.mkdir()
    C.registrar_procedencia(str(otra), "9.9.9", "aaa", 30, 5, avisos.append)
    assert len(avisos) == 2 and "antes de la 3.14.0" in avisos[1]

# test_modelo.py -- pruebas rapidas del modelo actual (v2.4.2).
#
# Uso:  python -m pytest test_modelo.py      (unos segundos)
#
# Complementan a test_fase0.py (que comprueba la Fase 0 frente al modelo
# de un punto). Aqui se prueba lo que se puede comprobar SIN simular un
# clima completo (que tarda ~20 minutos): el formato de exportacion, el
# calendario, la clave de la cache, el puente con C3N y la conservacion
# de la energia del calor del oceano profundo.

import json
import math

import numpy as np
import pytest

import cache_simulacion as C
import exportar_clima as X
import fase2b_atmosfera as F
from orbita import info_estaciones
from parametros import ORBITA_PERIODO, ROTACION_PERIODO


# ---- formato de exportacion (Fase 4) ----

def test_serie_ida_y_vuelta():
    x = np.array([[-41.237, 0.0, 15.5, 327.67], [np.nan, 32.123, -0.004, -327.67]])
    y = X.leer_serie(X.serie(x))
    assert np.isnan(y[1, 0])
    finito = np.isfinite(x)
    assert np.max(np.abs(y[finito] - x[finito])) <= 0.005 + 1e-9


def test_serie_fuera_de_rango_avisa():
    with pytest.raises(ValueError):
        X.serie(np.array([327.68]))
    with pytest.raises(ValueError):
        X.serie(np.array([-327.68]))   # antes se convertia en "sin dato" sin avisar


def test_huella_mapa_estable_y_sensible():
    tipo = np.zeros((36, 72), dtype=int); tipo[10:20, 10:20] = 1
    alt = np.where(tipo == 1, 500.0, 0.0)
    h1 = X.huella_mapa(tipo, alt)
    assert h1 == X.huella_mapa(tipo.copy(), alt.copy())
    alt2 = alt.copy(); alt2[15, 15] += 1
    assert h1 != X.huella_mapa(tipo, alt2)


# ---- calendario ----

def test_hitos_del_año_en_el_dia_que_los_contiene():
    num_dias = int(ORBITA_PERIODO // ROTACION_PERIODO)
    hitos = {h["nombre"]: h["dia"] for h in X.fechas_clave(num_dias)}
    for est in info_estaciones():
        # el dia real d (el dia 1 empieza en 1.0) esta en el registro floor(d - 1)
        indice = int(math.floor(est["dia_inicio"] - 1)) % num_dias
        assert indice in hitos.values()


# ---- cache ----

def _clave(**cambios):
    tipo = np.zeros((36, 72), dtype=int)
    args = dict(datos_orbita=[(1.0, 0.0, 0.0)], tipo_superficie=tipo, altitud_metros=tipo * 0.0,
                emisividad=0.77, albedo_por_tipo={0: 0.1, 1: 0.3}, inercia_por_tipo={0: 1.0, 1: 2.0},
                profundidad_optica=0.3, D=0.55, paso_tiempo=900, tolerancia_convergencia=0.015,
                max_anos=50, acelerar=True, estado_inicial="libre")
    args.update(cambios)
    return C._hash_de_objeto(C.clave_fase2b(**args))


def test_clave_cache_determinista():
    assert _clave() == _clave()


def test_clave_cache_cambia_con_los_argumentos():
    assert _clave() != _clave(max_anos=10)
    assert _clave() != _clave(acelerar=False)
    assert _clave() != _clave(estado_inicial="helado")


def test_clave_cache_incluye_el_codigo_y_no_depurar(monkeypatch):
    base = _clave()
    monkeypatch.setattr(F, "DEPURAR", not F.DEPURAR)
    assert _clave() == base, "DEPURAR solo cambia lo que se imprime: no debe invalidar la cache"
    monkeypatch.setattr(C, "huella_codigo", lambda *a, **k: "otro-codigo")
    assert _clave() != base, "un cambio en el codigo de la fisica debe invalidar la cache"


# ---- puente con C3N ----

def test_puente_pone_a_cero_la_altitud_del_agua(tmp_path, monkeypatch):
    import puente_c3n
    agua_tierra = [0] * (36 * 72); agua_tierra[100] = 1
    altitud = [11.0] * (36 * 72); altitud[100] = 800.0
    ruta = tmp_path / "puente.json"
    ruta.write_text(json.dumps({"filas": 36, "columnas": 72, "agua_tierra": agua_tierra,
                                "altitud_metros": altitud, "nombre_mapa_origen": "prueba"}))
    monkeypatch.setattr(puente_c3n, "RUTA_PUENTE", str(ruta))
    tipo, alt, nombre = puente_c3n.cargar_mapa_activo_de_c3n()
    assert nombre == "prueba"
    assert alt.ravel()[100] == 800.0
    assert np.all(alt.ravel()[np.arange(36 * 72) != 100] == 0.0)


# ---- conservacion de la energia del oceano profundo ----

def test_compensacion_oceano_profundo_conserva_la_energia():
    """Lo que entra por la base del hielo es exactamente lo que se retira
    del oceano (formula de fase2b_atmosfera.py, v2.4.1), tambien en los
    casos limite: todo el oceano helado, o casi."""
    rng = np.random.default_rng(0)
    peso = np.cos(np.radians(np.linspace(-87.5, 87.5, 36))).reshape(-1, 1) * np.ones((1, 72))
    es_agua = rng.random((36, 72)) < 0.7
    for fraccion in (0.0, 0.1, 0.99, 1.0):
        hielo = es_agua & (rng.random((36, 72)) < fraccion)
        if fraccion == 1.0:
            hielo = es_agua.copy()
        F_prof = F.FLUJO_OCEANO_PROFUNDO
        compensacion = F_prof * (peso * hielo).sum() / max((peso * es_agua).sum(), 1e-12)
        aporte = np.where(hielo, F_prof, 0.0) - np.where(es_agua, compensacion, 0.0)
        assert abs((aporte * peso).sum()) < 1e-9
        assert compensacion <= F_prof + 1e-12   # nunca una retirada desproporcionada


# ---- estaciones ----

def test_estacion_no_depende_de_donde_empiece_el_invierno(monkeypatch):
    """orbita.estacion (v2.4.1) funciona aunque el invierno empiece cerca
    de 360 grados; la version anterior lo suponia por debajo de 270."""
    import orbita
    for inicio in (0.3, 4.9, 6.0):
        monkeypatch.setattr(orbita, "AV_INVIERNO", inicio)
        nombres = [orbita.estacion((inicio + k * math.pi / 2 + 0.01) % (2 * math.pi)) for k in range(4)]
        assert nombres == ["Invierno", "Primavera", "Verano", "Otoño"]


# ---- duracion del dia (v2.4.2) ----

def test_angulo_horario_identico_con_dia_de_24h():
    import geometria as G
    for h in np.arange(0, 24, 0.25):
        for lon in (-177.5, -0.5, 0.0, 92.5):
            antes = (h - 12) * math.pi / 12 + lon * math.pi / 180   # formula de la v2.4.1
            assert G.angulo_horario(h, lon) == antes                 # igualdad exacta, bit a bit


@pytest.mark.parametrize("dia_s", [64800, 61200, 90000])   # 18 h, 17 h, 25 h
def test_angulo_horario_da_una_vuelta_por_dia(monkeypatch, dia_s):
    import geometria as G
    monkeypatch.setattr(G, "MEDIO_DIA_H", dia_s / 3600 / 2)
    horas = dia_s / 3600
    assert G.angulo_horario(0.0) == pytest.approx(-math.pi)
    assert G.angulo_horario(horas / 2) == pytest.approx(0.0)          # mediodia en el meridiano 0
    assert G.angulo_horario(horas) == pytest.approx(math.pi)
    # el sol avanza 360/horas grados por hora
    assert G.angulo_horario(1.0) - G.angulo_horario(0.0) == pytest.approx(2 * math.pi / horas)


def test_simulacion_rechaza_dia_que_no_es_multiplo_del_paso(monkeypatch):
    monkeypatch.setattr(F, "ROTACION_PERIODO", 63360)   # 17,6 h: 70,4 pasos de 900 s
    with pytest.raises(ValueError, match="multiplo exacto del paso"):
        F.simular_fase2b(None, None, None, None, None, None, None, None, paso_tiempo=900)

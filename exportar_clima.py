# exportar_clima.py -- Fase 4 de M3N (v2.4): exportacion estructurada del
# clima para H3N y otras herramientas.
#
# Uso:  python exportar_clima.py
#
# Simula (o recupera de la cache) el clima del MAPA ACTIVO de C3N y lo
# escribe en:
#     ~/Documentos/B3N/clima_activo/clima_activo_m3n.json
# (simetrico al puente de entrada, ~/Documentos/B3N/mapa_activo/).
#
# ================================================================
# FORMATO "m3n-clima", version_formato 1
# ================================================================
# Un unico JSON (decision de Carlos, 03/10/2026). Las series grandes van
# COMPRIMIDAS dentro del propio JSON, porque en texto plano los datos
# horarios pesarian cientos de MB:
#
#   {"codificacion": "int16-centesimas-zlib-base64",
#    "forma": [dias o horas, 36, 72], "datos": "<base64>"}
#
#   -> base64 -> zlib -> enteros de 16 bits (little-endian, orden C:
#      ultimo indice = columna) -> dividir entre 100 = grados C (o metros
#      en el caso del hielo). El valor -32768 significa "sin dato".
#
# Indices: fila 0 = la mas al norte (87.5 N), columna 0 = 177.5 O (igual
# que la rejilla de M3N). Dia d = indice 0..num_dias-1 (dia d+1 del año; 326 dias con el dia de 19,84 h).
# Hora k = dia*24 + hora, instantanea a la hora en punto del MERIDIANO 0.
# v2.4.3: son HORAS DE P3N (1/24 del dia solar, DURACION_HORA s).
# Hora solar local de una celda = hora del meridiano 0 + longitud/15.
#
# La documentacion completa del formato esta en DISENO_FASE_4.md.
#
# 3.16.0 (DISENO_FASE_4.md §6; DISENO_FASE_6.3.md §6.18): exportacion de la CLIMATOLOGIA con el nucleo
# dinamico (I16), resultado de "python clima_dinamico.py":
#     python exportar_clima.py --dinamico            # P3N, carpeta outputs/clima_dinamico/p3n_dos_capas
#     python exportar_clima.py --dinamico --carpeta RUTA
# Mismo formato (m3n-clima v1) con CLAVES NUEVAS (las existentes no cambian): medias de N años, lluvia,
# nieve, evaporacion, viento y humedad de la capa baja, extremos absolutos de los N años, desviacion entre
# años, campos anuales y la climatologia (años, error, deriva, procedencia). Mientras no sea el clima oficial
# (4.0.0, Fase 6.4), se escribe en ~/Documentos/B3N/clima_activo/clima_dinamico_m3n.json, NO en el que lee H3N;
# con --activo se escribe en el oficial (lo que hara la 4.0.0).

import hashlib
import json
import os
import subprocess
import sys
import base64
import zlib
from datetime import datetime

import numpy as np

from parametros import *
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS, FILAS, COLUMNAS
from fase1_geografia import TIERRA, ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from orbita import info_estaciones, anomalia_media
from cache_simulacion import precalcular_orbita_cacheada, simular_fase2b_cacheada
import fase2b_atmosfera as F2B

FORMATO = "m3n-clima"
VERSION_FORMATO = 1
SIN_DATO = -32768

CARPETA_SALIDA = os.path.expanduser(os.path.join("~", "Documentos", "B3N", "clima_activo"))
RUTA_SALIDA = os.path.join(CARPETA_SALIDA, "clima_activo_m3n.json")
RUTA_SALIDA_DINAMICO = os.path.join(CARPETA_SALIDA, "clima_dinamico_m3n.json")     # 3.16.0, hasta la 4.0.0

# 3.16.0: unidades de cada clave de las series (las de antes: C, y m en el hielo)
UNIDADES = {
    "aire2m_min": "C", "aire2m_media": "C", "aire2m_max": "C", "superficie_min": "C", "superficie_media": "C",
    "superficie_max": "C", "hielo_espesor_m": "m", "aire2m": "C", "superficie": "C",
    "precipitacion_mm_dia": "mm/dia (kg/m2 por dia terrestre de 86400 s)", "nieve_mm_dia": "mm/dia de agua",
    "evaporacion_mm_dia": "mm/dia", "viento_u_baja_m_s": "m/s (hacia el este +)",
    "viento_v_baja_m_s": "m/s (hacia el norte +)", "viento_rapidez_baja_m_s": "m/s",
    "humedad_especifica_baja_g_kg": "g/kg", "humedad_relativa_baja_pct": "% (respecto al agua liquida)",
    "aire2m_min_absoluta": "C", "aire2m_max_absoluta": "C", "superficie_min_absoluta": "C",
    "superficie_max_absoluta": "C", "viento_rapidez_baja_max_m_s": "m/s",
    "aire2m_media_desviacion_entre_anos": "C", "precipitacion_desviacion_entre_anos_mm_dia": "mm/dia",
    "agua_precipitable_kg_m2": "kg/m2",
    "nieve_suelo_cm_agua": "cm de agua equivalente (= kg/m2 / 10; la nieve en tierra llega hasta 1000 kg/m2)",
}


def serie(array):
    """Codifica un array en centesimas, int16, comprimido y en base64."""
    a = np.asarray(array, dtype=np.float64)
    finito = np.isfinite(a)
    centesimas = np.round(np.where(finito, a, 0.0) * 100)
    if np.any(finito & (np.abs(centesimas) > 32767)):
        raise ValueError("valor fuera del rango de int16 en centesimas (|x| > 327.67)")
    enteros = np.where(finito, centesimas, SIN_DATO)
    datos = enteros.astype("<i2").tobytes()
    return {
        "codificacion": "int16-centesimas-zlib-base64",
        "forma": list(a.shape),
        "datos": base64.b64encode(zlib.compress(datos, 6)).decode("ascii"),
    }


def leer_serie(s):
    """Inversa de serie() -- la usan las pruebas y cualquier lector en Python."""
    crudo = np.frombuffer(zlib.decompress(base64.b64decode(s["datos"])), dtype="<i2").reshape(s["forma"])
    return np.where(crudo == SIN_DATO, np.nan, crudo / 100.0)


def version_git():
    try:
        r = subprocess.run(["git", "describe", "--tags", "--always", "--dirty"],
                           capture_output=True, text=True, timeout=5, cwd=os.path.dirname(os.path.abspath(__file__)))
        return r.stdout.strip() or "sin git"
    except (OSError, subprocess.SubprocessError):
        return "sin git"


def huella_mapa(tipo, altitud):
    """Identifica la geografia exacta con la que se simulo."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(tipo, dtype=np.int8).tobytes())
    h.update(np.ascontiguousarray(np.round(altitud, 3), dtype=np.float64).tobytes())
    return h.hexdigest()[:16]


def astronomia_diaria(datos_orbita, num_dias, paso_tiempo):
    """Declinacion solar y flujo en lo alto de la atmosfera a mediodia
    (hora 12 del meridiano 0) de cada dia, y distancia a S3N."""
    pasos_dia = round(ROTACION_PERIODO / paso_tiempo)
    decl, flujo, dist = [], [], []
    for d in range(num_dias):
        toa, de, _ = datos_orbita[d * pasos_dia + pasos_dia // 2]
        decl.append(round(float(np.degrees(de)), 4))
        flujo.append(round(float(toa), 2))
        dist.append(round(float(np.sqrt(S3N_LUMINOSIDAD / (4 * PI * toa)) / 1.495978707e11), 6))
    return decl, flujo, dist


def fechas_clave(num_dias):
    am_inicio = anomalia_media(1, 0)
    t_peri = ((2 * PI - am_inicio) % (2 * PI)) / (2 * PI) * ORBITA_PERIODO
    t_afe = ((PI - am_inicio) % (2 * PI)) / (2 * PI) * ORBITA_PERIODO
    hitos = [
        {"nombre": "Perihelio", "dia": int(t_peri // ROTACION_PERIODO) % num_dias},
        {"nombre": "Afelio", "dia": int(t_afe // ROTACION_PERIODO) % num_dias},
    ]
    # Las estaciones de M3N son las del hemisferio norte: su inicio es un
    # solsticio o un equinoccio.
    nombres = {"Invierno": "Solsticio de invierno (HN)", "Primavera": "Equinoccio de primavera (HN)",
               "Verano": "Solsticio de verano (HN)", "Otoño": "Equinoccio de otoño (HN)"}
    for est in info_estaciones():
        hitos.append({"nombre": nombres.get(est["estacion"], est["estacion"]),
                      "dia": int(np.floor(est["dia_inicio"] - 1)) % num_dias})
    return sorted(hitos, key=lambda h: h["dia"])


def construir_exportacion(r, tipo, altitud, nombre_mapa, datos_orbita, paso_tiempo):
    num_dias = r["reg_media"].shape[0]
    if r.get("horario_aire2m") is None:
        raise ValueError("La simulacion no tiene registro horario (el paso de tiempo no divide una hora de P3N)")
    horas_dia = HORAS_POR_DIA                       # v2.4.3: horas de P3N (1/24 del dia)
    grados_hora = 360 / horas_dia
    decl, flujo, dist = astronomia_diaria(datos_orbita, num_dias, paso_tiempo)
    interruptores = {k: bool(v) for k, v in F2B.INTERRUPTORES_FASE2B.items()}
    return {
        "formato": FORMATO,
        "version_formato": VERSION_FORMATO,
        "generado": datetime.now().isoformat(timespec="seconds"),
        "m3n": {"version": version_git(), "interruptores": interruptores},
        "mapa": {"nombre": nombre_mapa, "huella": huella_mapa(tipo, altitud)},
        "rejilla": {
            "filas": FILAS, "columnas": COLUMNAS,
            "latitudes": [float(x) for x in LATITUDES_GRADOS],
            "longitudes": [float(x) for x in LONGITUDES_GRADOS],
            "nota": "fila 0 = norte; valores en el centro de cada celda",
        },
        "tiempo": {
            "dias": num_dias, "horas_por_dia": horas_dia,
            "duracion_dia_s": ROTACION_PERIODO, "duracion_hora_s": DURACION_HORA,
            "duracion_año_dias": ORBITA_PERIODO / ROTACION_PERIODO,
            "referencia_hora": f"meridiano 0 (hora solar local = hora + longitud/{grados_hora:g})",
            "fechas_clave": fechas_clave(num_dias),
        },
        "astronomia": {
            "declinacion_solar_grados": decl,
            "flujo_toa_W_m2": flujo,
            "distancia_S3N_UA": dist,
            "nota": f"valores a la hora {horas_dia / 2:g} del meridiano 0 de cada dia; "
                    f"longitud subsolar = -(hora - {horas_dia / 2:g}) * {grados_hora:g}",
        },
        "parametros": {
            "masa_S3N_soles": S3N_MASAS_SOLARES, "luminosidad_W": S3N_LUMINOSIDAD,
            "semieje_m": SEMIEJE_MAYOR, "excentricidad": ORBITA_EXCENTRICIDAD,
            "inclinacion_axial_grados": INCLINACION_AXIAL,
            # v3.1: los D que ha usado de verdad la simulacion (escalados con la rotacion)
            "D_atmosfera": r.get("D_usados", {}).get("atmosfera", F2B.D_ATMOSFERA),
            "D_oceano": r.get("D_usados", {}).get("oceano", F2B.D_OCEANO),
            "D_esquema": r.get("D_usados", {}).get("esquema", "desconocido"),
        },
        "simulacion": {
            "anos_hasta_convergencia": int(r["anos"]), "saltos": int(r["saltos"]),
            "balance_energia": float(r["energia"]["diferencia_relativa"]),
            "paso_tiempo_s": paso_tiempo,
        },
        "celdas": {
            "tipo": ["tierra" if t == TIERRA else "agua" for t in np.asarray(tipo).ravel()],
            "altitud_m": [round(float(a), 1) for a in np.asarray(altitud).ravel()],
            "nieve_permanente_posible": [bool(x) for x in np.asarray(r["nieve_permanente_posible"]).ravel()],
        },
        "diario": {
            "aire2m_min": serie(r["reg_min"]), "aire2m_media": serie(r["reg_media"]), "aire2m_max": serie(r["reg_max"]),
            "superficie_min": serie(r["suelo_min"]), "superficie_media": serie(r["suelo_media"]),
            "superficie_max": serie(r["suelo_max"]),
            "hielo_espesor_m": serie(r["hielo_espesor"]),
        },
        "horario": {
            "aire2m": serie(r["horario_aire2m"]),
            "superficie": serie(r["horario_superficie"]),
        },
    }


def construir_exportacion_dinamica(rc, tipo, altitud, nombre_mapa, datos_orbita, paso_tiempo):
    """3.16.0: la climatologia de clima_dinamico.py (resultado.pkl, media de N años con I16) en m3n-clima v1.
    Las claves de siempre tienen el mismo significado (ahora, medias de N años); las nuevas, en 'diario'
    (lluvia, viento, humedad), 'extremos', 'variabilidad', 'anual', 'climatologia' y 'unidades'."""
    c = rc["climatologia"]
    if rc.get("mapa") is not None and rc["mapa"]["huella"] != huella_mapa(tipo, altitud):
        raise ValueError(f"La climatologia se hizo con otro mapa ({rc['mapa']['nombre']}, huella "
                         f"{rc['mapa']['huella']}); el mapa activo es {nombre_mapa} ({huella_mapa(tipo, altitud)}). "
                         "Vuelve a lanzar clima_dinamico.py con el mapa activo.")
    es_tierra = np.asarray(tipo) == TIERRA
    r = dict(rc)
    r["anos"] = int(c["equilibrio"]["anos"])
    r["saltos"] = 0
    r["nieve_permanente_posible"] = F2B.nieve_permanente_posible(np.asarray(rc["reg_media"]), es_tierra)
    e = rc.get("energia")
    if not isinstance(e, dict) or "diferencia_relativa" not in e:
        ultimo = c["por_ano"][-1]["energia"] if c.get("por_ano") else {}
        r["energia"] = {"diferencia_relativa": float(abs(ultimo.get("cierre_relativo", float("nan"))))}
    out = construir_exportacion(r, tipo, altitud, nombre_mapa, datos_orbita, paso_tiempo)
    # con I16 la atmosfera no difunde (la mueve el nucleo); el oceano, si (D de la Fase 2.3, hasta la 6.5)
    out["parametros"].update({"D_atmosfera": None, "D_oceano": float(F2B.D_OCEANO_V30),
                              "D_esquema": "nucleo dinamico (I16); oceano por difusion"})
    a = rc["agua"]
    out["simulacion"].update({
        "nivel_modelo": "completo: nucleo dinamico (I16), 20 capas, ciclo del agua; sin nubes (Fase 5.2)",
        "anos_hasta_convergencia_significado": "años del modelo con nucleo hasta el equilibrio (fase6_equilibrio)",
    })
    diario = out["diario"]
    for clave, origen in (("precipitacion_mm_dia", a.get("diario_precipitacion")),
                          ("nieve_mm_dia", a.get("diario_nieve")),
                          ("evaporacion_mm_dia", a.get("diario_evaporacion")),
                          ("viento_u_baja_m_s", rc.get("viento_u_baja")), ("viento_v_baja_m_s", rc.get("viento_v_baja")),
                          ("viento_rapidez_baja_m_s", rc.get("viento_rapidez_baja")),
                          ("humedad_especifica_baja_g_kg", rc.get("humedad_especifica_baja")),
                          ("humedad_relativa_baja_pct", rc.get("humedad_relativa_baja"))):
        if origen is not None:
            diario[clave] = serie(origen)
    ext = rc.get("extremos", {})
    out["extremos"] = {k: serie(ext[o]) for k, o in (
        ("aire2m_min_absoluta", "reg_min"), ("aire2m_max_absoluta", "reg_max"),
        ("superficie_min_absoluta", "suelo_min"), ("superficie_max_absoluta", "suelo_max"),
        ("viento_rapidez_baja_max_m_s", "viento_rapidez_baja")) if o in ext}
    out["extremos"]["nota"] = "por dia del año, el minimo de los minimos y el maximo de los maximos de los N años"
    d = rc.get("desviacion_entre_anos", {})
    out["variabilidad"] = {"nota": "desviacion tipica entre años (n - 1)"}
    if "reg_media" in d:
        out["variabilidad"]["aire2m_media_desviacion_entre_anos"] = serie(d["reg_media"])
    if "precipitacion" in d:
        out["variabilidad"]["precipitacion_desviacion_entre_anos_mm_dia"] = serie(d["precipitacion"])
    out["anual"] = {k: serie(a[o]) for k, o in (("precipitacion_mm_dia", "precipitacion"), ("evaporacion_mm_dia", "evaporacion"),
                                                  ("agua_precipitable_kg_m2", "agua_precipitable")) if o in a}
    if "nieve_media" in a:          # en cm de agua: en kg/m2 pasaria del rango de las centesimas (|x| <= 327,67)
        out["anual"]["nieve_suelo_cm_agua"] = serie(np.asarray(a["nieve_media"]) / 10.0)
    dv = c.get("deriva") or {}
    out["climatologia"] = {
        "anos_promediados": int(c["anos_promediados"]), "anos_objetivo": int(c.get("anos_objetivo", c["anos_promediados"])),
        "exploratoria": bool(c.get("exploratoria", False)),
        "error_media_aire2m_max_K": float(c["error_T_max_K"]), "error_media_precipitacion_max_mm_dia": float(c["error_P_max_mm_dia"]),
        "deriva": {"hay_deriva": bool(dv.get("hay_deriva", False)), "metodo": dv.get("metodo"),
                   "aire2m_global_K_ano": dv.get("global_T", {}).get("pendiente")} if dv else None,
        "equilibrio": {"anos": int(c["equilibrio"]["anos"]), "convergido": bool(c["equilibrio"]["convergido"])},
        "procedencia": c.get("procedencia", []),
    }
    out["unidades"] = UNIDADES
    return out


def escribir_json_atomico(ruta, contenido):
    os.makedirs(os.path.dirname(ruta), exist_ok=True)
    temporal = ruta + ".tmp"
    with open(temporal, "w", encoding="utf-8") as f:
        json.dump(contenido, f, ensure_ascii=False)
    os.replace(temporal, ruta)


if __name__ == "__main__":
    import argparse
    import pickle
    from puente_c3n import cargar_mapa_activo_de_c3n
    from temperatura import PASO_TIEMPO

    ap = argparse.ArgumentParser(description="Exporta el clima de M3N para H3N (formato m3n-clima v1)")
    ap.add_argument("--dinamico", action="store_true", help="3.16.0: la climatologia con nucleo dinamico (I16)")
    ap.add_argument("--carpeta", default=None, help="carpeta de clima_dinamico.py (por defecto la de P3N)")
    ap.add_argument("--activo", action="store_true", help="con --dinamico: escribir en el clima que lee H3N")
    args = ap.parse_args()
    if MODO_TIERRA:
        sys.exit("exportar_clima.py es para P3N (quita M3N_MODO=tierra)")
    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo, altitud, nombre_mapa = cargar_mapa_activo_de_c3n()
    if args.dinamico:
        carpeta = args.carpeta or os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "clima_dinamico",
                                               "p3n_dos_capas")
        f_res = os.path.join(carpeta, "resultado.pkl")
        if not os.path.exists(f_res):
            sys.exit(f"Falta {f_res}: hay que terminar antes 'python clima_dinamico.py'.")
        with open(f_res, "rb") as f:
            rc = pickle.load(f)
        if rc.get("mapa") is None:
            sys.exit("Esta climatologia no dice con que mapa se hizo (es anterior a la 3.16.0): no se puede "
                     "comprobar que sea la del mapa activo. Vuelve a lanzar clima_dinamico.py (si esta terminada, "
                     "solo rehace el resumen y el resultado).")
        exportacion = construir_exportacion_dinamica(rc, tipo, altitud, nombre_mapa, datos_orbita, PASO_TIEMPO)
        ruta = RUTA_SALIDA if args.activo else RUTA_SALIDA_DINAMICO
        escribir_json_atomico(ruta, exportacion)
        print(f"Climatologia con nucleo dinamico exportada en {ruta} ({os.path.getsize(ruta) / 1e6:.1f} MB)"
              + ("" if args.activo else " (no es el clima que lee H3N: eso sera la 4.0.0)"))
        sys.exit(0)
    print(f"Simulando mapa '{nombre_mapa}' (si ya esta en cache, es instantaneo)...")
    r = simular_fase2b_cacheada(datos_orbita, tipo, altitud, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
                                PROFUNDIDAD_OPTICA, float(D_DIFUSION_REFERENCIA), nombre_mapa=nombre_mapa)
    exportacion = construir_exportacion(r, tipo, altitud, nombre_mapa, datos_orbita, PASO_TIEMPO)
    escribir_json_atomico(RUTA_SALIDA, exportacion)
    print(f"Clima exportado en {RUTA_SALIDA} ({os.path.getsize(RUTA_SALIDA) / 1e6:.1f} MB)")

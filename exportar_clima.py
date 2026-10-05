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
# La documentacion completa del formato esta en DISENO_FASE4.md.

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


def escribir_json_atomico(ruta, contenido):
    os.makedirs(os.path.dirname(ruta), exist_ok=True)
    temporal = ruta + ".tmp"
    with open(temporal, "w", encoding="utf-8") as f:
        json.dump(contenido, f, ensure_ascii=False)
    os.replace(temporal, ruta)


if __name__ == "__main__":
    from puente_c3n import cargar_mapa_activo_de_c3n
    from temperatura import PASO_TIEMPO

    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo, altitud, nombre_mapa = cargar_mapa_activo_de_c3n()
    print(f"Simulando mapa '{nombre_mapa}' (si ya esta en cache, es instantaneo)...")
    r = simular_fase2b_cacheada(datos_orbita, tipo, altitud, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
                                PROFUNDIDAD_OPTICA, float(D_DIFUSION_REFERENCIA), nombre_mapa=nombre_mapa)
    exportacion = construir_exportacion(r, tipo, altitud, nombre_mapa, datos_orbita, PASO_TIEMPO)
    escribir_json_atomico(RUTA_SALIDA, exportacion)
    print(f"Clima exportado en {RUTA_SALIDA} ({os.path.getsize(RUTA_SALIDA) / 1e6:.1f} MB)")

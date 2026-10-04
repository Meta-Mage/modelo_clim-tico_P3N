# consulta_punto.py -- clima de un punto concreto de P3N a lo largo del año.
#
# Uso:  python consulta_punto.py LAT LON      (sin argumentos, las pregunta)
#
# Busca la celda de la rejilla (5 x 5 grados) que contiene el punto y
# muestra, para 12 dias repartidos por el año, la minima, media y maxima
# del AIRE A 2 m, la temperatura media de la superficie y, en el oceano,
# el espesor del hielo marino (Fase 3). Usa la misma simulacion (y la
# misma cache) que el resto de herramientas.

import os
import sys
from datetime import date

import numpy as np

from parametros import *
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS
from fase1_geografia import TIERRA, ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from cache_simulacion import precalcular_orbita_cacheada, simular_fase2b_cacheada

NUM_DIAS_MUESTRA = 12
CARPETA_RESULTADOS = os.path.join("outputs", "consultas_punto")


def celda_mas_cercana(latitud_grados, longitud_grados):
    fila = int(np.argmin(np.abs(LATITUDES_GRADOS - latitud_grados)))
    dlon = (LONGITUDES_GRADOS - longitud_grados + 180) % 360 - 180   # con envolvido en +-180
    columna = int(np.argmin(np.abs(dlon)))
    return fila, columna


def consultar_punto(latitud_grados, longitud_grados, r, tipo_superficie, altitud_metros):
    fila, columna = celda_mas_cercana(latitud_grados, longitud_grados)
    es_tierra = tipo_superficie[fila, columna] == TIERRA
    num_dias = r["reg_media"].shape[0]
    indices = np.linspace(0, num_dias - 1, NUM_DIAS_MUESTRA).astype(int)

    lineas = [
        f"Coordenadas pedidas: lat {latitud_grados}, lon {longitud_grados}",
        f"Celda: fila {fila}, columna {columna} -> centro en lat {LATITUDES_GRADOS[fila]:.1f}, "
        f"lon {LONGITUDES_GRADOS[columna]:.1f} | "
        + (f"tierra, {altitud_metros[fila, columna]:.0f} m de altitud media" if es_tierra else "oceano"),
        "",
        "Aire a 2 m (minima / media / maxima del dia) y superficie (media del dia), en C"
        + ("" if es_tierra else "; hielo marino en metros"),
        f"{'Dia':>5} | {'Minima':>8} | {'Media':>8} | {'Maxima':>8} | {'Superficie':>10}" + ("" if es_tierra else f" | {'Hielo':>6}"),
    ]
    lineas.append("-" * len(lineas[-1]))
    for d in indices:
        linea = (f"{d + 1:>5} | {r['reg_min'][d, fila, columna]:>8.2f} | {r['reg_media'][d, fila, columna]:>8.2f} | "
                 f"{r['reg_max'][d, fila, columna]:>8.2f} | {r['suelo_media'][d, fila, columna]:>10.2f}")
        if not es_tierra:
            linea += f" | {r['hielo_espesor'][d, fila, columna]:>6.2f}"
        lineas.append(linea)
    lineas.append("-" * len(lineas[-2]))
    aire = r["reg_media"][:, fila, columna]
    lineas.append(f"Año completo (aire a 2 m): media {aire.mean():.2f} C | minima absoluta "
                  f"{r['reg_min'][:, fila, columna].min():.2f} C | maxima absoluta {r['reg_max'][:, fila, columna].max():.2f} C")
    texto = "\n".join(lineas)
    print(texto)
    return texto


def guardar_resultado(texto, latitud_grados, longitud_grados, nombre_mapa):
    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)
    nombre_archivo = f"consulta_lat{latitud_grados}_lon{longitud_grados}_{nombre_mapa}_{date.today().isoformat()}.txt"
    ruta = os.path.join(CARPETA_RESULTADOS, nombre_archivo)
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    print(f"\nGuardado en {ruta}")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        latitud_grados = float(sys.argv[1])
        longitud_grados = float(sys.argv[2])
    else:
        latitud_grados = float(input("Latitud (grados, -90 a 90): "))
        longitud_grados = float(input("Longitud (grados, -180 a 180): "))

    from puente_c3n import cargar_mapa_activo_de_c3n
    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print(f"Simulando mapa '{nombre_mapa}' (si ya esta en cache, es instantaneo)...")
    r = simular_fase2b_cacheada(datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, ALBEDO_POR_TIPO,
                                INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA, float(D_DIFUSION_REFERENCIA), nombre_mapa=nombre_mapa)
    print(f"Convergencia: {r['anos']} año(s)")
    if r["anos"] >= 50:
        print("AVISO: se alcanzo el limite de 50 años sin confirmar convergencia -- puede no ser el equilibrio.")
    print()
    texto = consultar_punto(latitud_grados, longitud_grados, r, tipo_superficie, altitud_metros)
    guardar_resultado(texto, latitud_grados, longitud_grados, nombre_mapa)

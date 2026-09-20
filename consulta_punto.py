import os
import sys
from datetime import date

import numpy as np

from parametros import *
from temperatura import precalcular_orbita
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS
from fase1_geografia import (
    mapa_falso_todo_tierra, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
    simular_rejilla_geografia_con_registro,
)

NUM_DIAS_MUESTRA = 12
CARPETA_RESULTADOS = os.path.join("outputs", "consultas_punto")


def celda_mas_cercana(latitud_grados, longitud_grados):
    fila = int(np.argmin(np.abs(LATITUDES_GRADOS - latitud_grados)))
    columna = int(np.argmin(np.abs(LONGITUDES_GRADOS - longitud_grados)))
    return fila, columna


def consultar_punto(latitud_grados, longitud_grados, registro_minima, registro_media, registro_maxima):
    fila, columna = celda_mas_cercana(latitud_grados, longitud_grados)
    lat_real = LATITUDES_GRADOS[fila]
    lon_real = LONGITUDES_GRADOS[columna]

    num_dias_totales = registro_media.shape[0]
    indices_muestra = np.linspace(0, num_dias_totales - 1, NUM_DIAS_MUESTRA).astype(int)

    lineas = []
    lineas.append(f"Coordenadas pedidas: lat {latitud_grados}, lon {longitud_grados}")
    lineas.append(f"Celda mas cercana: fila {fila}, columna {columna} "
                   f"-> lat real {lat_real:.2f}, lon real {lon_real:.2f}")
    lineas.append("")
    lineas.append(f"{'Dia del año':>12} | {'% del año':>10} | "
                   f"{'Minima (C)':>11} | {'Media (C)':>10} | {'Maxima (C)':>11}")
    lineas.append("-" * 68)

    minimas_muestra = []
    medias_muestra = []
    maximas_muestra = []
    for indice_dia in indices_muestra:
        minima = registro_minima[indice_dia, fila, columna]
        media = registro_media[indice_dia, fila, columna]
        maxima = registro_maxima[indice_dia, fila, columna]
        minimas_muestra.append(minima)
        medias_muestra.append(media)
        maximas_muestra.append(maxima)
        porcentaje = 100 * indice_dia / (num_dias_totales - 1)
        lineas.append(f"{indice_dia:>12} | {porcentaje:>9.1f}% | "
                       f"{minima:>11.2f} | {media:>10.2f} | {maxima:>11.2f}")

    minimas_muestra = np.array(minimas_muestra)
    medias_muestra = np.array(medias_muestra)
    maximas_muestra = np.array(maximas_muestra)
    lineas.append("-" * 68)
    lineas.append(f"De los {NUM_DIAS_MUESTRA} dias mostrados: "
                   f"minima absoluta {minimas_muestra.min():.2f} C | "
                   f"media de las medias {medias_muestra.mean():.2f} C | "
                   f"maxima absoluta {maximas_muestra.max():.2f} C")

    texto = "\n".join(lineas)
    print(texto)
    return texto


def guardar_resultado(texto, latitud_grados, longitud_grados, nombre_mapa):
    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)
    fecha = date.today().isoformat()
    nombre_archivo = f"consulta_lat{latitud_grados}_lon{longitud_grados}_{nombre_mapa}_{fecha}.txt"
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

    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    from puente_c3n import cargar_mapa_activo_de_c3n
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print(f"Simulando mapa '{nombre_mapa}'...")
    _, anos_convergencia, registro_minima, registro_media, registro_maxima = simular_rejilla_geografia_con_registro(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
    )
    print(f"Convergencia: {anos_convergencia} año(s)")
    if anos_convergencia >= 50:
        print("AVISO: se alcanzo el limite de 50 años sin confirmar convergencia real "
              "(probablemente por la enorme inercia termica del agua) -- este resultado "
              "puede no ser el equilibrio verdadero.")
    print()

    texto = consultar_punto(latitud_grados, longitud_grados, registro_minima, registro_media, registro_maxima)
    guardar_resultado(texto, latitud_grados, longitud_grados, nombre_mapa)

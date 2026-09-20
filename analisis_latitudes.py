import os
from datetime import date

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametros import *
from temperatura import precalcular_orbita
from orbita import info_estaciones
from rejilla import LATITUDES_GRADOS
from fase1_geografia import (
    TIERRA, mapa_falso_todo_tierra, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
    simular_rejilla_geografia_con_registro,
)

LATITUDES_REFERENCIA_GRADOS = [60, 30, 0, -30, -60]
CARPETA_TABLAS = os.path.join("outputs", "tablas_latitud")
CARPETA_GRAFICOS = os.path.join("outputs", "evolucion_por_latitud")


def fila_para_latitud(latitud_objetivo):
    return int(np.argmin(np.abs(LATITUDES_GRADOS - latitud_objetivo)))


def obtener_dias_especiales(num_dias):
    """
    8 dias de referencia: inicio de cada estacion (solsticio/equinoccio,
    calculado por info_estaciones() respetando la excentricidad orbital)
    y el punto medio (en tiempo real, no en angulo) de cada estacion.
    """
    info = info_estaciones()
    dias_especiales = {}
    for estacion in info:
        nombre = estacion["estacion"]
        dia_inicio = estacion["dia_inicio"]
        duracion = estacion["duracion"]
        dia_medio = dia_inicio + duracion / 2

        indice_inicio = int(round(dia_inicio - 1)) % num_dias
        indice_medio = int(round(dia_medio - 1)) % num_dias

        dias_especiales[f"Inicio de {nombre}"] = indice_inicio
        dias_especiales[f"Mitad de {nombre}"] = indice_medio
    return dias_especiales


def estadisticas_fila(valores_fila, tipo_superficie_fila):
    """
    valores_fila: array (columnas,) de temperaturas de una fila, un dia.
    tipo_superficie_fila: array (columnas,) TIERRA/AGUA de esa misma fila.
    Usa solo celdas de tierra si hay alguna en la fila; si la fila es
    toda agua, usa las celdas de agua en su lugar.
    """
    mascara_tierra = tipo_superficie_fila == TIERRA
    if mascara_tierra.any():
        return valores_fila[mascara_tierra], "tierra"
    return valores_fila, "agua"


def media_diaria_por_latitud(registro_media, tipo_superficie):
    """
    Para cada latitud de referencia, la media (tierra si hay, si no agua)
    a lo largo de todos los dias registrados. Forma: (latitudes, dias).
    """
    resultado = np.zeros((len(LATITUDES_REFERENCIA_GRADOS), registro_media.shape[0]))
    for i, latitud_objetivo in enumerate(LATITUDES_REFERENCIA_GRADOS):
        fila = fila_para_latitud(latitud_objetivo)
        tipo_fila = tipo_superficie[fila, :]
        for dia in range(registro_media.shape[0]):
            valores, _ = estadisticas_fila(registro_media[dia, fila, :], tipo_fila)
            resultado[i, dia] = valores.mean()
    return resultado


def generar_tabla(registro_minima, registro_media, registro_maxima, tipo_superficie, dias_especiales, nombre_mapa):
    encabezado = f"{'Latitud':>8} | {'Dia de referencia':<20} | {'Tipo':<7} | {'Minima (C)':>11} | {'Media (C)':>10} | {'Maxima (C)':>11}"
    lineas = [encabezado, "-" * len(encabezado)]

    for latitud_objetivo in LATITUDES_REFERENCIA_GRADOS:
        fila = fila_para_latitud(latitud_objetivo)
        tipo_fila = tipo_superficie[fila, :]
        for nombre_dia, indice_dia in dias_especiales.items():
            minima_vals, tipo_usado = estadisticas_fila(registro_minima[indice_dia, fila, :], tipo_fila)
            media_vals, _ = estadisticas_fila(registro_media[indice_dia, fila, :], tipo_fila)
            maxima_vals, _ = estadisticas_fila(registro_maxima[indice_dia, fila, :], tipo_fila)
            lineas.append(
                f"{latitud_objetivo:>8} | {nombre_dia:<20} | {tipo_usado:<7} | "
                f"{minima_vals.min():>11.2f} | {media_vals.mean():>10.2f} | {maxima_vals.max():>11.2f}"
            )
        lineas.append("-" * len(encabezado))

    texto = "\n".join(lineas)
    print(texto)

    os.makedirs(CARPETA_TABLAS, exist_ok=True)
    fecha = date.today().isoformat()
    ruta = os.path.join(CARPETA_TABLAS, f"tabla_latitud_{nombre_mapa}_{fecha}.txt")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    print(f"\nTabla guardada en {ruta}")


def generar_grafico(registro_media, tipo_superficie, dias_especiales, nombre_mapa):
    medias_por_latitud = media_diaria_por_latitud(registro_media, tipo_superficie)
    dias = np.arange(medias_por_latitud.shape[1])

    colores = {
        60: "#01579b", 30: "#4fc3f7", 0: "#616161", -30: "#ff8a65", -60: "#bf360c",
    }

    plt.figure(figsize=(12, 6))
    for i, latitud_objetivo in enumerate(LATITUDES_REFERENCIA_GRADOS):
        plt.plot(dias, medias_por_latitud[i], color=colores[latitud_objetivo],
                  linewidth=2, label=f"Latitud {latitud_objetivo}°")

    for indice_dia in sorted(set(dias_especiales.values())):
        plt.axvline(indice_dia, color="gray", linestyle="--", linewidth=0.6, alpha=0.5)

    plt.xlabel("Dia del año")
    plt.ylabel("Temperatura media (°C)")
    plt.title(f"Evolucion anual de temperatura media por latitud ({nombre_mapa})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    os.makedirs(CARPETA_GRAFICOS, exist_ok=True)
    fecha = date.today().isoformat()
    ruta = os.path.join(CARPETA_GRAFICOS, f"evolucion_por_latitud_{nombre_mapa}_{fecha}.png")
    plt.savefig(ruta, dpi=150)
    print(f"Grafico guardado en {ruta}")


if __name__ == "__main__":
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    from puente_c3n import cargar_mapa_activo_de_c3n
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print("Simulando (mapa falso todo tierra, PROVISIONAL hasta el puente con C3N)...")
    _, anos_convergencia, registro_minima, registro_media, registro_maxima = simular_rejilla_geografia_con_registro(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
    )
    print(f"Convergencia: {anos_convergencia} año(s)")
    print()

    dias_especiales = obtener_dias_especiales(registro_media.shape[0])
    print("Dias de referencia detectados:")
    for nombre, indice in dias_especiales.items():
        print(f"  {nombre}: dia {indice}")
    print()

    generar_tabla(registro_minima, registro_media, registro_maxima, tipo_superficie, dias_especiales, nombre_mapa)
    generar_grafico(registro_media, tipo_superficie, dias_especiales, nombre_mapa)


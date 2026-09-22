import os
from datetime import date

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametros import *
from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from cache_simulacion import precalcular_orbita_cacheada, simular_rejilla_combinada_cacheada
from puente_c3n import cargar_mapa_activo_de_c3n
from rejilla import FILAS, COLUMNAS

# D de la difusion horizontal activa en esta herramienta. Pon esto a
# np.zeros((FILAS, COLUMNAS)) para volver a la fisica de Fase 1 pura
# (sin transporte de calor entre celdas), por ejemplo para comparar.
D_GRID_ACTIVO = np.full((FILAS, COLUMNAS), D_DIFUSION_REFERENCIA)

CARPETA_RESULTADOS = os.path.join("outputs", "mapas_calor")


if __name__ == "__main__":
    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print(f"Simulando mapa '{nombre_mapa}' (media anual por celda, promediando dia/noche y estaciones)...")
    _, anos_convergencia, _, registro_media, _ = simular_rejilla_combinada_cacheada(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA, D_GRID_ACTIVO,
        nombre_mapa=nombre_mapa,
    )
    print(f"Convergencia: {anos_convergencia} año(s)")
    if anos_convergencia >= 50:
        print("AVISO: se alcanzo el limite de 50 años sin confirmar convergencia real "
              "(probablemente por la enorme inercia termica del agua) -- este resultado "
              "puede no ser el equilibrio verdadero.")

    T_media_anual = registro_media.mean(axis=0)

    print(f"Temperatura media anual -- minima: {T_media_anual.min():.1f} C | maxima: {T_media_anual.max():.1f} C")
    print(f"Celdas por debajo de -100 C: {int((T_media_anual < -100).sum())} de {T_media_anual.size}")

    plt.figure(figsize=(14, 7))
    im = plt.imshow(T_media_anual, extent=[-180, 180, -90, 90], origin="upper", cmap="turbo",
                     aspect="auto", vmin=-60, vmax=60)
    plt.colorbar(im, label="Temperatura media anual (C), escala recortada a [-60, 60]")
    plt.xlabel("Longitud (grados)")
    plt.ylabel("Latitud (grados)")
    plt.title(f"Mapa de calor -- media anual por celda -- {nombre_mapa} (72x36)")
    plt.tight_layout()

    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)
    fecha = date.today().isoformat()
    ruta = os.path.join(CARPETA_RESULTADOS, f"mapa_calor_{nombre_mapa}_{fecha}.png")
    plt.savefig(ruta, dpi=150)
    print(f"Mapa guardado en {ruta}")

# mapa_calor.py -- mapa de la temperatura media anual del aire a 2 m.
#
# Uso:  python mapa_calor.py
#
# v2.4.1: escala de color DIVERGENTE azul-rojo centrada en 0 C y fija
# (-50 a +50 C), la misma que usa H3N; antes era "turbo" (un arcoiris,
# que engaña al ojo y no se lee con daltonismo -- Crameri et al. 2020).
# El hielo marino (celdas heladas al menos medio año) se marca rayado.

import os
from datetime import date

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametros import *
from fase1_geografia import TIERRA, ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from cache_simulacion import precalcular_orbita_cacheada, simular_fase2b_cacheada
from puente_c3n import cargar_mapa_activo_de_c3n

CARPETA_RESULTADOS = os.path.join("outputs", "mapas_calor")


if __name__ == "__main__":
    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print(f"Simulando mapa '{nombre_mapa}' (si ya esta en cache, es instantaneo)...")
    r = simular_fase2b_cacheada(datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, ALBEDO_POR_TIPO,
                                INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA, float(D_DIFUSION_REFERENCIA), nombre_mapa=nombre_mapa)
    print(f"Convergencia: {r['anos']} año(s)")
    if r["anos"] >= 50:
        print("AVISO: se alcanzo el limite de 50 años sin confirmar convergencia -- puede no ser el equilibrio.")

    T_media_anual = r["reg_media"].mean(axis=0)
    helado = (r["hielo_espesor"] > 0).mean(axis=0) >= 0.5
    print(f"Aire a 2 m, media anual por celda -- minima: {T_media_anual.min():.1f} C | maxima: {T_media_anual.max():.1f} C")

    fig, ax = plt.subplots(figsize=(14, 7))
    im = ax.imshow(T_media_anual, extent=[-180, 180, -90, 90], origin="upper", cmap="RdBu_r",
                   aspect="auto", vmin=-50, vmax=50, interpolation="nearest")
    fig.colorbar(im, ax=ax, label="Aire a 2 m, media anual (C); escala fija de -50 a +50 C")
    # costa (limite tierra/agua) y hielo marino
    tierra = (tipo_superficie == TIERRA).astype(float)
    ax.contour(tierra, levels=[0.5], colors="black", linewidths=0.8, extent=[-180, 180, -90, 90], origin="upper")
    if helado.any():
        ax.contourf(helado.astype(float), levels=[0.5, 1.5], colors="none", hatches=["//"],
                    extent=[-180, 180, -90, 90], origin="upper")
    ax.set_xlabel("Longitud (grados)")
    ax.set_ylabel("Latitud (grados)")
    ax.set_title(f"Temperatura media anual del aire a 2 m -- {nombre_mapa} (72x36). Rayado: hielo marino al menos medio año")
    fig.tight_layout()

    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)
    ruta = os.path.join(CARPETA_RESULTADOS, f"mapa_calor_{nombre_mapa}_{date.today().isoformat()}.png")
    fig.savefig(ruta, dpi=150)
    print(f"Mapa guardado en {ruta}")

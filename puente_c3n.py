import os
import json

import numpy as np

from fase1_geografia import TIERRA, AGUA
from rejilla import FILAS, COLUMNAS

RUTA_PUENTE = os.path.expanduser(os.path.join("~", "Documentos", "B3N", "mapa_activo", "mapa_activo_m3n.json"))


def cargar_mapa_activo_de_c3n():
    if not os.path.exists(RUTA_PUENTE):
        raise FileNotFoundError(
            f"No se encontro el mapa puente en {RUTA_PUENTE}. "
            "Activa un mapa en C3N primero (la exportacion es automatica al activar)."
        )
    with open(RUTA_PUENTE, "r", encoding="utf-8") as f:
        datos = json.load(f)

    filas = datos["filas"]
    columnas = datos["columnas"]
    if filas != FILAS or columnas != COLUMNAS:
        raise ValueError(
            f"El mapa puente tiene resolucion {filas}x{columnas}, "
            f"pero M3N espera {FILAS}x{COLUMNAS}."
        )

    tipo_superficie = np.array(datos["agua_tierra"], dtype=int).reshape(FILAS, COLUMNAS)
    altitud_metros = np.array(datos["altitud_metros"], dtype=float).reshape(FILAS, COLUMNAS)

    valores_validos = {TIERRA, AGUA}
    if not set(np.unique(tipo_superficie)).issubset(valores_validos):
        raise ValueError(f"El mapa puente contiene valores de agua_tierra fuera de {valores_validos}")

    nombre_mapa = datos.get("nombre_mapa_origen", "desconocido")
    return tipo_superficie, altitud_metros, nombre_mapa


if __name__ == "__main__":
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()
    num_tierra = int(np.sum(tipo_superficie == TIERRA))
    num_agua = int(np.sum(tipo_superficie == AGUA))
    print(f"Mapa cargado: '{nombre_mapa}'")
    print(f"Celdas de tierra: {num_tierra} | Celdas de agua: {num_agua} (de {FILAS * COLUMNAS} totales)")
    print(f"Altitud minima: {altitud_metros.min():.1f} m | Altitud maxima: {altitud_metros.max():.1f} m | "
          f"Altitud media: {altitud_metros.mean():.1f} m")

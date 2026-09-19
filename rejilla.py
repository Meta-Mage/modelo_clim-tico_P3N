# rejilla.py — Fase 0: rejilla de trabajo lat/lon para el modelo vectorizado.
#
# Resolución de trabajo: 72 columnas x 36 filas (5° por celda). Usa la MISMA
# convención de "centro de celda" que ya usa C3N en mapa.js:
#   lat = 90 - (fila + 0.5) * grados_por_fila
#   lon = (col + 0.5) * grados_por_col - 180
# para que el futuro puente C3N -> M3N no tenga que traducir nada entre los
# dos sistemas de coordenadas.

import numpy as np

FILAS = 36
COLUMNAS = 72

GRADOS_POR_FILA = 180 / FILAS    # 5°
GRADOS_POR_COL = 360 / COLUMNAS  # 5°


def lat_de_fila(fila):
    return 90 - (fila + 0.5) * GRADOS_POR_FILA


def lon_de_columna(columna):
    return (columna + 0.5) * GRADOS_POR_COL - 180


LATITUDES_GRADOS = lat_de_fila(np.arange(FILAS))
LONGITUDES_GRADOS = lon_de_columna(np.arange(COLUMNAS))


if __name__ == "__main__":
    print(f"Rejilla: {FILAS} filas x {COLUMNAS} columnas ({GRADOS_POR_FILA}° por celda)")
    print(f"Latitudes: de {LATITUDES_GRADOS[0]:.1f}° a {LATITUDES_GRADOS[-1]:.1f}°")
    print(f"Longitudes: de {LONGITUDES_GRADOS[0]:.1f}° a {LONGITUDES_GRADOS[-1]:.1f}°")
    print(f"Fila 17 (justo al norte del ecuador): {LATITUDES_GRADOS[17]:.1f}°")
    print(f"Fila 18 (justo al sur del ecuador):    {LATITUDES_GRADOS[18]:.1f}°")

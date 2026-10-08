# modo_tierra.py -- v3.0: mapa de la Tierra en la rejilla de M3N (72 x 36, 5 grados)
#
# Solo para el MODO TIERRA (calibracion y validacion; ver parametros.py,
# variable de entorno M3N_MODO=tierra). Nunca se usa para P3N.
#
# Tierra/agua: paquete global-land-mask (mascara GLOBE de la NOAA a ~1 km),
# muestreada con 20 x 20 puntos por celda; la celda es tierra si al menos
# la mitad de sus puntos lo es (mismo criterio que el mapa de P3N).
# Altitud: la misma simplificacion que la calibracion de la v2.2.2, para
# que los resultados sean comparables: Antartida 2300 m y Groenlandia
# 2000 m (altitudes medias aproximadas de sus casquetes); el resto de la
# tierra a 0 m. LIMITACION conocida: las demas montañas (Himalaya, Andes,
# Rocosas) no estan; en la v3.0 cambian la presion en superficie, asi que
# su ausencia hace la Tierra simulada algo mas calida sobre los
# continentes. Mejora pendiente: altitud media real por celda (ETOPO).
#
# Requiere:  pip install global-land-mask

import numpy as np
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS, GRADOS_POR_FILA, GRADOS_POR_COL
from fase1_geografia import TIERRA, AGUA

ALTITUD_ANTARTIDA = 2300.0
ALTITUD_GROENLANDIA = 2000.0
MUESTRAS_POR_LADO = 20


def mapa_tierra():
    """Devuelve (tipo_superficie, altitud_metros), arrays (FILAS, COLUMNAS)."""
    from global_land_mask import globe
    n = MUESTRAS_POR_LADO
    sub = (np.arange(n) + 0.5) / n
    lat = 90 - (np.arange(FILAS)[:, None] + sub[None, :]).ravel() * GRADOS_POR_FILA
    lon = (np.arange(COLUMNAS)[:, None] + sub[None, :]).ravel() * GRADOS_POR_COL - 180
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    es = globe.is_land(np.clip(LAT, -89.999, 89.999), np.clip(LON, -179.999, 179.999))
    fraccion = es.reshape(FILAS, n, COLUMNAS, n).mean(axis=(1, 3))
    tierra = fraccion >= 0.5
    tipo = np.where(tierra, TIERRA, AGUA)
    lat_c = LATITUDES_GRADOS[:, None] * np.ones((1, COLUMNAS))
    lon_c = np.ones((FILAS, 1)) * LONGITUDES_GRADOS[None, :]
    altitud = np.zeros((FILAS, COLUMNAS))
    altitud[tierra & (lat_c < -60)] = ALTITUD_ANTARTIDA
    # v3.1: la caja de Groenlandia ya no incluye Baffin ni Ellesmere (revision del 05/10/2026):
    # al sur de 75 N, Groenlandia esta al este de -60; al norte, al este de -70 (Thule, -68,7)
    groenlandia = tierra & (lat_c > 59) & (lon_c < -10) & np.where(lat_c < 75, lon_c > -60, lon_c > -70)
    altitud[groenlandia] = ALTITUD_GROENLANDIA
    return tipo, altitud


if __name__ == "__main__":
    tipo, alt = mapa_tierra()
    peso = np.cos(np.radians(LATITUDES_GRADOS))[:, None] * np.ones((1, COLUMNAS))
    print(f"Tierra: {100 * (peso * (tipo == TIERRA)).sum() / peso.sum():.1f} % de la superficie (real: 29,2 %)")
    for f in range(FILAS):
        print("".join("A" if alt[f, c] > 1500 else ("#" if tipo[f, c] == TIERRA else ".") for c in range(COLUMNAS)))

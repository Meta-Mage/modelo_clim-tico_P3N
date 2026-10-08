# fase6_equilibrio.py -- v3.10.0: criterio de equilibrio del modelo con nucleo dinamico (I16).
# Diseno: DISENO_FASE_6.3.md §6.3 (decisiones de Carlos del 06/10) y §6.10.
#
# Con el nucleo dinamico hay tiempo meteorologico: el criterio de la v3.1 (que ninguna celda cambie su media
# anual mas de 0,015 K) no se cumpliria nunca. El equilibrio se juzga con MEDIAS GLOBALES ANUALES en una
# ventana de los ultimos 5 años, con cuatro condiciones a la vez (umbrales de partida 🔶, se fijan en firme
# midiendo el ruido de la primera simulacion en modo Tierra):
#   - balance en el tope de la atmosfera: |N medio| < 0,2 W/m2;
#   - tendencia del aire a 2 m: |pendiente| < 0,02 K/año;
#   - tendencia del area de hielo marino: |pendiente| < 0,1 % del oceano por año;
#   - tendencia del agua del suelo: |pendiente| < 1 % de su capacidad por año.
# La pendiente es la de minimos cuadrados de las medias anuales frente al año.

import numpy as np

VENTANA = 5
UMBRALES = {"N": 0.2, "T2m": 0.02, "hielo": 0.001, "agua_suelo": 0.01}
CLAVES = ("N", "T2m", "hielo", "agua_suelo")


def pendiente(y):
    """Pendiente de minimos cuadrados de y frente a 0, 1, ..., n-1 (unidades de y por año)."""
    y = np.asarray(y, dtype=float)
    x = np.arange(len(y), dtype=float)
    x = x - x.mean()
    return float((x * (y - y.mean())).sum() / (x * x).sum())


def evaluar(serie, ventana=VENTANA, umbrales=None):
    """serie: dict con listas de medias anuales (una por año simulado) de las claves de CLAVES; una lista
    vacia o None (p. ej., sin tierra o sin oceano) no cuenta. Devuelve (en_equilibrio, valores), con valores
    = {clave: (valor, umbral, cumple)} para las ultimas `ventana` medias. Con menos años, (False, {})."""
    u = dict(UMBRALES if umbrales is None else umbrales)
    n = len(serie["N"])
    if n < ventana:
        return False, {}
    valores = {}
    m = float(np.mean(serie["N"][-ventana:]))
    valores["N"] = (m, u["N"], abs(m) < u["N"])
    for k in ("T2m", "hielo", "agua_suelo"):
        y = serie.get(k)
        if y is None or len(y) == 0 or any(v is None for v in y[-ventana:]):
            continue
        p = pendiente(y[-ventana:])
        valores[k] = (p, u[k], abs(p) < u[k])
    return all(v[2] for v in valores.values()), valores

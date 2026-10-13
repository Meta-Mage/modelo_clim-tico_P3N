# fase52_hitran.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §7): lector PROPIO de las lineas de HITRAN
# (formato .par de 160 caracteres, el de los .data que guarda herramientas/descargar_hitran.py) y reparto del espectro
# en BANDAS. Asi M3N no depende de hitran-api para leer los datos (decision T3: datos, no librerias); la prueba lo
# compara linea a linea con el lector de hitran-api.
#
# Columnas del formato .par (HITRAN 2004 en adelante; Rothman et al. 2005, JQSRT 96, 139, tabla 1 ⚠️):
#   molecula [0:2], isotopologo [2:3] (1..9; '0' = 10; 'A' = 11, 'B' = 12...), nu [3:15] (cm^-1),
#   S [15:25] (cm^-1/(molecula cm^-2), a 296 K y con la abundancia terrestre del isotopologo), A de Einstein [25:35],
#   gamma_air [35:40], gamma_self [40:45] (cm^-1/atm, HWHM a 296 K), E'' [45:55] (cm^-1), n_air [55:59],
#   delta_air [59:67] (cm^-1/atm).
# Las intensidades llevan la abundancia de los isotopos de la Tierra: P3N se supone con la misma composicion isotopica.

import re

import numpy as np

COLUMNAS = {"molec": (0, 2), "iso": (2, 3), "nu0": (3, 15), "S": (15, 25), "A": (25, 35), "g_air": (35, 40),
            "g_self": (40, 45), "E": (45, 55), "n_air": (55, 59), "d_air": (59, 67)}

# Bandas (cm^-1), decision tecnica T8 (PROVISIONAL): los limites de RRTMG (Mlawer et al. 1997, J. Geophys. Res. 102,
# 16663, infrarrojo, 16 bandas; Iacono et al. 2008, luz, 14 bandas). Comprobados el 12/10/2026 con los del codigo de
# RRTMG en CESM (radconstants.F90: wavenumber1/2_longwave y wavenum_low/high) ✅; la primera version ponia 2380 en
# lugar de 2390 entre las bandas 14 y 15 del infrarrojo. Son limites elegidos por donde absorbe cada gas, publicados
# y muy usados; se adoptan los LIMITES, no las tablas de RRTMG (que se calculan aqui con HITRAN y para la atmosfera de P3N). Revisables al medir el error (§4).
BANDAS_IR = ((10, 350), (350, 500), (500, 630), (630, 700), (700, 820), (820, 980), (980, 1080), (1080, 1180),
             (1180, 1390), (1390, 1480), (1480, 1800), (1800, 2080), (2080, 2250), (2250, 2390), (2390, 2600),
             (2600, 3250))
BANDAS_LUZ = ((2600, 3250), (3250, 4000), (4000, 4650), (4650, 5150), (5150, 6150), (6150, 7700), (7700, 8050),
              (8050, 12850), (12850, 16000), (16000, 22650), (22650, 29000), (29000, 38000), (38000, 50000),
              (820, 2600))


def _iso(c):
    c = c.strip()
    if c == "0":
        return 10
    if c.isdigit():
        return int(c)
    return 11 + ord(c) - ord("A")


def _flotante(s):
    try:
        return float(s)
    except ValueError:
        # formato raro de algunas intensidades muy pequeñas, p. ej. "2.700-164" (sin la E): lo trata igual hitran-api
        m = re.search(r"(\d\.\d+)\-(\d+)", s)
        if m:
            return float(m.group(1) + "E-" + m.group(2))
        raise


def leer_par(ruta, nu_min=None, nu_max=None):
    """Lee un archivo .par/.data de HITRAN. Devuelve un dict de arrays (las claves de COLUMNAS), solo con las lineas
    con nu0 en [nu_min, nu_max] si se dan, en el orden del archivo (HITRAN las da ordenadas por nu)."""
    with open(ruta, "r", encoding="ascii") as f:
        lineas = [l for l in f.read().split("\n") if len(l) >= 67]
    nu = np.array([float(l[3:15]) for l in lineas])
    sel = np.ones(len(nu), dtype=bool)
    if nu_min is not None:
        sel &= nu >= nu_min
    if nu_max is not None:
        sel &= nu <= nu_max
    elegidas = [l for l, s in zip(lineas, sel) if s]
    out = {"nu0": nu[sel]}
    for k, (a, b) in COLUMNAS.items():
        if k == "nu0":
            continue
        if k == "molec":
            out[k] = np.array([int(l[a:b]) for l in elegidas], dtype=int)
        elif k == "iso":
            out[k] = np.array([_iso(l[a:b]) for l in elegidas], dtype=int)
        else:
            out[k] = np.array([_flotante(l[a:b]) for l in elegidas], dtype=float)
    return out


def lineas_de_isotopologo(lineas, iso):
    """Las lineas de un isotopologo (con las claves que usa fase52_lineas.coeficiente_absorcion)."""
    m = lineas["iso"] == iso
    return {k: v[m] for k, v in lineas.items()}

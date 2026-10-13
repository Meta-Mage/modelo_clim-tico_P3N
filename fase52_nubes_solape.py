# fase52_nubes_solape.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §3.4, decision tecnica T23): como se
# SUPERPONEN las nubes de capas distintas y como entran en la radiacion por bandas (McICA). No se usa todavia en el
# modelo.
#
# El problema: una celda de M3N mide 400-550 km; en cada capa solo una fraccion f esta nublada. Para la radiacion
# importa como se apilan esas fracciones: si las nubes de capas vecinas estan unas encima de otras (solapamiento
# "maximo") cubren menos cielo que si estan repartidas al azar ("aleatorio").
#
# T23 (provisional):
#   - SOLAPAMIENTO EXPONENCIAL-ALEATORIO (Hogan e Illingworth 2000, Q. J. R. Meteorol. Soc. 126, 2903 ✅t): entre dos
#     capas nubladas seguidas, la cobertura conjunta es  C = alfa max(f1, f2) + (1 - alfa) (f1 + f2 - f1 f2),
#     con alfa = exp(-dz / L) (dz, distancia entre los centros de las capas; L, "longitud de decorrelacion"); las nubes
#     separadas por una capa despejada se superponen al azar. Formulas como las implementa ecRad
#     (radiation_cloud_cover.F90, cum_cloud_cover_exp_ran; Apache 2.0, leido el 13/10/2026 ✅).
#   - L segun la latitud, la forma suavizada de Shonk et al. 2010 (Q. J. R. Meteorol. Soc. 136, 1191 ⚠️) que usa el
#     IFS: L = 0,75 + 2,149 cos^2(latitud) km (ecRad, ifs/cloud_overlap_decorr_len.F90 ✅). Es una CALIBRACION TERRESTRE
#     (c): se mantiene en el modo Tierra y en P3N mientras no haya nada mejor (limitacion declarada).
#   - McICA (Pincus, Barker y Morcrette 2003, J. Geophys. Res. 108(D13) ✅t): cada punto g de cada banda ve UNA
#     subcolumna (nublada o despejada en cada capa) sacada al azar con esas reglas; en media es el calculo exacto, con
#     ruido aleatorio sin sesgo. El generador de subcolumnas es el de ecRad (radiation_cloud_generator.F90,
#     generate_column_exp_ran ✅): primero se sortea la capa mas alta nublada segun la cobertura acumulada desde
#     arriba, y despues, capa a capa hacia abajo, si la nube sigue o empieza, con las probabilidades condicionadas
#     exactas. Aqui sin heterogeneidad del agua dentro de la nube (cada capa nublada, con su agua media): limitacion
#     declarada ⚠️ (ecRad la añade con una distribucion gamma o lognormal).

import numpy as np

F_MAX = 1.0 - 10.0 * np.finfo(float).eps


def longitud_decorrelacion_km(latitud_grados):
    """L (km) de Shonk et al. 2010 suavizada en el ecuador (la del IFS)."""
    c = np.cos(np.radians(latitud_grados))
    return 0.75 + 2.149 * c * c


def alfa_capas(z_centros_km, L_km):
    """Parametro de solapamiento entre cada par de capas seguidas: exp(-|dz| / L) (n-1 valores)."""
    z = np.asarray(z_centros_km, dtype=float)
    return np.exp(-np.abs(np.diff(z)) / L_km)


def cobertura_acumulada(f, alfa):
    """(acumulada, par): cobertura vista desde arriba hasta la base de cada capa, y cobertura de cada par de capas
    seguidas (exponencial-aleatorio; con alfa = 1, maximo-aleatorio). f de arriba abajo."""
    f = np.asarray(f, dtype=float)
    alfa = np.asarray(alfa, dtype=float)
    n = len(f)
    par = alfa * np.maximum(f[:-1], f[1:]) + (1.0 - alfa) * (f[:-1] + f[1:] - f[:-1] * f[1:])
    acum = np.empty(n)
    acum[0] = f[0]
    prod = 1.0 - f[0]
    for k in range(n - 1):
        if f[k] >= F_MAX:
            prod = 0.0
        else:
            prod = prod * (1.0 - par[k]) / (1.0 - f[k])
        acum[k + 1] = 1.0 - prod
    return acum, par


def cobertura_total(f, alfa):
    return float(cobertura_acumulada(f, alfa)[0][-1])


def subcolumnas(f, alfa, n, rng):
    """(n, capas) de booleanos: n subcolumnas nubladas o no en cada capa, con el solapamiento de arriba (el algoritmo de
    ecRad, generate_column_exp_ran, sin heterogeneidad). rng: numpy.random.Generator."""
    f = np.asarray(f, dtype=float)
    nl = len(f)
    acum, par = cobertura_acumulada(f, alfa)
    total = acum[-1]
    out = np.zeros((n, nl), dtype=bool)
    if total <= 1e-9:
        return out
    nubladas = np.nonzero(f > 0)[0]
    ini, fin = nubladas[0], nubladas[-1]
    sobre = np.diff(acum)                      # lo que aumenta la cobertura acumulada en la capa siguiente
    nubosas = rng.random(n) < total            # que subcolumnas tienen alguna nube
    for i in np.nonzero(nubosas)[0]:
        disparo = rng.random() * total
        k = ini
        while disparo > acum[k] and k < fin:
            k += 1
        out[i, k] = True
        for j in range(k + 1, fin + 1):
            u = rng.random()
            if out[i, j - 1]:
                # hay nube encima: sigue con probabilidad (f_j + f_{j-1} - par) / f_{j-1}
                if u * f[j - 1] < f[j] + f[j - 1] - par[j - 1]:
                    out[i, j] = True
            else:
                # despejado encima: nube nueva con probabilidad (par - sobre - f_{j-1}) / (acum_{j-1} - f_{j-1})
                den = acum[j - 1] - f[j - 1]
                if den > 0 and u * den < par[j - 1] - sobre[j - 1] - f[j - 1]:
                    out[i, j] = True
    return out

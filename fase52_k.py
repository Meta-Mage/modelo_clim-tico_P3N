# fase52_k.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §3.1 y §7, paso 1 de "lo siguiente"): de los
# coeficientes de absorcion linea a linea (fase52_lineas.py) a TABLAS k por banda, el metodo de la distribucion k
# correlacionada de la opcion (B). No se usa todavia en el modelo.
#
# Idea (Lacis y Oinas 1991, J. Geophys. Res. 96, 9027 ⚠️; Fu y Liou 1992, J. Atmos. Sci. 49, 2139 ⚠️; Goody et al.
# 1989, J. Quant. Spectrosc. Radiat. Transf. 42, 539 ⚠️): dentro de una banda, la transmision media
#     T(u) = (1/Δnu) ∫ exp(-k(nu) u) dnu
# solo depende de CUANTAS veces aparece cada valor de k, no de DONDE. Ordenando k(nu) de menor a mayor se obtiene una
# funcion suave y monotona k(g), g en [0, 1] (la fraccion de la banda con un coeficiente menor), y
#     T(u) = ∫_0^1 exp(-k(g) u) dg  ≈  Σ_i w_i exp(-k(g_i) u)
# con pocos puntos de cuadratura en g (decenas) en vez de cientos de miles de numeros de onda. En UNA capa homogenea es
# exacto (solo es reordenar la integral); el error viene de la cuadratura. Con VARIAS capas se supone que el orden de
# los numeros de onda es el mismo en todas (de ahi "correlacionada"): exacto si k(nu) de cada capa es funcion monotona
# de la de otra (p. ej., una sola linea a distinta presion: lo comprueba la prueba), aproximado si no (lineas que
# se cruzan con la altura, gases que se solapan: el error se medira contra el calculo linea a linea, §4).
#
# Cuadratura en g (decision tecnica T6, PROVISIONAL): Gauss-Legendre en tramos de g, con los tramos cada vez mas
# estrechos hacia g = 1, donde estan los centros de las lineas (pocos numeros de onda con k muy grande que dominan las
# capas altas y los caminos cortos). Es la practica habitual (los modelos operativos reparten 16 puntos asi ⚠️); el
# reparto exacto de M3N se fijara midiendo el error de los flujos frente al calculo linea a linea (§4), no copiando.

import numpy as np

# tramos de g y puntos de Gauss en cada uno: 16 puntos (T6, provisional)
TRAMOS_G = ((0.0, 0.3, 2), (0.3, 0.6, 2), (0.6, 0.85, 3), (0.85, 0.95, 3), (0.95, 0.99, 3), (0.99, 1.0, 3))


def nodos_g(tramos=TRAMOS_G):
    """Puntos y pesos de la cuadratura en g: Gauss-Legendre en cada tramo (a, b, n). Los pesos suman 1."""
    g, w = [], []
    for a, b, n in tramos:
        x, p = np.polynomial.legendre.leggauss(int(n))
        g.append(0.5 * (b - a) * x + 0.5 * (a + b))
        w.append(0.5 * (b - a) * p)
    return np.concatenate(g), np.concatenate(w)


def distribucion_k(k_nu, dnu=None):
    """k(nu) de una banda -> (k ordenado de menor a mayor, g acumulada al final de cada valor). dnu: anchura de cada
    numero de onda (por defecto, todos iguales); g llega a 1."""
    k_nu = np.asarray(k_nu, dtype=float)
    w = np.ones_like(k_nu) if dnu is None else np.broadcast_to(np.asarray(dnu, dtype=float), k_nu.shape)
    orden = np.argsort(k_nu, kind="stable")
    G = np.cumsum(w[orden])
    return k_nu[orden], G / G[-1]


def k_en_g(k_ord, G, g):
    """La inversa de la distribucion: el menor k cuya fraccion acumulada alcanza g (exacta para la rejilla dada)."""
    i = np.searchsorted(G, np.asarray(g, dtype=float), side="left")
    return k_ord[np.clip(i, 0, len(k_ord) - 1)]


def tabla_k(k_nu, dnu=None, tramos=TRAMOS_G):
    """La tabla k de una banda: (g, pesos, k(g))."""
    g, w = nodos_g(tramos)
    k_ord, G = distribucion_k(k_nu, dnu)
    return g, w, k_en_g(k_ord, G, g)


def tabla_k_capas(k_nu_capas, dnu=None, tramos=TRAMOS_G):
    """Tablas k de varias capas con la hipotesis correlacionada: cada capa se ordena por separado y se evalua en los
    MISMOS puntos g. k_nu_capas: (capas, nu). Devuelve (g, pesos, k (capas, puntos g))."""
    g, w = nodos_g(tramos)
    k = np.array([k_en_g(*distribucion_k(kc, dnu), g) for kc in np.atleast_2d(k_nu_capas)])
    return g, w, k


def transmision_lbl(k_nu, u, dnu=None):
    """Transmision media de la banda linea a linea: k_nu (capas, nu) o (nu,), u (capas,) o escalar."""
    k_nu = np.atleast_2d(np.asarray(k_nu, dtype=float))
    u = np.atleast_1d(np.asarray(u, dtype=float))
    tau = (k_nu * u[:, None]).sum(axis=0)
    w = np.ones(k_nu.shape[1]) if dnu is None else np.broadcast_to(np.asarray(dnu, dtype=float), (k_nu.shape[1],))
    return float((np.exp(-tau) * w).sum() / w.sum())


def transmision_ck(pesos, k_g, u):
    """La misma transmision con la tabla k: k_g (capas, puntos g) o (puntos g,), u (capas,) o escalar."""
    k_g = np.atleast_2d(np.asarray(k_g, dtype=float))
    u = np.atleast_1d(np.asarray(u, dtype=float))
    return float((np.asarray(pesos) * np.exp(-(k_g * u[:, None]).sum(axis=0))).sum())

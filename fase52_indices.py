# fase52_indices.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §7, decision tecnica T17): los INDICES DE
# REFRACCION complejos m = n + i k del agua liquida y del hielo, a cualquier longitud de onda (y, el agua, a cualquier
# temperatura), para la teoria de Mie de las nubes (fase52_mie.py). No se usa todavia en el modelo.
#
# Datos: datos/indices_refraccion_agua_hielo.npz (de dominio publico; lo prepara
# herramientas/preparar_indices_refraccion.py, donde estan las fuentes y sus huellas).
#
# T17 (provisional):
#   - AGUA LIQUIDA. Rowe, Fergoda y Neshyba 2020 (240, 253, 263 y 273 K, de 0,667 µm en adelante), porque en el
#     infrarrojo el indice cambia mucho con la temperatura y las nubes de agua subfundida son comunes (medido con los
#     datos: k a 12 µm es 0,31 a 240 K y 0,20 a 298 K). Como nodo de 298 K, Hale y Querry 1973 (25 °C), que es el dato
#     clasico "a temperatura ambiente" del que parten esos trabajos; asi la temperatura se interpola entre 240 y
#     298 K. Por debajo de 0,667 µm (donde el agua apenas absorbe), Hale y Querry, sin temperatura.
#     Mas alla de 200 µm (fin de Hale y Querry), solo Rowe et al.: por encima de 273 K se usa la de 273 K.
#     Fuera de 240-298 K, el valor del extremo (el agua liquida no existe por debajo de ~235 K).
#   - HIELO. Warren y Brandt 2008, a 266 K, sin dependencia con la temperatura (limitacion declarada ⚠️).
#   - Interpolacion: lineal en log(lambda) para n y para log(k) (k cambia varios ordenes de magnitud), y lineal en T.

import os

import numpy as np

RUTA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "datos", "indices_refraccion_agua_hielo.npz")
T_HQ = 298.15                      # 25 °C
_DATOS = None


def datos():
    global _DATOS
    if _DATOS is None:
        with np.load(RUTA) as d:
            _DATOS = {k: d[k] for k in d.files}
    return _DATOS


def _interp_log(lam, tabla_lam, n, k):
    """n y k en las longitudes de onda lam (µm), lineales en log(lambda); k en escala logaritmica."""
    x = np.log(lam)
    xt = np.log(tabla_lam)
    return np.interp(x, xt, n), np.exp(np.interp(x, xt, np.log(k)))


def _comprobar_rango(lam, lo, hi, que):
    if np.any(lam < lo * (1 - 1e-12)) or np.any(lam > hi * (1 + 1e-12)):
        raise ValueError(f"{que}: longitud de onda fuera de los datos ({lo:g}-{hi:g} µm)")


def indice_hielo(lam_um):
    """m = n + i k del hielo (Warren y Brandt 2008, 266 K) en lam_um (µm)."""
    d = datos()
    lam = np.asarray(lam_um, dtype=float)
    _comprobar_rango(lam, d["wb_lam"][0], d["wb_lam"][-1], "hielo")
    n, k = _interp_log(lam, d["wb_lam"], d["wb_n"], d["wb_k"])
    return n + 1j * k


def indice_agua(lam_um, T):
    """m = n + i k del agua liquida en lam_um (µm) a la temperatura T (K) (ver T17 arriba)."""
    d = datos()
    lam = np.atleast_1d(np.asarray(lam_um, dtype=float))
    _comprobar_rango(lam, d["hq_lam"][0], d["rowe_lam"][-1], "agua")
    T = float(T)
    n = np.empty_like(lam)
    k = np.empty_like(lam)
    corto = lam < d["rowe_lam"][0]                             # solo Hale y Querry
    medio = ~corto & (lam <= d["hq_lam"][-1])                   # Rowe (240-273 K) y Hale y Querry (298 K)
    largo = lam > d["hq_lam"][-1]                               # solo Rowe
    if corto.any():
        n[corto], k[corto] = _interp_log(lam[corto], d["hq_lam"], d["hq_n"], d["hq_k"])
    if medio.any() or largo.any():
        nodos_T = list(d["rowe_T"])
        n_nodos, k_nodos = [], []
        sub = lam[medio | largo]
        for i in range(len(nodos_T)):
            a, b = _interp_log(sub, d["rowe_lam"], d["rowe_n"][i], d["rowe_k"][i])
            n_nodos.append(a)
            k_nodos.append(np.log(b))
        n_nodos, k_nodos = np.array(n_nodos), np.array(k_nodos)
        m_sub = medio[medio | largo]
        if m_sub.any():                                         # nodo de 298 K solo donde hay Hale y Querry
            a, b = _interp_log(sub[m_sub], d["hq_lam"], d["hq_n"], d["hq_k"])
            nodos = np.array(nodos_T + [T_HQ])
            Tc = min(max(T, nodos[0]), nodos[-1])
            nn = np.vstack([n_nodos[:, m_sub], a])
            kk = np.vstack([k_nodos[:, m_sub], np.log(b)])
            n[medio] = [np.interp(Tc, nodos, nn[:, j]) for j in range(nn.shape[1])]
            k[medio] = np.exp([np.interp(Tc, nodos, kk[:, j]) for j in range(kk.shape[1])])
        l_sub = ~m_sub
        if l_sub.any():
            nodos = np.array(nodos_T)
            Tc = min(max(T, nodos[0]), nodos[-1])
            n[largo] = [np.interp(Tc, nodos, n_nodos[:, l_sub][:, j]) for j in range(l_sub.sum())]
            k[largo] = np.exp([np.interp(Tc, nodos, k_nodos[:, l_sub][:, j]) for j in range(l_sub.sum())])
    m = n + 1j * k
    return m if np.ndim(lam_um) else m[0]

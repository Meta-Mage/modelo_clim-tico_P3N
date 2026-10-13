# fase52_referencia.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §7): COLUMNAS DE REFERENCIA de P3N y el
# calculo LINEA A LINEA de una columna (espesor optico de cada capa y numero de onda, de todos los gases), con las
# mismas funciones que el generador de tablas (fase52_generador: rejilla adaptada T12, nucleo y alas T13, pedestal T4,
# continuo MT_CKD). Lo usan herramientas/medir_solapamiento.py (decision T9) y servira de referencia en la validacion
# (§4). No se usa en el modelo.
#
# Columnas de prueba (solo para medir; NO son el clima de P3N):
#   - 30 capas con la receta de niveles de M3N (z = 40 km (k/N)^1,6, sigma = e^(-z/7,5 km); tope ~4,8 hPa) y la gravedad
#     de P3N;
#   - temperatura: la del suelo (288 K "ref", 300 K "tropical", 255 K "polar") bajando 6,5 K/km hasta 217 K y
#     constante por encima (forma de la atmosfera estandar);
#   - vapor: humedad relativa de Manabe y Wetherald (1967) ⚠️, 0,77 (p/ps - 0,02)/0,98, con un minimo de 3 ppm;
#   - CO2 174,7 ppm, CH4 358 ppb, N2O 201,9 ppb (FICHA_COMPOSICION_P3N.md §6.6);
#   - ozono: una campana en log p con maximo de 5 ppm a 10 hPa (forma aproximada; la cantidad de P3N esta por decidir).

import json
import os

import numpy as np

import fase52_continuo as C
import fase52_generador as G
import fase52_hitran as H

GASES = {"H2O": 1, "CO2": 2, "O3": 3, "N2O": 4, "CH4": 6}
VMR = {"CO2": 174.7e-6, "CH4": 358e-9, "N2O": 201.9e-9}
M_AIRE = 28.964e-3 / 6.02214076e23          # kg por molecula de aire seco
CP = 1004.64
G_P3N = 9.05                                # ~parametros.P3N_GRAVEDAD (solo para la columna de prueba)
UMBRAL_LBL = 1e-9                           # descarte de lineas en el linea a linea (mas estricto que T14)


COLUMNAS = {"ref": 288.0, "tropical": 300.0, "polar": 255.0}       # temperatura del suelo de cada columna de prueba


def columna(n=30, ps=1000.0, Ts=288.0):
    zk = 40.0 * (np.arange(n, -1, -1) / n) ** 1.6                    # la receta de niveles de M3N (DISENO_FASE_2.3 §1)
    ph = ps * np.exp(-zk / 7.5)                                       # seminiveles, de arriba abajo (hPa)
    ph[0] = ps * np.exp(-40.0 / 7.5)
    pm = np.sqrt(ph[1:] * ph[:-1])
    z = lambda p: -7.0 * np.log(p / ps)                               # km aproximados (solo para el perfil)
    Tf = lambda p: np.maximum(Ts - 6.5 * z(p), 217.0)
    Th, Tm = Tf(ph), Tf(pm)
    es = 6.112 * np.exp(17.67 * (Tm - 273.15) / (Tm - 29.65))        # hPa (Bolton 1980 ⚠️)
    rh = np.clip(0.77 * (pm / ps - 0.02) / 0.98, 0.0, 1.0)
    vmr = {"H2O": np.maximum(rh * es / pm, 3e-6)}
    for g, v in VMR.items():
        vmr[g] = np.full(n, v)
    vmr["O3"] = 5e-6 * np.exp(-0.5 * (np.log(pm / 10.0) / 1.2) ** 2)
    n_aire = np.diff(ph) * 100.0 / (G_P3N * M_AIRE) / 1e4             # moleculas de aire por cm^2 en cada capa
    return {"ph": ph, "pm": pm, "Th": Th, "Tm": Tm, "Ts": Th[-1], "vmr": vmr, "n_aire": n_aire,
            "dp": np.diff(ph) * 100.0}


def cargar(hitran, gas, a, b):
    with open(os.path.join(hitran, "tips.json"), encoding="utf-8") as f:
        tips = json.load(f)
    with open(os.path.join(hitran, "isotopologos.json"), encoding="utf-8") as f:
        iso = json.load(f)
    m = GASES[gas]
    isotopos = {}
    for clave, inf in iso.items():
        mm, ii = (int(x) for x in clave.split("_"))
        if mm == m and tips["Q"].get(clave):
            isotopos[ii] = {"masa_u": inf["masa_u"], "T": np.array(tips["T_K"], float), "Q": np.array(tips["Q"][clave])}
    ruta = os.path.join(hitran, gas + ".data")
    if not os.path.exists(ruta):
        return None, isotopos
    return H.leer_par(ruta, a - 25.0, b + 25.0), isotopos


def k_linea_a_linea(lineas, isotopos, nu, p_hPa, T, vmr, u_col, es_vapor, mtckd):
    """k(nu) (cm^2/molecula) de un gas en una capa, sobre la rejilla comun nu."""
    p = p_hPa / 1013.25
    k = np.zeros_like(nu)
    for iso, inf in isotopos.items():
        m = lineas["iso"] == iso
        if not m.any():
            continue
        sub = {kk: v[m] for kk, v in lineas.items()}
        Q_T, Q_ref = G.q_de(inf["T"], inf["Q"], T), G.q_de(inf["T"], inf["Q"], G.T_REF)
        S, c, gL, aD = G.parametros(sub, T, p, vmr * p, Q_T, Q_ref, inf["masa_u"])
        q, _ = G.descartar(S, u_col, UMBRAL_LBL)
        if len(q):
            k += G.k_rapido(nu, S[q], c[q], gL[q], aD[q], quitar_pedestal=es_vapor)
    if es_vapor:
        k = np.maximum(k, 0.0)
        if mtckd is not None:
            k += C.k_continuo(nu, T, p_hPa, vmr, mtckd)
    return k


def rejilla_comun(datos, banda, col):
    """T12 con las anchuras de la capa de ARRIBA (las mas estrechas) de todos los gases: resuelve todas las capas."""
    a, b = banda
    cs, gls, ads = [], [], []
    p, T = col["pm"][0] / 1013.25, col["Tm"][0]
    for gas, (lineas, isotopos) in datos.items():
        if lineas is None:
            continue
        for iso, inf in isotopos.items():
            m = (lineas["iso"] == iso) & (lineas["nu0"] > a - 1) & (lineas["nu0"] < b + 1)
            if not m.any():
                continue
            sub = {kk: v[m] for kk, v in lineas.items()}
            S, c, gL, aD = G.parametros(sub, T, p, 0.0, 1.0, 1.0, inf["masa_u"])
            u_col = float((col["vmr"][gas] * col["n_aire"]).sum())
            fuerte = sub["S"] * u_col >= UMBRAL_LBL
            cs.append(c[fuerte]); gls.append(gL[fuerte]); ads.append(aD[fuerte])
    if cs:
        return G.rejilla_adaptada(a, b, np.concatenate(cs), np.concatenate(gls), np.concatenate(ads))
    return G.rejilla_adaptada(a, b, np.zeros(0), np.zeros(0), np.zeros(0))


def espesores(datos, nu, col, mtckd):
    """{gas: (n, n_nu)} espesor optico de cada capa de la columna, linea a linea."""
    n = len(col["Tm"])
    taus = {}
    for gas, (lineas, isotopos) in datos.items():
        if lineas is None:
            continue
        u = col["vmr"][gas] * col["n_aire"]
        k = np.zeros((n, len(nu)))
        for c in range(n):
            k[c] = k_linea_a_linea(lineas, isotopos, nu, col["pm"][c], col["Tm"][c], col["vmr"][gas][c],
                                   float(u.sum()), gas == "H2O", mtckd)
        taus[gas] = k * u[:, None]
    return taus

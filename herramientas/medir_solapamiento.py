# medir_solapamiento.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7, decision T9): MIDE, banda a banda del infrarrojo,
# el error de cada forma de tratar el SOLAPAMIENTO de los gases frente al calculo linea a linea, en una columna de
# referencia de P3N. Diagnostico: no cambia nada del modelo.
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python herramientas/medir_solapamiento.py --hilos 4             # las 16 bandas del infrarrojo (horas)
#     python herramientas/medir_solapamiento.py --bandas 3 4 5        # solo algunas
# Columnas: 'ref' (suelo a 288 K) y, con el MISMO agrupamiento del mapa hecho con 'ref', 'tropical' (300 K, mas vapor) y
# 'polar' (255 K, casi seco): miden si un agrupamiento fijo sirve para atmosferas muy distintas (--columnas ...).
# Opciones: --hitran CARPETA (por defecto ~/Documentos/B3N/datos/hitran; la de descargar_hitran.py), --gases ...,
# --capas N (30), --hilos N (bandas a la vez; 1), --interpolar, --salida CARPETA (outputs/solapamiento).
#
# Columnas de prueba (solo para medir; NO son el clima de P3N): las de fase52_referencia.py (niveles de M3N, composicion
# sorteada de P3N, vapor de Manabe y Wetherald y una campana de ozono).
# Metodos que compara (los flujos con el mismo dos flujos de fase52_infrarrojo):
#   LBL      linea a linea (cada numero de onda con su Planck): la referencia;
#   PM       k de la MEZCLA ordenado capa a capa (k correlacionado "premezclado"): el error propio de las tablas k;
#   PM+fr    igual, con la fraccion de Planck de cada punto g (en vez de repartir la emision de la banda por igual);
#   RORR     cada gas con su tabla y solapamiento aleatorio con reordenacion y reagrupamiento;
#   SUMA     cada gas con su tabla, sumando como si estuvieran correlacionados (para ver por que hace falta tratarlo);
#   MAPA     agrupamiento fijo de numeros de onda (como ecCKD) con media aritmetica de cada gas en cada grupo, y su
#            fraccion de Planck exacta; MAPA-log, con media logaritmica; MAPA-tr, con la media de la transmision de cada
#            gas pesada con Planck (la de ecCKD por defecto); MAPA-trR, igual pero con la cantidad de gas de la columna
#            de referencia en cada capa, como lo haria una tabla (en la columna 'ref' coincide con MAPA-tr);
#            MAPA-trV, como trR pero el vapor con su propia cantidad (una tabla con la dimension del vapor); MAPA32-trV,
#            lo mismo con 32 grupos;
#            con --interpolar, TABLA*: como MAPA-trR pero interpolando en la rejilla (p, T) de las tablas (T15), con
#            varias formas de interpolar (ver tabla_interpolada).
# Salida: informe.txt (y resultados.json) en la carpeta de salida: por banda, errores de OLR, del infrarrojo que llega
# al suelo y del calentamiento (K/dia de 86 400 s) en la troposfera (p > 100 hPa) y por encima.

import argparse
import json
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, AQUI)

import fase52_continuo as C          # noqa: E402,F401
import fase52_generador as G         # noqa: E402,F401
import fase52_hitran as H            # noqa: E402
import fase52_infrarrojo as F        # noqa: E402
import fase52_k as K                 # noqa: E402
from fase52_referencia import (GASES, VMR, M_AIRE, CP, G_P3N, UMBRAL_LBL, COLUMNAS, columna, cargar,  # noqa: E402
                               k_linea_a_linea, rejilla_comun, espesores)

def flujos_lbl(tau, nu, dnu, col, trozo=100000):
    """Linea a linea: cada numero de onda con su Planck (W/m^2/cm^-1 * dnu)."""
    n = tau.shape[0]
    U = np.zeros(n + 1)
    D = np.zeros(n + 1)
    for i0 in range(0, len(nu), trozo):
        s = slice(i0, i0 + trozo)
        x = nu[s][:, None]
        Bh = np.pi * F.C1 * x ** 3 / np.expm1(F.C2 * x / col["Th"][None, :])
        u, d = F.dos_flujos(tau[:, s].T, dnu[s], Bh, Bh[:, -1])
        U += u
        D += d
    return U, D


def flujos_terminos(tau_g, pesos, Bh_banda, fr=None):
    """Con tablas k: tau_g (n, puntos g); emision de la banda repartida por los pesos (fr = None) o por las fracciones
    de Planck fr (puntos g, n+1)."""
    if fr is None:
        Bh = Bh_banda[None, :]
    else:
        Bh = fr * Bh_banda[None, :] / np.asarray(pesos)[:, None]
    return F.dos_flujos(tau_g.T, pesos, Bh, Bh[:, -1] if Bh.shape[0] > 1 else Bh[0, -1])


def fracciones_planck_mezcla(tau, nu, dnu, col, pesos):
    """Fraccion de la emision de la banda en cada punto g, con el orden de la mezcla de la capa de cada seminivel."""
    bordes = np.concatenate([[0.0], np.cumsum(pesos)])
    n = tau.shape[0]
    fr = np.zeros((len(pesos), n + 1))
    for h in range(n + 1):
        capa = min(max(h - 1, 0), n - 1) if h > 0 else 0
        orden = np.argsort(tau[capa], kind="stable")
        Gc = np.cumsum(dnu[orden]); Gc = (Gc - 0.5 * dnu[orden]) / Gc[-1]
        grupo = np.clip(np.searchsorted(bordes, Gc, side="right") - 1, 0, len(pesos) - 1)
        B = nu[orden] ** 3 / np.expm1(F.C2 * nu[orden] / col["Th"][h]) * dnu[orden]
        fr[:, h] = np.bincount(grupo, weights=B, minlength=len(pesos)) / B.sum()
    return fr


def medir_banda(datos, banda, col, mtckd, col_ref=None, rejilla=None, taus_ref=None, taus=None, interpolar=False):
    """Errores de cada metodo en la columna col. El mapa (MAPA*) se hace con la columna de referencia col_ref (por
    defecto, la misma col): con otra columna se mide si un agrupamiento fijo sirve para atmosferas distintas."""
    nu, dnu = rejilla if rejilla is not None else rejilla_comun(datos, banda, col_ref or col)
    n = len(col["Tm"])
    if taus is None:
        taus = espesores(datos, nu, col, mtckd)
    tau = sum(taus.values())
    g, pesos = K.nodos_g()
    Bh_banda = F.planck_integral(banda[0], banda[1], col["Th"])
    res = {"puntos": int(len(nu))}
    res["LBL"] = flujos_lbl(tau, nu, dnu, col)
    # PM: la mezcla ordenada capa a capa
    tau_pm = np.array([K.k_en_g(*K.distribucion_k(tau[c], dnu), g) for c in range(n)])
    res["PM"] = flujos_terminos(tau_pm, pesos, Bh_banda)
    res["PM+fr"] = flujos_terminos(tau_pm, pesos, Bh_banda, fracciones_planck_mezcla(tau, nu, dnu, col, pesos))
    # cada gas con su tabla
    por_gas = [np.array([K.k_en_g(*K.distribucion_k(t[c], dnu), g) for c in range(n)]) for t in taus.values()]
    res["RORR"] = flujos_terminos(F.combinar_rorr(por_gas, pesos), pesos, Bh_banda)
    res["SUMA"] = flujos_terminos(sum(por_gas), pesos, Bh_banda)
    # mapa espectral fijo
    if col_ref is None:
        tau_ref = tau
    else:
        tau_ref = sum((taus_ref if taus_ref is not None else espesores(datos, nu, col_ref, mtckd)).values())
    acum = np.vstack([np.zeros(len(nu)), np.cumsum(tau_ref, axis=0)])
    grupo = F.tabla_mapa(acum, dnu, pesos)
    B = nu[None, :] ** 3 / np.expm1(F.C2 * nu[None, :] / col["Th"][:, None]) * dnu[None, :]
    fr = np.array([B[:, grupo == j].sum(axis=1) for j in range(len(pesos))]) / B.sum(axis=1)[None, :]
    pesos_mapa = np.array([dnu[grupo == j].sum() for j in range(len(pesos))]) / dnu.sum()
    for tipo, nombre in (("aritmetica", "MAPA"), ("logaritmica", "MAPA-log")):
        tg = sum(F.media_por_grupo(t, dnu, grupo, len(pesos), tipo) for t in taus.values())
        ok = pesos_mapa > 0
        res[nombre] = flujos_terminos(tg[:, ok], pesos_mapa[ok], Bh_banda, fr[ok])
    # MAPA-tr: cada gas con su transmision media en cada grupo (pesada con Planck de la capa), y se suman
    Bm = nu[None, :] ** 3 / np.expm1(F.C2 * nu[None, :] / col["Tm"][:, None]) * dnu[None, :]
    tg = sum(F.media_transmision_por_grupo(t, Bm, grupo, len(pesos)) for t in taus.values())
    ok = pesos_mapa > 0
    res["MAPA-tr"] = flujos_terminos(tg[:, ok], pesos_mapa[ok], Bh_banda, fr[ok])
    # MAPA-trR: como lo haria una tabla: la media de la transmision se hace con la cantidad de gas de la columna de
    # REFERENCIA en esa capa (u_ref), sobre el k(nu) de esta columna, y despues se multiplica por la cantidad de esta
    if col_ref is not None:
        tg = 0.0
        for gas, t in taus.items():
            u = col["vmr"][gas] * col["n_aire"]
            u_ref = col_ref["vmr"][gas] * col_ref["n_aire"]
            f = (u_ref / u)[:, None]
            tg = tg + F.media_transmision_por_grupo(t * f, Bm, grupo, len(pesos)) / f
        res["MAPA-trR"] = flujos_terminos(tg[:, ok], pesos_mapa[ok], Bh_banda, fr[ok])
        # MAPA-trV: como trR, pero el vapor promediado con SU cantidad (como haria una tabla con la dimension de la
        # fraccion de vapor en la que la cantidad del promedio es la de ese vapor: exacto en sus nodos)
        tg = 0.0
        for gas, t in taus.items():
            if gas == "H2O":
                tg = tg + F.media_transmision_por_grupo(t, Bm, grupo, len(pesos))
            else:
                f = ((col_ref["vmr"][gas] * col_ref["n_aire"]) / (col["vmr"][gas] * col["n_aire"]))[:, None]
                tg = tg + F.media_transmision_por_grupo(t * f, Bm, grupo, len(pesos)) / f
        res["MAPA-trV"] = flujos_terminos(tg[:, ok], pesos_mapa[ok], Bh_banda, fr[ok])
        # MAPA32-trV: igual con 32 grupos (el doble de puntos de Gauss en cada tramo de T6)
        g32, p32 = K.nodos_g(tuple((a_, b_, 2 * n_) for a_, b_, n_ in K.TRAMOS_G))
        grupo32 = F.tabla_mapa(acum, dnu, p32)
        pm32 = np.bincount(grupo32, weights=dnu, minlength=len(p32)) / dnu.sum()
        fr32 = np.array([B[:, grupo32 == j].sum(axis=1) for j in range(len(p32))]) / B.sum(axis=1)[None, :]
        tg = 0.0
        for gas, t in taus.items():
            if gas == "H2O":
                tg = tg + F.media_transmision_por_grupo(t, Bm, grupo32, len(p32))
            else:
                f = ((col_ref["vmr"][gas] * col_ref["n_aire"]) / (col["vmr"][gas] * col["n_aire"]))[:, None]
                tg = tg + F.media_transmision_por_grupo(t * f, Bm, grupo32, len(p32)) / f
        ok32 = pm32 > 0
        res["MAPA32-trV"] = flujos_terminos(tg[:, ok32], pm32[ok32], Bh_banda, fr32[ok32])
        if interpolar:
            for nombre, tg in tabla_interpolada(datos, nu, dnu, grupo, len(pesos), col, col_ref, mtckd).items():
                res[nombre] = flujos_terminos(tg[:, ok], pesos_mapa[ok], Bh_banda, fr[ok])
    return res


def tabla_interpolada(datos, nu, dnu, grupo, n_grupos, col, col_ref, mtckd, variantes=None):
    """Como MAPA-trR, pero cada k sale de INTERPOLAR en una rejilla de tabla: en cada capa se calcula el k de cada grupo
    en las esquinas (p, T) que la rodean (con la cantidad de referencia de esa capa y Planck a la T de la esquina) y se
    interpola en ln p y en la temperatura. Variantes (nombre: (temperaturas de la rejilla, k o ln k, eje de T)):
      TABLA       presiones de la tabla (generar_tablas_k.P_HPA) y temperaturas cada 30 K, lineal en k y en T;
      TABLA-log   la misma, lineal en ln k;
      TABLA-1/T   la misma, ln k lineal en 1/T;
      TABLA-15K   temperaturas cada 15 K (150-330, la rejilla de las tablas desde la 3.23.0), ln k lineal en T.
    Sin la dimension del vapor (se usa la fraccion molar de la capa)."""
    import generar_tablas_k as GT
    P = GT.P_HPA[::-1]                                  # presion creciente
    T30 = np.arange(150.0, 331.0, 30.0)                 # la rejilla T15 de la 3.22.0
    T15 = np.arange(150.0, 331.0, 15.0)                 # la de la 3.23.0 (generar_tablas_k.T_K)
    if variantes is None:
        variantes = {"TABLA": (T30, "lin", "T"), "TABLA-log": (T30, "log", "T"),
                     "TABLA-1/T": (T30, "log", "1/T"), "TABLA-15K": (T15, "log", "T")}
    n = len(col["Tm"])
    out = {v: np.zeros((n, n_grupos)) for v in variantes}
    for gas, (lineas, isotopos) in datos.items():
        if lineas is None:
            continue
        u = col["vmr"][gas] * col["n_aire"]
        u_ref = col_ref["vmr"][gas] * col_ref["n_aire"]
        for c in range(n):
            ip = int(np.clip(np.searchsorted(P, col["pm"][c]) - 1, 0, len(P) - 2))
            fp = (np.log(col["pm"][c]) - np.log(P[ip])) / (np.log(P[ip + 1]) - np.log(P[ip]))
            cache = {}

            def esquina(pp, Tc):
                if (pp, Tc) not in cache:
                    k = k_linea_a_linea(lineas, isotopos, nu, pp, Tc, col["vmr"][gas][c], float(u.sum()),
                                        gas == "H2O", mtckd)
                    peso = nu ** 3 / np.expm1(F.C2 * nu / Tc) * dnu
                    cache[pp, Tc] = F.media_transmision_por_grupo((k * u_ref[c])[None, :], peso, grupo,
                                                                  n_grupos)[0] / u_ref[c]
                return cache[pp, Tc]
            for nombre, (TT, modo, eje) in variantes.items():
                it = int(np.clip(np.searchsorted(TT, col["Tm"][c]) - 1, 0, len(TT) - 2))
                T0, T1, Tx = TT[it], TT[it + 1], col["Tm"][c]
                ft = (1 / Tx - 1 / T0) / (1 / T1 - 1 / T0) if eje == "1/T" else (Tx - T0) / (T1 - T0)
                w = {(0, 0): (1 - fp) * (1 - ft), (0, 1): (1 - fp) * ft, (1, 0): fp * (1 - ft), (1, 1): fp * ft}
                e = {(a, b): esquina(P[ip + a], TT[it + b]) for a in (0, 1) for b in (0, 1)}
                if modo == "lin":
                    val = sum(w[q] * e[q] for q in w)
                else:
                    val = np.exp(sum(w[q] * np.log(np.maximum(e[q], 1e-300)) for q in w))
                out[nombre][c] += val * u[c]
    return out


def resumen(res, col):
    U0, D0 = res["LBL"]
    h0 = F.calentamiento(U0, D0) / (col["dp"] / G_P3N * CP) * 86400.0
    trop = col["pm"] > 100.0
    out = {}
    for m, v in res.items():
        if m == "puntos":
            continue
        U, D = v
        h = F.calentamiento(U, D) / (col["dp"] / G_P3N * CP) * 86400.0
        out[m] = {"OLR": float(U[0]), "dOLR": float(U[0] - U0[0]), "DLR": float(D[-1]), "dDLR": float(D[-1] - D0[-1]),
                  "dQ_trop_max": float(np.max(np.abs(h - h0)[trop])), "dQ_estr_max": float(np.max(np.abs(h - h0)[~trop])),
                  "dQ_rms": float(np.sqrt(np.mean((h - h0) ** 2)))}
    return out


def procesar_banda(args):
    """Una banda con todas las columnas (para repartir las bandas entre procesos)."""
    ib, hitran, gases, capas, columnas, interpolar = args
    t0 = time.time()
    banda = H.BANDAS_IR[ib]
    ref = columna(capas, Ts=COLUMNAS["ref"])
    ruta = os.path.join(hitran, "absco-ref_wv-mt-ckd.nc")
    mtckd = C.leer_mtckd(ruta) if os.path.exists(ruta) else None
    datos = {gas: cargar(hitran, gas, *banda) for gas in gases}
    rejilla = rejilla_comun(datos, banda, ref)
    taus_ref = espesores(datos, rejilla[0], ref, mtckd)
    resultados, texto = {}, []
    for nc in columnas:
        col = ref if nc == "ref" else columna(capas, Ts=COLUMNAS[nc])
        res = medir_banda(datos, banda, col, mtckd, col_ref=ref, rejilla=rejilla, taus_ref=taus_ref,
                          taus=taus_ref if nc == "ref" else None, interpolar=interpolar)
        r = resumen(res, col)
        resultados[nc] = {"puntos": res["puntos"], "metodos": r}
        texto += [f"\nBanda {ib} ({banda[0]}-{banda[1]} cm-1), columna {nc} (suelo {COLUMNAS[nc]:.0f} K): "
                  f"{res['puntos']} numeros de onda. LBL: OLR {r['LBL']['OLR']:.3f}, al suelo {r['LBL']['DLR']:.3f} W/m2",
                  f"  {'metodo':9s} {'dOLR':>8s} {'dSuelo':>8s} {'dQ trop':>8s} {'dQ arr.':>8s} {'dQ rms':>8s}"
                  "   (W/m2; K/dia)"]
        for m, v in r.items():
            if m != "LBL":
                texto.append(f"  {m:9s} {v['dOLR']:8.3f} {v['dDLR']:8.3f} {v['dQ_trop_max']:8.3f} "
                             f"{v['dQ_estr_max']:8.3f} {v['dQ_rms']:8.3f}")
    texto.append(f"  ({time.time() - t0:.0f} s)")
    return ib, resultados, texto, mtckd is not None


def main():
    ap = argparse.ArgumentParser(description="Error del solapamiento de gases frente al linea a linea (Fase 5.2, T9)")
    ap.add_argument("--hitran", default=os.path.expanduser(os.path.join("~", "Documentos", "B3N", "datos", "hitran")))
    ap.add_argument("--bandas", type=int, nargs="+", default=list(range(len(H.BANDAS_IR))))
    ap.add_argument("--gases", nargs="+", default=list(GASES))
    ap.add_argument("--columnas", nargs="+", default=list(COLUMNAS), choices=list(COLUMNAS))
    ap.add_argument("--capas", type=int, default=30)
    ap.add_argument("--interpolar", action="store_true",
                    help="mide tambien el error de interpolar en la rejilla de las tablas (T15); 4 veces mas lento")
    ap.add_argument("--hilos", type=int, default=1,
                    help="bandas a la vez (cada una puede necesitar 2-4 GB de memoria)")
    ap.add_argument("--salida", default=os.path.join("outputs", "solapamiento"))
    a = ap.parse_args()
    os.makedirs(a.salida, exist_ok=True)
    tareas = [(ib, a.hitran, a.gases, a.capas, a.columnas, a.interpolar) for ib in a.bandas]
    todo, por_banda, hubo_mtckd = {}, {}, False
    print(f"Midiendo {len(tareas)} bandas con {a.hilos} proceso(s); cada banda tarda de minutos a mas de una hora.",
          flush=True)

    def guardar(ib, resultados, texto, con_mtckd):
        nonlocal hubo_mtckd
        hubo_mtckd = hubo_mtckd or con_mtckd
        banda = H.BANDAS_IR[ib]
        for nc, v in resultados.items():
            todo.setdefault(nc, {})[f"{banda[0]}-{banda[1]}"] = v
        por_banda[ib] = texto
        print("\n".join(texto), flush=True)
        with open(os.path.join(a.salida, "resultados.json"), "w", encoding="utf-8") as f:
            json.dump(todo, f, indent=1)

    if a.hilos > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.hilos) as pool:
            for r in pool.imap_unordered(procesar_banda, tareas):
                guardar(*r)
    else:
        for t in tareas:
            guardar(*procesar_banda(t))
    inf = [l for ib in a.bandas for l in por_banda[ib]]
    inf.append("\nSuma de las bandas medidas (errores de flujo, W/m2) y peor error de calentamiento de una banda (K/dia):")
    for nc in a.columnas:
        inf.append(f"  columna {nc}:")
        tot = {}
        for v in todo.get(nc, {}).values():
            for m, x in v["metodos"].items():
                t = tot.setdefault(m, {"dOLR": 0.0, "dDLR": 0.0, "dQ": 0.0})
                t["dOLR"] += x["dOLR"]; t["dDLR"] += x["dDLR"]
                t["dQ"] = max(t["dQ"], x["dQ_trop_max"], x["dQ_estr_max"])
        for m, t in tot.items():
            if m != "LBL":
                inf.append(f"    {m:9s} dOLR {t['dOLR']:8.3f}   dSuelo {t['dDLR']:8.3f}   peor dQ {t['dQ']:6.3f}")
    with open(os.path.join(a.salida, "informe.txt"), "w", encoding="utf-8") as f:
        f.write("SOLAPAMIENTO DE GASES (Fase 5.2, T9) -- errores frente al linea a linea\n")
        f.write(f"HITRAN: {a.hitran}; gases: {a.gases}; capas: {a.capas}; MT_CKD: {'si' if hubo_mtckd else 'no'}\n")
        f.write("El mapa (MAPA*) se hace siempre con la columna 'ref'; las demas columnas miden si sirve fuera de ella.\n")
        f.write("\n".join(inf) + "\n")
    print(f"\nGuardado en {os.path.join(a.salida, 'informe.txt')}: pasale ese archivo a Claude.")


if __name__ == "__main__":
    main()

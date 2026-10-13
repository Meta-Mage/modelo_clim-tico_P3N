# generar_tablas_k.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7): calcula en el PC las TABLAS k del INFRARROJO de cada
# gas a partir de las lineas de HITRAN descargadas (herramientas/descargar_hitran.py). Calculo largo (horas): con barra
# de avance, usa los hilos del PC y se puede cortar y continuar (guarda cada casilla terminada).
#
# 3.23.0 (decision tecnica T9): tablas con MAPA ESPECTRAL FIJO por banda (fase52_mapa, como ecCKD), en vez de una
# distribucion k ordenada por separado en cada casilla. Primero se hace el mapa de cada banda con una columna de
# referencia de P3N (mapas/mapa_bNN.npz); despues, cada casilla promedia cada gas en los grupos del mapa. Las tablas de
# la luz se haran con el espectro de S3N (pendiente): esta version solo calcula las 16 bandas del infrarrojo.
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python herramientas/generar_tablas_k.py --prueba        # minutos: 4 bandas y 3 alturas; estima el total
#     python herramientas/generar_tablas_k.py                 # todo (se puede cortar con Ctrl+C y relanzar)
# Opciones: --hitran CARPETA (por defecto ~/Documentos/B3N/datos/hitran), --salida CARPETA (por defecto
# ~/Documentos/B3N/datos/tablas_k), --gases H2O CO2 ..., --bandas 0 1 ... (solo esas), --hilos N (por defecto, todos
# menos uno).
#
# Rejilla de la tabla (decision tecnica T15, PROVISIONAL; el error de interpolar en ella se mide con
# herramientas/medir_solapamiento.py --interpolar):
#   presion: 36 niveles de 1050 a 0,01 hPa, equiespaciados en log p;
#   temperatura: 150 a 330 K cada 15 K (13). Medido con medir_solapamiento.py --interpolar (banda de 820-980 cm^-1,
#   la del continuo del vapor): cada 30 K, 0,28-0,29 W/m2 de error al suelo interpolando ln k en T; cada 15 K, 0,05-0,13;
#   fraccion molar del vapor (su ensanchamiento consigo mismo, su continuo propio y la cantidad con la que se promedia
#   en cada grupo, fase52_mapa): 1e-6, 1e-5, 1e-4, 1e-3, 4e-3, 1,5e-2 y 6e-2 (se interpola en ln); por encima de
#   100 hPa solo las tres primeras (alli el vapor es < 1e-4) y el resto se copia;
#   bandas: las 16 del infrarrojo de fase52_hitran (T8).
# Vapor: lineas sin el pedestal (T4) MAS el continuo MT_CKD 4.3 (fase52_continuo), en la misma tabla.
# Lineas descartadas (T14): umbral 1e-7 cm^-1 por linea con la columna maxima del gas (U_MAX); la suma de lo
# descartado en cada casilla queda en la tabla.
# Salida: tabla_k_<GAS>.npz con k (bandas, presiones, temperaturas, vapor, grupos) en cm^2/molecula, y
# mapas_infrarrojo.npz con el peso de cada grupo y su fraccion de Planck de 100 a 350 K (cada 1 K).

import argparse
import json
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, AQUI)

GASES = {"H2O": 1, "CO2": 2, "O3": 3, "N2O": 4, "CH4": 6, "O2": 7}
P_HPA = np.logspace(np.log10(1050.0), np.log10(0.01), 36)
T_K = np.arange(150.0, 331.0, 15.0)        # 13 (T15 revisada el 13/10/2026: antes cada 30 K)
VMR_H2O = (1e-6, 1e-5, 1e-4, 1e-3, 4e-3, 1.5e-2, 6e-2)
P_VMR_MIN = 100.0   # hPa
VMR_ARRIBA = 3       # por encima de P_VMR_MIN solo las 3 primeras (hasta 1e-4); el resto se copia
T_FR = np.arange(100.0, 350.5, 1.0)       # temperaturas de la tabla de fracciones de Planck
# columna maxima de cada gas (moleculas/cm^2) para el descarte: aire de 1050 hPa con g = 8 m/s^2 (cota por arriba de
# la columna: P3N tiene g = 9,05, parametros.P3N_GRAVEDAD) por una fraccion molar maxima generosa (mucho mayor que la
# sorteada, para que la tabla sirva tambien si se cambia la composicion por worldbuilding)
N_AIRE = 1050e2 / (28.96e-3 / 6.02214076e23 * 8.0) / 1e4
U_MAX = {"H2O": 0.05 * N_AIRE, "CO2": 2e-3 * N_AIRE, "O3": 2e-5 * N_AIRE, "N2O": 2e-6 * N_AIRE,
         "CH4": 5e-6 * N_AIRE, "O2": 0.3 * N_AIRE}
UMBRAL = 1e-7
BANDAS_PRUEBA = (0, 3, 6, 13)
P_PRUEBA = (0, 18, 35)


def tareas(gases, bandas, prueba):
    out = []
    for g in gases:
        vmrs = VMR_H2O if g == "H2O" else (0.0,)
        for ib in bandas:
            for ip, p in enumerate(P_HPA):
                if prueba and (ib not in BANDAS_PRUEBA or ip not in P_PRUEBA):   # una muestra de bandas y alturas
                    continue
                out.append((g, ib, ip, float(p), vmrs))
    return out


def _isotopos(hitran, gas):
    with open(os.path.join(hitran, "tips.json"), encoding="utf-8") as f:
        tips = json.load(f)
    with open(os.path.join(hitran, "isotopologos.json"), encoding="utf-8") as f:
        iso = json.load(f)
    m = GASES[gas]
    out = {}
    for clave, inf in iso.items():
        mm, ii = (int(x) for x in clave.split("_"))
        if mm == m and tips["Q"].get(clave):
            out[ii] = {"masa_u": inf["masa_u"], "T": np.array(tips["T_K"], float), "Q": np.array(tips["Q"][clave])}
    return out


_LINEAS = {}


def _lineas_banda(hitran, gas, banda):
    """Lineas del gas que afectan a la banda (su centro a menos de 25 cm^-1), leidas una vez por proceso."""
    import fase52_hitran as H
    clave = (gas, banda)
    if clave not in _LINEAS:
        ruta = os.path.join(hitran, gas + ".data")
        lineas = H.leer_par(ruta, banda[0] - 25.0, banda[1] + 25.0) if os.path.exists(ruta) else None
        if lineas is not None and len(lineas["nu0"]) == 0:
            lineas = None
        _LINEAS.clear()                      # una banda cada vez en memoria
        _LINEAS[clave] = (lineas, _isotopos(hitran, gas))
    return _LINEAS[clave]


def _mtckd(hitran):
    import fase52_continuo as C
    ruta = os.path.join(hitran, "absco-ref_wv-mt-ckd.nc")
    return C.leer_mtckd(ruta) if os.path.exists(ruta) else None


def hacer_mapa(args_tarea):
    """Fase 1: el mapa de una banda (todos los gases, columna de referencia). Se guarda en mapas/mapa_bNN.npz."""
    hitran, salida, gases, ib = args_tarea
    import fase52_hitran as H
    import fase52_mapa as MP
    fichero = os.path.join(salida, "mapas", f"mapa_b{ib:02d}.npz")
    if os.path.exists(fichero):
        return fichero, 0.0
    t0 = time.time()
    banda = H.BANDAS_IR[ib]
    datos = {}
    for g in gases:
        if g in ("O2",):                       # sin lineas que cuenten en el infrarrojo termico de P3N
            continue
        lineas, iso = _lineas_banda(hitran, g, banda)
        datos[g] = (lineas, iso)
        _LINEAS.clear()
    mp_ = MP.mapa_banda(datos, banda, _mtckd(hitran))
    np.savez_compressed(fichero + ".tmp.npz", nu=mp_["nu"], dnu=mp_["dnu"], grupo=mp_["grupo"], pesos=mp_["pesos"],
                        pm_ref=mp_["pm_ref"], banda=np.array(banda, float),
                        **{f"u_ref_{g}": v for g, v in mp_["u_ref"].items()})
    os.replace(fichero + ".tmp.npz", fichero)
    return fichero, time.time() - t0


def cargar_mapa(salida, ib):
    d = np.load(os.path.join(salida, "mapas", f"mapa_b{ib:02d}.npz"))
    return {"banda": tuple(d["banda"]), "nu": d["nu"], "dnu": d["dnu"], "grupo": d["grupo"], "pesos": d["pesos"],
            "pm_ref": d["pm_ref"], "u_ref": {k[6:]: d[k] for k in d.files if k.startswith("u_ref_")}}


def calcular(args_tarea):
    """Fase 2: una casilla (gas, banda, presion) con todas las temperaturas y fracciones de vapor."""
    hitran, salida, (gas, ib, ip, p, vmrs) = args_tarea
    import fase52_hitran as H
    import fase52_mapa as MP
    fichero = os.path.join(salida, "casillas_mapa", f"{gas}_b{ib:02d}_p{ip:02d}.npz")
    if os.path.exists(fichero):
        return fichero, 0.0
    t0 = time.time()
    mapa = cargar_mapa(salida, ib)
    lineas, isotopos = _lineas_banda(hitran, gas, H.BANDAS_IR[ib])
    mtckd = _mtckd(hitran) if gas == "H2O" else None
    n_g = len(mapa["pesos"])
    k = np.zeros((len(T_K), len(vmrs), n_g))
    perdida = np.zeros((len(T_K), len(vmrs)))
    for it, T in enumerate(T_K):
        for iv, v in enumerate(vmrs):
            if iv >= VMR_ARRIBA and p < P_VMR_MIN:   # copia de la ultima calculada (ver la rejilla, arriba)
                k[it, iv], perdida[it, iv] = k[it, VMR_ARRIBA - 1], perdida[it, VMR_ARRIBA - 1]
                continue
            k[it, iv], perdida[it, iv] = MP.k_grupos(lineas, isotopos, mapa, gas, p, float(T), v, U_MAX[gas], UMBRAL,
                                                     mtckd)
    np.savez(fichero + ".tmp.npz", k=k, perdida=perdida)
    os.replace(fichero + ".tmp.npz", fichero)
    return fichero, time.time() - t0


def _barra(hechas, total, t0, hechas_ahora, que):
    transc = time.time() - t0
    quedan = (total - hechas) * transc / max(hechas_ahora, 1)
    b = "#" * int(25 * hechas / max(total, 1))
    sys.stdout.write(f"\r{que} [{b:25s}] {hechas}/{total} | {transc / 60:6.1f} min | quedan ~{quedan / 3600:5.1f} h ")
    sys.stdout.flush()


def main():
    import fase52_hitran as H
    ap = argparse.ArgumentParser(description="Tablas k del infrarrojo de M3N (Fase 5.2) a partir de HITRAN")
    ap.add_argument("--hitran", default=os.path.expanduser("~/Documentos/B3N/datos/hitran"))
    ap.add_argument("--salida", default=os.path.expanduser("~/Documentos/B3N/datos/tablas_k"))
    ap.add_argument("--gases", nargs="+", default=[g for g in GASES if g != "O2"])
    ap.add_argument("--bandas", type=int, nargs="+", default=list(range(len(H.BANDAS_IR))))
    ap.add_argument("--hilos", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--prueba", action="store_true", help="solo 4 bandas y 3 presiones, para medir el tiempo")
    a = ap.parse_args()
    for falta in ("isotopologos.json", "tips.json", "absco-ref_wv-mt-ckd.nc"):
        if not os.path.exists(os.path.join(a.hitran, falta)):
            sys.exit(f"Falta {falta}: vuelve a lanzar herramientas/descargar_hitran.py (no repite lo ya descargado).")
    bandas = [b for b in a.bandas if not a.prueba or b in BANDAS_PRUEBA]
    os.makedirs(os.path.join(a.salida, "mapas"), exist_ok=True)
    os.makedirs(os.path.join(a.salida, "casillas_mapa"), exist_ok=True)
    import multiprocessing as mp
    # fase 1: mapas
    t0, hechas, ahora = time.time(), 0, 0
    # los mapas necesitan mas memoria (hasta ~2-3 GB por banda): como mucho 3 a la vez
    with mp.get_context("spawn").Pool(max(1, min(a.hilos, len(bandas), 3))) as pool:
        for f, dt in pool.imap_unordered(hacer_mapa, [(a.hitran, a.salida, a.gases, ib) for ib in bandas]):
            hechas += 1
            ahora += dt > 0
            _barra(hechas, len(bandas), t0, ahora, "mapas   ")
    print(flush=True)
    # fase 2: casillas
    lista = tareas(a.gases, bandas, a.prueba)
    hechas = sum(os.path.exists(os.path.join(a.salida, "casillas_mapa", f"{g}_b{ib:02d}_p{ip:02d}.npz"))
                 for g, ib, ip, _, _ in lista)
    print(f"{len(lista)} casillas ({hechas} ya hechas); {a.hilos} hilos", flush=True)
    t0, ahora, tiempos = time.time(), 0, []
    with mp.get_context("spawn").Pool(a.hilos) as pool:
        for f, dt in pool.imap_unordered(calcular, [(a.hitran, a.salida, t) for t in lista]):
            hechas += 1
            if dt > 0:
                ahora += 1
                tiempos.append(dt)
            _barra(hechas, len(lista), t0, ahora, "casillas")
    print()
    if a.prueba:
        n_total = len(tareas(a.gases, list(range(len(H.BANDAS_IR))), False))
        media = float(np.mean(tiempos)) if tiempos else 0.0
        print(f"PRUEBA: {media:.1f} s por casilla de media con {a.hilos} hilos en paralelo; el total ({n_total} casillas"
              f" y 16 mapas) llevaria del orden de {n_total * media / a.hilos / 3600:.1f} h (orientativo: las bandas "
              "difieren).")
        return
    # juntar
    import fase52_mapa as MP
    mapas = [cargar_mapa(a.salida, ib) for ib in bandas]
    n_g = len(mapas[0]["pesos"])
    pesos = np.array([m["pesos"] for m in mapas])
    fr = np.array([MP.fracciones_planck(m, T_FR) for m in mapas])
    np.savez_compressed(os.path.join(a.salida, "mapas_infrarrojo.npz"), bandas=np.array([H.BANDAS_IR[i] for i in bandas]),
                        indices_banda=np.array(bandas), pesos=pesos, T_fr=T_FR, fraccion_planck=fr, metodo="mapa (T9)")
    for g in a.gases:
        vmrs = VMR_H2O if g == "H2O" else (0.0,)
        k = np.zeros((len(bandas), len(P_HPA), len(T_K), len(vmrs), n_g), dtype=np.float64)
        perd = np.zeros((len(bandas), len(P_HPA), len(T_K), len(vmrs)))
        for jb, ib in enumerate(bandas):
            for ip in range(len(P_HPA)):
                d = np.load(os.path.join(a.salida, "casillas_mapa", f"{g}_b{ib:02d}_p{ip:02d}.npz"))
                k[jb, ip], perd[jb, ip] = d["k"], d["perdida"]
        np.savez_compressed(os.path.join(a.salida, f"tabla_k_{g}.npz"), k=k, perdida_cm1=perd,
                            bandas=np.array([H.BANDAS_IR[i] for i in bandas]), indices_banda=np.array(bandas),
                            p_hPa=P_HPA, T_K=T_K, vmr=np.array(vmrs), umbral=UMBRAL, u_max=U_MAX[g],
                            metodo="mapa (T9)")
        print(f"{g}: tabla_k_{g}.npz | absorcion descartada maxima en una casilla {perd.max():.2e} cm^-1", flush=True)
    print(f"Listo en {a.salida}. Pasale a Claude lo que ha salido en la terminal.")


if __name__ == "__main__":
    main()

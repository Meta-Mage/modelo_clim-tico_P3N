# generar_tablas_k.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7): calcula en el PC las TABLAS k de cada gas a partir
# de las lineas de HITRAN descargadas (herramientas/descargar_hitran.py). Calculo largo (horas): con barra de avance,
# usa los hilos del PC y se puede cortar y continuar (guarda cada casilla terminada).
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python herramientas/generar_tablas_k.py --prueba        # minutos: 4 bandas y 3 alturas; estima el total
#     python herramientas/generar_tablas_k.py                 # todo (se puede cortar con Ctrl+C y relanzar)
# Opciones: --hitran CARPETA (por defecto ~/Documentos/B3N/datos/hitran), --salida CARPETA (por defecto
# ~/Documentos/B3N/datos/tablas_k), --gases H2O CO2 ..., --hilos N (por defecto, todos menos uno).
#
# Rejilla de la tabla (decision tecnica T15, PROVISIONAL; el error de interpolar en ella se medira, §4):
#   presion: 36 niveles de 1050 a 0,01 hPa, equiespaciados en log p;
#   temperatura: 150 a 330 K cada 30 K (7);
#   fraccion molar propia (solo el vapor, por el ensanchamiento consigo mismo): 0, 0,01 y 0,04; por encima de 100 hPa
#   solo se calcula la de 0 y se copia (alli el vapor es < 1e-4 y su ensanchamiento propio, despreciable);
#   bandas: las 16 del infrarrojo y las 14 de la luz de fase52_hitran (T8).
# Vapor: lineas sin el pedestal (T4) MAS el continuo MT_CKD 4.3 (fase52_continuo), en la misma tabla.
# Lineas descartadas (T14): umbral 1e-7 cm^-1 por linea con la columna maxima del gas (U_MAX); la suma de lo
# descartado en cada casilla queda en la tabla.

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
T_K = np.arange(150.0, 331.0, 30.0)
VMR_H2O = (0.0, 0.01, 0.04)
P_VMR_MIN = 100.0   # hPa
# columna maxima de cada gas (moleculas/cm^2) para el descarte: aire de 1050 hPa con g = 8 m/s^2 (cota por arriba de
# la columna: P3N tiene g = 9,05, parametros.P3N_GRAVEDAD) por una fraccion molar maxima generosa (mucho mayor que la sorteada, para que la tabla sirva
# tambien si se cambia la composicion por worldbuilding)
N_AIRE = 1050e2 / (28.96e-3 / 6.02214076e23 * 8.0) / 1e4
U_MAX = {"H2O": 0.05 * N_AIRE, "CO2": 2e-3 * N_AIRE, "O3": 2e-5 * N_AIRE, "N2O": 2e-6 * N_AIRE,
         "CH4": 5e-6 * N_AIRE, "O2": 0.3 * N_AIRE}
UMBRAL = 1e-7


def tareas(gases, bandas, prueba):
    out = []
    for g in gases:
        vmrs = VMR_H2O if g == "H2O" else (0.0,)
        for ib, banda in enumerate(bandas):
            for ip, p in enumerate(P_HPA):
                if prueba and (ib not in (0, 4, 10, 17) or ip not in (0, 18, 35)):   # una muestra de bandas y alturas
                    continue
                out.append((g, ib, banda, ip, float(p), vmrs))
    return out


_DATOS = {}


def _cargar(hitran, gas):
    if gas not in _DATOS:
        import fase52_hitran as H
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
        _DATOS[gas] = (H.leer_par(os.path.join(hitran, gas + ".data")), isotopos)
    return _DATOS[gas]


def calcular(args_tarea):
    hitran, salida, (gas, ib, banda, ip, p, vmrs) = args_tarea
    import fase52_generador as G
    import fase52_continuo as C
    fichero = os.path.join(salida, "casillas", f"{gas}_b{ib:02d}_p{ip:02d}.npz")
    if os.path.exists(fichero):
        return fichero, 0.0
    t0 = time.time()
    lineas, isotopos = _cargar(hitran, gas)
    ruta_mtckd = os.path.join(hitran, "absco-ref_wv-mt-ckd.nc")
    continuo = C.leer_mtckd(ruta_mtckd) if gas == "H2O" and os.path.exists(ruta_mtckd) else None
    k = np.zeros((len(T_K), len(vmrs), 16))
    perdida = np.zeros((len(T_K), len(vmrs)))
    usadas = np.zeros((len(T_K), len(vmrs)), dtype=int)
    for it, T in enumerate(T_K):
        for iv, v in enumerate(vmrs):
            if iv > 0 and p < P_VMR_MIN:            # copia de la de vmr 0 (ver la rejilla, arriba)
                k[it, iv], perdida[it, iv], usadas[it, iv] = k[it, 0], perdida[it, 0], usadas[it, 0]
                continue
            extra = None
            if gas == "H2O" and continuo is not None:                # el continuo MT_CKD se suma a las lineas
                extra = (lambda nu, T=float(T), p=p, v=v: C.k_continuo(nu, T, p, v, continuo))
            kg, per, n, _ = G.k_gas_banda(lineas, isotopos, banda, p / 1013.25, float(T), v, U_MAX[gas], UMBRAL,
                                          quitar_pedestal=(gas == "H2O"), k_extra=extra)
            k[it, iv], perdida[it, iv], usadas[it, iv] = kg, per, n
    np.savez(fichero + ".tmp.npz", k=k, perdida=perdida, usadas=usadas)
    os.replace(fichero + ".tmp.npz", fichero)
    return fichero, time.time() - t0


def main():
    import fase52_hitran as H
    ap = argparse.ArgumentParser(description="Tablas k de M3N (Fase 5.2) a partir de HITRAN")
    ap.add_argument("--hitran", default=os.path.expanduser("~/Documentos/B3N/datos/hitran"))
    ap.add_argument("--salida", default=os.path.expanduser("~/Documentos/B3N/datos/tablas_k"))
    ap.add_argument("--gases", nargs="+", default=list(GASES))
    ap.add_argument("--hilos", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--prueba", action="store_true", help="solo 4 bandas y 3 presiones, para medir el tiempo")
    a = ap.parse_args()
    for falta in ("isotopologos.json", "absco-ref_wv-mt-ckd.nc"):
        if "H2O" in a.gases or falta == "isotopologos.json":
            if not os.path.exists(os.path.join(a.hitran, falta)):
                sys.exit(f"Falta {falta}: vuelve a lanzar herramientas/descargar_hitran.py (no repite lo ya descargado).")
    bandas = list(H.BANDAS_IR) + [b for b in H.BANDAS_LUZ if b not in H.BANDAS_IR]
    os.makedirs(os.path.join(a.salida, "casillas"), exist_ok=True)
    lista = tareas(a.gases, bandas, a.prueba)
    hechas = sum(os.path.exists(os.path.join(a.salida, "casillas", f"{g}_b{ib:02d}_p{ip:02d}.npz"))
                 for g, ib, _, ip, _, _ in lista)
    print(f"{len(lista)} casillas ({hechas} ya hechas); {a.hilos} hilos", flush=True)
    import multiprocessing as mp
    t0, hechas_ahora, tiempos = time.time(), 0, []
    with mp.Pool(a.hilos) as pool:
        for f, dt in pool.imap_unordered(calcular, [(a.hitran, a.salida, t) for t in lista]):
            hechas += 1
            if dt > 0:
                hechas_ahora += 1
                tiempos.append(dt)
            transc = time.time() - t0
            quedan = (len(lista) - hechas) * transc / max(hechas_ahora, 1)
            b = "#" * int(25 * hechas / len(lista))
            sys.stdout.write(f"\r[{b:25s}] {hechas}/{len(lista)} casillas | {transc / 60:6.1f} min | quedan ~{quedan / 3600:5.1f} h ")
            sys.stdout.flush()
    print()
    if a.prueba:
        n_total = len(tareas(a.gases, bandas, False))
        media = float(np.mean(tiempos)) if tiempos else 0.0
        print(f"PRUEBA: {media:.1f} s por casilla de media con {a.hilos} hilos en paralelo; el total ({n_total} casillas)"
              f" llevaria del orden de {n_total * media / a.hilos / 3600:.1f} h (orientativo: las bandas difieren).")
        return
    # juntar
    for g in a.gases:
        vmrs = VMR_H2O if g == "H2O" else (0.0,)
        k = np.zeros((len(bandas), len(P_HPA), len(T_K), len(vmrs), 16), dtype=np.float32)
        perd = np.zeros((len(bandas), len(P_HPA), len(T_K), len(vmrs)))
        for ib in range(len(bandas)):
            for ip in range(len(P_HPA)):
                d = np.load(os.path.join(a.salida, "casillas", f"{g}_b{ib:02d}_p{ip:02d}.npz"))
                k[ib, ip], perd[ib, ip] = d["k"], d["perdida"]
        import fase52_k as K
        gq, wq = K.nodos_g()
        np.savez_compressed(os.path.join(a.salida, f"tabla_k_{g}.npz"), k=k, perdida_cm1=perd, bandas=np.array(bandas),
                            p_hPa=P_HPA, T_K=T_K, vmr=np.array(vmrs), g=gq, pesos=wq, umbral=UMBRAL, u_max=U_MAX[g])
        print(f"{g}: tabla_k_{g}.npz | absorcion descartada maxima en una casilla {perd.max():.2e} cm^-1", flush=True)


if __name__ == "__main__":
    main()

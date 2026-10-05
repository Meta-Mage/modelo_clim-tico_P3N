# calibrar_v30.py -- v3.0: calibracion de la difusion atmosferica (D_ATMOSFERA_V30)
# y oceanica (D_OCEANO) en MODO TIERRA, frente a Trenberth y Caron (2001).
#
# Uso (terminal B, venv de M3N):
#     M3N_MODO=tierra python calibrar_v30.py            # rejilla completa
#     M3N_MODO=tierra python calibrar_v30.py --rapido   # 1 año, solo para ver que funciona
#
# Que hace: simula la Tierra (Sol, orbita, eje, dia, gravedad y mapa
# terrestres; parametros.py y modo_tierra.py) con la atmosfera de N capas
# (I10) para varias combinaciones de D, hasta el equilibrio, y mide el
# transporte de calor hacia los polos de la atmosfera y del oceano.
#
# Referencia (Trenberth y Caron 2001, J. Climate 14, 3433, resumen ✅):
#   - transporte atmosferico maximo 5,0 +- 0,14 PW a 43 N, parecido
#     cerca de 40 S;
#   - a 35 grados, la atmosfera lleva el 78 % del total en el norte y el
#     92 % en el sur (oceano: 22 % y 8 %).
#
# Igual que en la calibracion de la v2.2c (DISENO_FASE2B.md, 12.3), sin
# hielo marino: el hielo se valida despues, con D ya fijado, y asi cada
# simulacion converge mucho antes. OJO: la v3.0 todavia no tiene vapor
# de agua (Fase 5a): el D que salga incluye "de prestado" el transporte
# de calor latente y se recalibra en la v3.1.
#
# Las simulaciones van en paralelo (una por nucleo, hasta --nucleos).

import os
import sys
import json
import time
import argparse

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")   # v3.1: un hilo por simulacion (van varias en paralelo)

import numpy as np
import parametros as P

D_ATM_REJILLA = [0.40, 0.55, 0.70, 0.90]
D_OC_REJILLA = [0.12, 0.20]
OBS = {"pico_atm_N": 5.0, "lat_pico_N": 43, "pico_atm_S": 5.0, "lat_pico_S": -40,
       "oceano_35N": 22, "oceano_35S": 8}
CARPETA = os.path.join("outputs", "calibracion_v30")


def una(d_atm, d_oc, max_anos, acelerar, n_capas=None, v31=False):
    import fase2b_atmosfera as F
    import fase30_multicapa as M30
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
    from cache_simulacion import precalcular_orbita_cacheada
    from rejilla import LATITUDES_GRADOS
    from modo_tierra import mapa_tierra
    tipo, alt = mapa_tierra()
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    I = dict(F.INTERRUPTORES_FASE2B)
    I["hielo_marino"] = False
    I["atmosfera_multicapa"] = True
    if v31:      # v3.1: con el ciclo del agua (el transporte incluye el calor latente)
        for k in ("ciclo_agua", "conveccion_humeda", "suelo_termico_agua"):
            I[k] = True
    t = time.time()
    r = F.simular_fase2b(orbita, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
                         P.PROFUNDIDAD_OPTICA, 0.55, interruptores=I, d_atmosfera=d_atm, d_oceano=d_oc,
                         max_anos=max_anos, acelerar=acelerar, n_capas_atm=n_capas)
    bordes, atm = M30.transporte_meridional(r["flujos"]["conv_atmosfera"], P.P3N_RADIO)
    _, oc = M30.transporte_meridional(r["flujos"]["conv_oceano"], P.P3N_RADIO)
    peso = np.cos(np.radians(LATITUDES_GRADOS))
    zonal = r["reg_media"].mean(axis=(0, 2))
    n, s = bordes > 0, bordes < 0
    i35n, i35s = int(np.argmin(np.abs(bordes - 35))), int(np.argmin(np.abs(bordes + 35)))
    tot = atm + oc
    return {
        "D_atm": d_atm, "D_oc": d_oc, "N": n_capas, "anos": int(r["anos"]),
        "zonal_aire2m": zonal.tolist(), "segundos": round(time.time() - t),
        "aire2m_global": float((zonal * peso).sum() / peso.sum()),
        "aire2m_ecuador": float(zonal[17:19].mean()),
        "aire2m_polo_N": float(zonal[0]), "aire2m_polo_S": float(zonal[-1]),
        "pico_atm_N": float(atm[n].max()), "lat_pico_N": float(bordes[n][np.argmax(atm[n])]),
        "pico_atm_S": float(-atm[s].min()), "lat_pico_S": float(bordes[s][np.argmin(atm[s])]),
        "oceano_35N": float(100 * oc[i35n] / tot[i35n]), "oceano_35S": float(100 * oc[i35s] / tot[i35s]),
        "total_35N": float(tot[i35n]), "total_35S": float(-tot[i35s]),
        # la difusion conserva la energia: lo recibido por TODO el planeta debe ser ~0
        "cierre_polo_sur_PW": float(abs(sum((r["flujos"][k] * peso[:, None]).sum() for k in ("conv_atmosfera", "conv_oceano"))
                                        * P.P3N_RADIO ** 2 * np.radians(5) ** 2 / 1e15)),
        "bordes": bordes.tolist(), "atm_PW": atm.tolist(), "oc_PW": oc.tolist(),
        "balance": float(r["energia"]["diferencia_relativa"]),
    }


N_PRUEBA = [10, 20, 30, 40]


def informe_n(resultados):
    """Convergencia en el numero de capas: se elige el N mas pequeño cuyo
    resultado ya no cambia al pasar al siguiente (mismo criterio que la
    prueba del paso de tiempo: < 0,1 C en la media global y < 0,3 C en
    cada banda de latitud)."""
    resultados.sort(key=lambda r: r["N"])
    nombre = os.path.join(CARPETA, "convergencia_n_v30.json")
    with open(nombre, "w") as f:
        json.dump(resultados, f, indent=1)
    print()
    print("  N | años | aire 2 m global | pico atm N / S (PW) | cambio respecto al N siguiente: global / peor banda")
    for i, r in enumerate(resultados):
        txt = ""
        if i + 1 < len(resultados):
            s = resultados[i + 1]
            dg = abs(s["aire2m_global"] - r["aire2m_global"])
            db = float(np.max(np.abs(np.array(s["zonal_aire2m"]) - np.array(r["zonal_aire2m"]))))
            ok = "OK" if dg < 0.1 and db < 0.3 else "todavia cambia"
            txt = f"{dg:.3f} / {db:.3f} C  ({ok})"
        print(f"{r['N']:3d} | {r['anos']:4d} | {r['aire2m_global']:8.2f} C | {r['pico_atm_N']:5.2f} / {r['pico_atm_S']:5.2f} | {txt}")
    print(f"Guardado en {nombre}")


def _trabajo(args):
    return una(*args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rapido", action="store_true", help="1 año sin acelerar, solo para comprobar")
    ap.add_argument("--nucleos", type=int, default=max(1, min(4, (os.cpu_count() or 2) // 2)))
    ap.add_argument("--max-anos", type=int, default=60)
    ap.add_argument("--d-atm", type=float, nargs="+", default=None,
                    help="valores de D de la atmosfera a probar (por defecto, la rejilla de la tanda 1)")
    ap.add_argument("--d-oc", type=float, nargs="+", default=None, help="valores de D del oceano a probar")
    ap.add_argument("--v31", action="store_true", help="con el ciclo del agua (I11-I13)")
    ap.add_argument("--convergencia-n", action="store_true",
                    help="en vez de la rejilla de D: misma simulacion con N = 10, 20, 30 y 40 capas")
    a = ap.parse_args()
    if not P.MODO_TIERRA:
        sys.exit("Hay que lanzarlo en modo Tierra:  M3N_MODO=tierra python calibrar_v30.py")
    os.makedirs(CARPETA, exist_ok=True)
    # la orbita se calcula (y se guarda en la cache) una vez aqui, antes de repartir el trabajo,
    # para que las simulaciones en paralelo no la escriban a la vez
    from cache_simulacion import precalcular_orbita_cacheada
    precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    if a.convergencia_n:
        trabajos = [(0.55, 0.12, a.max_anos, True, n) for n in N_PRUEBA]
    elif a.rapido:
        trabajos = [(0.55, 0.12, 1, False)]
    else:
        # v3.0 (05/10/2026): tanda 1 = la rejilla por defecto; tandas 2 y 3 con --d-atm 1.0 1.15 1.3
        # --d-oc 0.16 0.22 y --d-atm 1.35 1.5 --d-oc 0.22 0.28 (DISENO_V3.0.md, seccion 13)
        trabajos = [(da, do, a.max_anos, True, None, a.v31) for do in (a.d_oc or D_OC_REJILLA)
                    for da in (a.d_atm or D_ATM_REJILLA)]
    print(f"{len(trabajos)} simulaciones en modo Tierra, {a.nucleos} a la vez...", flush=True)
    resultados = []
    if a.nucleos > 1 and len(trabajos) > 1:
        from multiprocessing import Pool
        with Pool(a.nucleos) as pool:
            for r in pool.imap_unordered(_trabajo, trabajos):
                resultados.append(r)
                print(f"  hecha D_atm={r['D_atm']}, D_oc={r['D_oc']} ({r['anos']} años, {r['segundos']} s)", flush=True)
    else:
        for t in trabajos:
            r = _trabajo(t)
            resultados.append(r)
            print(f"  hecha D_atm={r['D_atm']}, D_oc={r['D_oc']} ({r['anos']} años, {r['segundos']} s)", flush=True)
    if a.convergencia_n:
        informe_n(resultados)
        return
    resultados.sort(key=lambda r: (r["D_oc"], r["D_atm"]))
    nombre = os.path.join(CARPETA, ("calibracion_v31" if a.v31 else "calibracion_v30") + ("_rapida" if a.rapido else "")
                          + "_" + time.strftime("%Y%m%d_%H%M%S") + ".json")   # v3.1: nunca se sobrescribe
    with open(nombre, "w") as f:
        json.dump({"observado": OBS, "resultados": resultados}, f, indent=1)

    print()
    print("D atm | D oc | años | pico atm N (lat) | pico atm S (lat) | % oceano 35N/35S | total 35N/S | aire 2 m global / ecuador / polos N,S")
    for r in resultados:
        print(f"{r['D_atm']:5.2f} | {r['D_oc']:4.2f} | {r['anos']:4d} | {r['pico_atm_N']:5.2f} PW ({r['lat_pico_N']:+.0f}) | "
              f"{r['pico_atm_S']:5.2f} PW ({r['lat_pico_S']:+.0f}) | {r['oceano_35N']:5.1f} / {r['oceano_35S']:5.1f} | "
              f"{r['total_35N']:4.2f} / {r['total_35S']:4.2f} | {r['aire2m_global']:5.1f} / {r['aire2m_ecuador']:5.1f} / "
              f"{r['aire2m_polo_N']:6.1f}, {r['aire2m_polo_S']:6.1f}")
    print(f"Observado (Trenberth y Caron 2001): pico atm 5,0 PW a 43 N, parecido a 40 S; % oceano a 35: 22 N / 8 S")
    print(f"Comprobacion de cierre (transporte en el polo sur, debe ser ~0): "
          f"max {max(r['cierre_polo_sur_PW'] for r in resultados):.2e} PW")
    print(f"Guardado en {nombre}")


if __name__ == "__main__":
    main()

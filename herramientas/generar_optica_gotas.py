# generar_optica_gotas.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7, T18): calcula la TABLA de las propiedades opticas
# de las gotas de agua de las nubes en las 16 bandas del INFRARROJO (fase52_optica_nubes: Mie propio, indices de T17,
# promedio de banda pesado con Planck a la temperatura de la nube). No depende de ningun dato descargado: se calculo en
# el entorno de la IA y el resultado esta en el repositorio (datos/optica_gotas_infrarrojo.npz); este programa queda
# para repetirlo. Las bandas de la luz se haran cuando este el espectro de S3N (el peso es el de la estrella).
#
# Uso:  python herramientas/generar_optica_gotas.py [--hilos N]           (~30 min con 1 hilo)
#
# Rejilla: radio efectivo 2, 3, 4, 6, 8, 10, 12, 15, 20, 25, 30 y 40 µm (varianza efectiva 0,1, la de
# fase52_optica_nubes); temperatura 240, 253, 263, 273 y 298 K (los nodos de T17). 20 numeros de onda por banda
# (medido: con 40 cambia < 0,01 %) y 400 radios por poblacion (con 200 u 800, igual en la quinta cifra).

import argparse
import os
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, AQUI)

R_E = np.array([2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0])
T_NODOS = np.array([240.0, 253.0, 263.0, 273.0, 298.0])
SALIDA = os.path.join(AQUI, "datos", "optica_gotas_infrarrojo.npz")


def casilla(args):
    import fase52_optica_nubes as O
    ib, a, b, ir, it = args
    return ib, ir, it, O.optica_banda_agua(a, b, float(R_E[ir]), float(T_NODOS[it]), n_nu=20, n_radios=400)


def main():
    import fase52_hitran as H
    ap = argparse.ArgumentParser(description="Tabla de la optica de las gotas en el infrarrojo (Fase 5.2)")
    ap.add_argument("--hilos", type=int, default=1)
    ap.add_argument("--salida", default=SALIDA)
    a = ap.parse_args()
    bandas = list(H.BANDAS_IR)
    tareas = [(ib, ba, bb, ir, it) for ib, (ba, bb) in enumerate(bandas) for ir in range(len(R_E))
              for it in range(len(T_NODOS))]
    out = np.zeros((3, len(bandas), len(R_E), len(T_NODOS)))
    t0 = time.time()
    if a.hilos > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.hilos) as pool:
            res = list(pool.imap_unordered(casilla, tareas))
    else:
        res = [casilla(t) for t in tareas]
    for ib, ir, it, v in res:
        out[:, ib, ir, it] = v
    os.makedirs(os.path.dirname(a.salida), exist_ok=True)
    with open(a.salida, "wb") as f:
        np.savez_compressed(f, beta_m2_kg=out[0], omega=out[1], g=out[2], r_e_um=R_E, T_K=T_NODOS,
                            bandas=np.array(bandas, float), v_e=np.array(0.1))
    print(f"Escrito {a.salida} en {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()

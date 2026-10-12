# control_tierra_p3n.py -- 3.20.0, diagnostico de control (Fase 6.4, 12/10/2026; resultado en DISENO_FASE_6.4.md §8.3
# y DISENO_FASE_6.3.md §6.26). NO cambia nada del modelo ni simula.
# Pregunta: la estructura de P3N (sin celula polar, oeste en superficie hasta el polo, oeste en la alta troposfera
# ecuatorial) es de P3N o del CODIGO? Aplica a la validacion del modo Tierra (3.18.1, mismo codigo que la 3.19.0
# con los valores por defecto) EXACTAMENTE las mismas funciones de analizar_p3n.py que se usaron con P3N.
#
# Uso, en la carpeta de M3N con el venv activado:
#   python herramientas/control_tierra_p3n.py outputs/validacion_i16_3.18.1/campos.npz outputs/analisis_p3n/campos.npz
import os
import sys
import numpy as np

sys.path.insert(0, os.getcwd())             # importar M3N desde la carpeta actual (el script vive en herramientas/)

f_tierra, f_p3n = sys.argv[1], sys.argv[2]
sys.argv = sys.argv[:1]                      # analizar_p3n lee argumentos al importarse: que no vea los nuestros

import analizar_p3n as A
from fase30_multicapa import sigma_seminiveles
from rejilla import FILAS, LATITUDES_GRADOS as LAT

sig = np.asarray(sigma_seminiveles(), dtype=float)          # el mismo que usa el nucleo (fase2b_atmosfera)
sm = 0.5 * (sig[1:] + sig[:-1])
k5 = int(np.argmin(np.abs(sig - 0.5)))
k25 = int(np.argmin(np.abs(sm - 0.25)))
lat_c = 90.0 - np.arange(1, FILAS) * (180.0 / FILAS)        # caras interiores, de norte a sur


def describir(nombre, f):
    c = dict(np.load(f))
    print(f"=== {nombre} ({f})")
    for clave in [k[4:] for k in c if k.startswith("psi_")]:
        psi, u = c["psi_" + clave], c.get("u_" + clave)
        cs = A.celulas(psi, lat_c, sig, k5)
        print(f"  celulas {clave}: {len(cs)}")
        for x in cs:
            print(f"    {x['extremo']:+7.2f} a {x['lat']:+5.1f} | de {x['norte']:+5.1f} a {x['sur']:+5.1f}"
                  + (" (debil)" if x["debil"] else ""))
        # una celula polar poco profunda no cambiaria el signo a sigma 0,5: se mira tambien a sigma ~0,8
        k8 = int(np.argmin(np.abs(sig - 0.8)))
        pol = np.abs(lat_c) >= 50
        print(f"  psi a sigma ~0,8 desde 50 grados hacia los polos ({clave}): "
              + " ".join(f"{la:+.0f}:{x:+.2f}" for la, x in zip(lat_c[pol], psi[k8, pol])))
        if u is None:
            continue
        print(f"  bandas de viento en la capa baja ({clave}):")
        for b in A.bandas_viento(u[-1], LAT):
            print(f"    {'oeste' if b['signo'] > 0 else 'este ':5s} {b['extremo']:+5.1f} m/s a {b['lat']:+5.1f}"
                  f" | de {b['norte']:+5.1f} a {b['sur']:+5.1f}")
        ec = np.abs(LAT) < 10
        print(f"  viento zonal a sigma ~0,25 entre 10 S y 10 N ({clave}): "
              + " ".join(f"{la:+.1f}:{x:+.1f}" for la, x in zip(LAT[ec], u[k25, ec])) + " m/s")
    print()


describir("MODO TIERRA 3.18.1", f_tierra)
describir("P3N 3.19.0", f_p3n)

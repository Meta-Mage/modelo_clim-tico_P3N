# herramientas/diagnostico_tierra.py -- DIAGNOSTICO (no forma parte del modelo; v3.10.0-v3.10.2, DISENO_FASE_6.3.md
# §6.11-6.12): repite el año 1 del modo Tierra con el
# nucleo dinamico desde el mismo arranque caliente que clima_dinamico.py y, si aparece un valor no valido,
# guarda lo necesario para encontrar la causa:
#   - outputs/diagnostico_tierra/ultimos_estados.npz : los ultimos 60 subpasos del nucleo (p_s, T, u, v) y la
#     tendencia de la fisica que se les aplico (float32)
#   - outputs/diagnostico_tierra/informe.txt          : donde y cuando se rompe
# Uso (DESDE la carpeta de M3N, con el venv activado):  python herramientas/diagnostico_tierra.py
import os, sys, time, pickle, warnings
os.environ["M3N_MODO"] = "tierra"
sys.path.insert(0, os.getcwd())          # los modulos de M3N, de la carpeta desde la que se lanza
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
warnings.simplefilter("ignore")
from collections import deque
import numpy as np
import parametros as P, fase2b_atmosfera as F, fase6_nucleo as N6, fase6_superficie as SUP
from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from cache_simulacion import precalcular_orbita_cacheada
from modo_tierra import mapa_tierra
from rejilla import LATITUDES_GRADOS as LAT, LONGITUDES_GRADOS as LON

base = os.getcwd()
ini = pickle.load(open(os.path.join(base, "outputs", "clima_dinamico", "tierra_dos_capas", "arranque_dos_capas.pkl"), "rb"))["estado"]
salida = os.path.join(base, "outputs", "diagnostico_tierra")
os.makedirs(salida, exist_ok=True)
orb = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
tipo, alt = mapa_tierra()
I = dict(F.INTERRUPTORES_FASE2B)
for k in ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral", "nucleo_dinamico"):
    I[k] = True

memoria = deque(maxlen=60)
c = {"n": 0, "t0": time.time()}
lineas = []


def guardar(motivo):
    est = list(memoria)
    np.savez_compressed(os.path.join(salida, "ultimos_estados.npz"),
                        sub=np.array([e[0] for e in est]),
                        ps=np.array([e[1] for e in est]), T=np.array([e[2] for e in est]),
                        u=np.array([e[3] for e in est]), v=np.array([e[4] for e in est]),
                        FT=np.array([e[5] for e in est]), Fu=np.array([e[6] for e in est]), Fv=np.array([e[7] for e in est]),
                        alt=alt, tipo=tipo)
    with open(os.path.join(salida, "informe.txt"), "w") as f:
        f.write(motivo + "\n" + "\n".join(lineas[-200:]) + "\n")
    print("\n" + motivo + f"\nGuardado en {salida} (informe.txt y ultimos_estados.npz): pasaselos a Claude.", flush=True)


_orig = N6.NucleoSeco.avanzar


def avanzar(self, ant, act, forz=None, **kw):
    f = forz(*ant) if forz is not None else None
    r = _orig(self, ant, act, forz, **kw)
    c["n"] += 1
    ps, T, u, v = r[1]
    z = lambda x: np.zeros_like(x, dtype=np.float32)
    memoria.append((c["n"], ps.astype(np.float32), T.astype(np.float32), u.astype(np.float32), v.astype(np.float32),
                    f[1].astype(np.float32) if f is not None else z(T), f[2].astype(np.float32) if f is not None else z(u),
                    f[3].astype(np.float32) if f is not None else z(v)))
    if c["n"] % 200 == 0:
        k, i, j = np.unravel_index(np.abs(u).argmax(), u.shape)
        lin = (f"subpaso {c['n']:6d} ({100 * c['n'] / (2 * len(orb)):5.1f} % del año) | max|u| {abs(u).max():6.1f} capa {k} "
               f"lat {LAT[i]:+.1f} lon {LON[j]:+.1f} | T {T.min():.0f}..{T.max():.0f} | ps {ps.min() / 100:.0f}..{ps.max() / 100:.0f}")
        lineas.append(lin)
        sys.stdout.write("\r" + lin + f" | {time.time() - c['t0']:.0f} s ")
        sys.stdout.flush()
    if not (np.isfinite(T).all() and np.isfinite(u).all() and np.isfinite(v).all() and np.isfinite(ps).all()):
        guardar(f"VALORES NO FINITOS EN EL NUCLEO en el subpaso {c['n']}")
        sys.exit(1)
    return r


_coef = SUP.coeficientes


def coeficientes(z_a, viento, th_a, th_s, *a, **k):
    d = _coef(z_a, viento, th_a, th_s, *a, **k)
    malo = ~np.isfinite(d["c_h"]) | ~np.isfinite(np.asarray(z_a)) | (np.asarray(z_a) <= 0)
    if malo.any():
        i, j = np.argwhere(malo)[0]
        guardar(f"COEFICIENTES DE SUPERFICIE NO VALIDOS tras el subpaso {c['n']} en lat {LAT[i]:+.1f} lon {LON[j]:+.1f} "
                f"(altitud {alt[i, j]:.0f} m): z_a {np.asarray(z_a)[i, j]}, viento {np.asarray(viento)[i, j]}, "
                f"T aire {np.asarray(th_a)[i, j]}, T superficie {np.asarray(th_s)[i, j]}")
        sys.exit(1)
    return d


N6.NucleoSeco.avanzar = avanzar
SUP.coeficientes = coeficientes
print("Diagnostico: año 1 del modo Tierra con el nucleo dinamico (se para solo si algo se rompe)...", flush=True)
F.simular_fase2b(orb, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55,
                 interruptores=I, max_anos=1, estado_inicial=ini)
print("\nEl año 1 ha terminado SIN romperse (no se reproduce el fallo).")

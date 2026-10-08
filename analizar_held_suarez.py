# analizar_held_suarez.py -- informe de la prueba de Held y Suarez (DISENO_FASE_6.2.md §10).
#
# Uso: python analizar_held_suarez.py      (lee outputs/held_suarez/hs_tau*_medias.npz)
# Imprime, para cada simulacion y cada hemisferio: el chorro en altura (velocidad, latitud, sigma), el
# viento del oeste en superficie, el maximo del transporte de calor y de momento por remolinos, y la
# incertidumbre (diferencia entre la primera y la segunda mitad del promedio). Guarda una figura por
# simulacion (viento zonal y temperatura medios, latitud x sigma) en outputs/held_suarez/.
import glob
import os

import numpy as np


def max_parabola(y, x):
    """Maximo de y(x) refinado con una parabola por los tres puntos de alrededor (latitud con decimales)."""
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        den = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return y1 - 0.25 * (y0 - y2) * d, x[i] + d * (x[i + 1] - x[i])
    return y[i], x[i]


def chorro(u, lat, sig, hemis):
    sel = lat > 0 if hemis == "N" else lat < 0
    la = lat[sel]; uu = u[:, sel]
    k = np.unravel_index(uu.argmax(), uu.shape)[0]
    v, l = max_parabola(uu[k], la)
    return v, abs(l), sig[k]


def main():
    carpeta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "held_suarez")
    archivos = sorted(glob.glob(os.path.join(carpeta, "hs_tau*_medias.npz")))
    if not archivos:
        print("No hay resultados: primero python held_suarez.py"); return
    for f in archivos:
        d = np.load(f)
        lat, sig = d["lat"], d["sigma"]
        res = f"{180 / int(d['filas']):g} grados, " if "filas" in d.files else "5 grados, "
        if "niveles" in d.files:                                                   # v3.4.0
            res += f"niveles {str(d['niveles'])}, "
        if int(d["n_medias"]) == 0 or np.ndim(d["u"]) != 2:                       # v3.4.0: corrida corta
            print(f"\n=== {os.path.basename(f)}: {res}{int(d['dias'])} dias simulados: sin dias promediados "
                  f"(la media empieza el dia 200); no hay informe ===")
            continue
        print(f"\n=== {os.path.basename(f)}: {res}tau = {float(d['tau']):g} dias, {int(d['dias'])} dias simulados, "
              f"media de {int(d['n_medias'])} dias ===")
        ks = int(np.argmax(sig))                                                   # capa mas baja
        k25 = int(np.argmin(abs(sig - 0.25))); k85 = int(np.argmin(abs(sig - 0.85)))
        for h in ("N", "S"):
            v, l, s = chorro(d["u"], lat, sig, h)
            if d["u1"].ndim == 2 and d["u2"].ndim == 2:
                v1, l1, _ = chorro(d["u1"], lat, sig, h)
                v2, l2, _ = chorro(d["u2"], lat, sig, h)
            else:                                   # simulacion demasiado corta para partirla en dos mitades
                v1 = l1 = v2 = l2 = float("nan")
            sel = lat > 0 if h == "N" else lat < 0
            us, ls = max_parabola(d["u"][ks, sel], lat[sel])
            sgn = 1 if h == "N" else -1
            vt, lvt = max_parabola(sgn * d["vT"][k85, sel], lat[sel])
            uv, luv = max_parabola(sgn * d["uv"][k25, sel], lat[sel])
            print(f" Hemisferio {h}: chorro {v:5.1f} m/s en {l:5.1f} grados, sigma {s:.2f} "
                  f"(1.a mitad {v1:.1f} m/s en {l1:.1f}; 2.a mitad {v2:.1f} m/s en {l2:.1f})")
            print(f"   viento del oeste en superficie: max {us:4.1f} m/s en {abs(ls):5.1f} | "
                  f"calor por remolinos (sigma 0,85): max en {abs(lvt):5.1f} | "
                  f"momento por remolinos hacia el polo (sigma 0,25): max {uv:5.1f} m2/s2 en {abs(luv):5.1f}")
        print(" Referencia (HS94): chorros hacia ~45 grados y ~250 hPa (verificado en Mendonca et al. 2016, THOR); "
              "la velocidad, ~30 m/s, sigue pendiente de cita textual.")
        alto = sig < 0.02                                                          # v3.4.0: vigilancia del techo
        if alto.any() and d["u1"].ndim == 2 and d["u2"].ndim == 2:
            ua1 = float(np.abs(d["u1"][alto]).max()); ua2 = float(np.abs(d["u2"][alto]).max())
            print(f" Capas altas (centro < 20 hPa): |u zonal medio| max {float(np.abs(d['u'][alto]).max()):.1f} m/s "
                  f"(1.a mitad {ua1:.1f}, 2.a mitad {ua2:.1f}; cambio {100 * (ua2 - ua1) / max(ua1, 1e-9):+.0f} %)")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
            c = ax[0].contourf(lat, sig, d["u"], levels=np.arange(-20, 45, 5), cmap="RdBu_r", extend="both")
            ax[0].contour(lat, sig, d["u"], levels=[0], colors="k", linewidths=0.8)
            fig.colorbar(c, ax=ax[0], label="u zonal medio (m/s)")
            c2 = ax[1].contourf(lat, sig, d["T"], levels=np.arange(180, 320, 10), cmap="viridis")
            fig.colorbar(c2, ax=ax[1], label="T zonal media (K)")
            for a_ in ax:
                a_.set_xlabel("latitud"); a_.invert_yaxis()
            ax[0].set_ylabel("sigma = p / p_s")
            fig.suptitle(f"Held y Suarez con el nucleo de M3N ({res}tau = {float(d['tau']):g} d, media de {int(d['n_medias'])} dias)")
            png = f.replace("_medias.npz", ".png")
            fig.savefig(png, dpi=110, bbox_inches="tight"); plt.close(fig)
            print(f" Figura: {png}")
        except Exception as e:                                                        # sin matplotlib: solo texto
            print(f" (sin figura: {e})")


if __name__ == "__main__":
    main()

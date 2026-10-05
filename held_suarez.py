# held_suarez.py -- Fase 6.2, prueba 6: Held y Suarez (1994) con el nucleo dinamico propio de M3N.
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python held_suarez.py                 # tau = 0,5 dias (el de CESM), 1200 dias
#     python held_suarez.py --tau 2         # sensibilidad: hiperdifusion 4 veces mas debil
#     python held_suarez.py --dias 300      # mas corto (solo para probar)
#
# Configuracion (DISENO_FASE6_2.md §5 y §10): Tierra (Williamson/HS94: a = 6,37122e6 m, Omega = 7,292e-5,
# g = 9,80616, R = 287, cp = 1004), rejilla de M3N (72 x 36, 5 grados), 20 capas sigma IGUALES (la
# especificacion de la prueba), dt = 450 s, semiimplicito + RAW, hiperdifusion nabla^4 implicita.
#
# - Muestra el PROGRESO en pantalla (barra, dias, velocidad, tiempo restante) y una linea cada 20 dias.
# - Guarda un PUNTO DE CONTROL cada 20 dias en outputs/held_suarez/: si se corta (Ctrl+C, apagon...),
#   se vuelve a lanzar la MISMA orden y continua exactamente donde iba.
# - Acumula la media de cada dia desde el dia 200 (u, T, transportes por remolinos) y al terminar escribe
#   outputs/held_suarez/hs_tauX_medias.npz; el informe lo da analizar_held_suarez.py.

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # con varios hilos, el algebra lineal va MAS lenta aqui
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import math
import sys
import time

import numpy as np

from fase6_nucleo import NucleoSeco
from fase6_forzamientos import HeldSuarez, sigma_capas
from fase6_aguas_someras import A_TIERRA, OMEGA_TIERRA, G_TIERRA

DT = 450.0
DIA_INICIO_MEDIA = 200
CADA_LINEA = 20            # dias entre lineas del registro y puntos de control


def main():
    ap = argparse.ArgumentParser(description="Prueba de Held y Suarez (1994) con el nucleo de la Fase 6")
    ap.add_argument("--tau", type=float, default=0.5, help="hiperdifusion: dias de amortiguamiento de la onda mas corta")
    ap.add_argument("--dias", type=int, default=1200)
    a = ap.parse_args()

    carpeta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "held_suarez")
    os.makedirs(carpeta, exist_ok=True)
    nombre = f"hs_tau{a.tau:g}"
    f_estado = os.path.join(carpeta, nombre + "_estado.npz")
    f_medias = os.path.join(carpeta, nombre + "_medias.npz")
    f_log = os.path.join(carpeta, nombre + ".log")

    m = NucleoSeco(A_TIERRA, OMEGA_TIERRA, G_TIERRA, 287.0, 1004.0, np.linspace(0, 1, 21))
    hs = HeldSuarez(m)
    N, Rj = m.N, m.Rj
    sig = sigma_capas(m)
    m.preparar_semiimplicito(DT)
    m.preparar_hiperdifusion(a.tau)
    pasos_dia = int(round(86400 / DT))
    total = a.dias * pasos_dia
    claves = ("u", "T", "uv", "vT", "TT", "uu", "vv", "u1")    # u1: media de la 1.a mitad (incertidumbre)

    if os.path.exists(f_estado):
        d = np.load(f_estado)
        n0 = int(d["n"])
        ant = tuple(d["ant_" + k] for k in ("ps", "T", "u", "v"))
        act = tuple(d["act_" + k] for k in ("ps", "T", "u", "v"))
        acum = {k: d["acum_" + k] for k in claves}
        n_med = int(d["n_medias"]); n1 = int(d["n_mitad1"])
        reanudar = (n0, ant, act)
        print(f"-> continua desde el dia {n0 / pasos_dia:.0f} (punto de control {f_estado})")
    else:
        rng = np.random.default_rng(1)
        ps = np.full((Rj.filas, Rj.columnas), 1e5)
        T = hs.T_equilibrio(ps) + 0.1 * rng.standard_normal((N, Rj.filas, Rj.columnas))
        acum = {k: 0.0 for k in claves}
        n_med = 0; n1 = 0
        reanudar = None
        n0 = 1
    print(f"Held y Suarez | tau = {a.tau:g} dias | {a.dias} dias = {total} pasos de {DT:.0f} s | "
          f"media desde el dia {DIA_INICIO_MEDIA}")

    t0 = time.time()
    estado = {"n_med": n_med, "n1": n1}
    mitad = 0.5 * (DIA_INICIO_MEDIA + a.dias)

    def al_dia(n, ps, T, u, v):
        dia = n / pasos_dia
        uc = 0.5 * (u + np.roll(u, -1, axis=-1))                                  # u en el centro
        vc = np.zeros_like(uc); vc[:, :-1] += 0.5 * v; vc[:, 1:] += 0.5 * v       # v en el centro
        if dia > DIA_INICIO_MEDIA:
            uz = uc.mean(-1); vz = vc.mean(-1); Tz = T.mean(-1)
            up = uc - uz[..., None]; vp = vc - vz[..., None]; Tp = T - Tz[..., None]
            acum["u"] = acum["u"] + uz; acum["T"] = acum["T"] + Tz
            acum["uv"] = acum["uv"] + (up * vp).mean(-1); acum["vT"] = acum["vT"] + (vp * Tp).mean(-1)
            acum["TT"] = acum["TT"] + (Tp * Tp).mean(-1)
            acum["uu"] = acum["uu"] + (up * up).mean(-1); acum["vv"] = acum["vv"] + (vp * vp).mean(-1)
            estado["n_med"] += 1
            if dia <= mitad:
                acum["u1"] = acum["u1"] + uz; estado["n1"] += 1
        # progreso
        ahora = time.time()
        hecho = (n - 1) / (total - 1)
        vel = (n - n0) / max(ahora - t0, 1e-9)                                    # pasos por segundo
        resta = (total - n) / max(vel, 1e-9)
        barra = "#" * int(30 * hecho) + "-" * (30 - int(30 * hecho))
        sys.stdout.write(f"\r[{barra}] {100 * hecho:5.1f} %  dia {dia:6.0f}/{a.dias}  "
                         f"{vel * DT / 86400 * 60:5.1f} dias/min  quedan {resta / 60:5.0f} min ")
        sys.stdout.flush()
        if round(dia) % CADA_LINEA == 0:
            uz = uc.mean(-1)
            k, i = np.unravel_index(uz.argmax(), uz.shape)
            linea = (f"{dia:6.0f} d  u_zonal max {uz.max():5.1f} m/s en {math.degrees(Rj.phi_c[i]):6.1f}, "
                     f"sigma {sig[k]:.2f} | max|v| {abs(v).max():5.1f} | T {T.min():.0f}..{T.max():.0f} K | "
                     f"ps {ps.min() / 100:.0f}..{ps.max() / 100:.0f} hPa")
            sys.stdout.write("\r" + " " * 100 + "\r" + linea + "\n")
            with open(f_log, "a") as fl:
                fl.write(linea + "\n")
            if not np.isfinite(T).all():
                print("ERROR: la simulacion ha divergido (valores no finitos). Se detiene.")
                sys.exit(2)

    def guardar(n, ant, act):
        datos = {"n": n, "n_medias": estado["n_med"], "n_mitad1": estado["n1"]}
        for nombre_k, x in zip(("ps", "T", "u", "v"), ant):
            datos["ant_" + nombre_k] = x
        for nombre_k, x in zip(("ps", "T", "u", "v"), act):
            datos["act_" + nombre_k] = x
        for k in claves:
            datos["acum_" + k] = np.asarray(acum[k])
        tmp = f_estado + ".tmp.npz"
        np.savez(tmp, **datos)
        os.replace(tmp, f_estado)

    if reanudar is None:
        z = (np.zeros((N, Rj.filas, Rj.columnas)), np.zeros((N, Rj.filas - 1, Rj.columnas)))
        m.integrar_si(ps, T, *z, total, forzamiento=hs, cada=pasos_dia, al_registrar=al_dia,
                      guardar=guardar, cada_guardar=CADA_LINEA * pasos_dia)
    else:
        m.integrar_si(None, None, None, None, total, forzamiento=hs, cada=pasos_dia, al_registrar=al_dia,
                      reanudar=reanudar, guardar=guardar, cada_guardar=CADA_LINEA * pasos_dia)

    nm = max(estado["n_med"], 1)
    med = {k: np.asarray(acum[k]) / nm for k in claves}
    med["u1"] = np.asarray(acum["u1"]) / max(estado["n1"], 1)
    n2 = estado["n_med"] - estado["n1"]
    med["u2"] = (np.asarray(acum["u"]) - np.asarray(acum["u1"])) / max(n2, 1)
    np.savez(f_medias, **med, n_medias=estado["n_med"],
             sigma=sig, lat=np.degrees(Rj.phi_c), tau=a.tau, dias=a.dias)
    print(f"\nTerminado en {(time.time() - t0) / 60:.0f} min. Medias de {estado['n_med']} dias en {f_medias}")
    print("Siguiente paso: python analizar_held_suarez.py")


if __name__ == "__main__":
    main()

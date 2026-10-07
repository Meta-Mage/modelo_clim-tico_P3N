# clima_dinamico.py -- v3.1-pre12: simulacion LARGA de M3N con el nucleo dinamico (I16), de principio a fin:
# arranque caliente -> equilibrio -> climatologia de N años (fase6_clima.py; DISENO_FASE6_3.md §6.10).
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python clima_dinamico.py --tierra        # modo Tierra con el mapa de la Tierra (validacion)
#     python clima_dinamico.py                 # P3N con el mapa activo de C3N
# Opciones:
#     --arranque v31    arranque caliente desde la v3.1 SIN nucleo en vez del modelo de 2 capas (alternativa,
#                       solo si la opcion por defecto se vuelve inestable)
#     --dias N          SOLO PARA PROBAR que todo funciona: "años" de N dias (no da un clima)
#     --carpeta RUTA    donde va todo (por defecto outputs/clima_dinamico/<tierra|p3n>_<arranque>)
#
# - Barra de progreso: fase, año, % del año, velocidad y tiempo restante del año.
# - Una linea por año en pantalla y en <carpeta>/registro.txt.
# - PUNTO DE CONTROL: si se corta (Ctrl+C, apagon...), la MISMA orden continua donde iba (el resultado es el
#   mismo bit a bit). Lo mas que se pierde es el año en curso.
# - Al terminar: <carpeta>/resumen.txt (lo que hay que pasarle a Claude) y <carpeta>/resultado.pkl.

import argparse
import os
import sys

ap = argparse.ArgumentParser(description="M3N con nucleo dinamico (I16): equilibrio y climatologia")
ap.add_argument("--tierra", action="store_true", help="modo Tierra (mapa y parametros terrestres)")
ap.add_argument("--arranque", choices=("dos_capas", "v31"), default="dos_capas")
ap.add_argument("--dias", type=int, default=None, help="solo para probar: años de N dias")
ap.add_argument("--carpeta", default=None)
args = ap.parse_args()
if args.tierra:
    os.environ["M3N_MODO"] = "tierra"            # antes de importar parametros
elif os.environ.get("M3N_MODO", "").strip().lower() == "tierra":
    sys.exit("M3N_MODO=tierra esta puesto en el entorno: usa --tierra (o quitalo) para que quede claro")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import pickle
import time
import numpy as np

import parametros as P
import fase2b_atmosfera as F
import fase6_clima as C
from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from cache_simulacion import precalcular_orbita_cacheada
from rejilla import LATITUDES_GRADOS


def main():
    assert P.MODO_TIERRA == args.tierra
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    from temperatura import PASO_TIEMPO
    pasos_dia = round(P.ROTACION_PERIODO / PASO_TIEMPO)
    if args.dias:
        orbita = orbita[:args.dias * pasos_dia]
    if args.tierra:
        from modo_tierra import mapa_tierra
        tipo, alt = mapa_tierra()
        nombre = "tierra"
    else:
        from puente_c3n import cargar_mapa_activo_de_c3n
        tipo, alt, nombre_mapa = cargar_mapa_activo_de_c3n()
        nombre = "p3n"
    carpeta = args.carpeta or os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "clima_dinamico",
                                           f"{nombre}_{args.arranque}" + (f"_prueba{args.dias}d" if args.dias else ""))
    os.makedirs(carpeta, exist_ok=True)
    f_reg = os.path.join(carpeta, "registro.txt")

    def informar(texto):
        sys.stdout.write("\r" + " " * 110 + "\r" + texto + "\n")
        sys.stdout.flush()
        with open(f_reg, "a", encoding="utf-8") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + texto + "\n")

    argumentos = (orbita, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55)
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral",
              "nucleo_dinamico"):
        I[k] = True

    t0 = time.time()
    barra = {"t": 0.0, "t_ano": time.time(), "ano": None}

    def al_paso(ano, p, total, registrado):
        ahora = time.time()
        if barra["ano"] != (ano, registrado):
            barra.update(ano=(ano, registrado), t_ano=ahora)
        if ahora - barra["t"] < 1.0 and p + 1 < total:
            return
        barra["t"] = ahora
        hecho = (p + 1) / total
        vel = (p + 1) / max(ahora - barra["t_ano"], 1e-9)
        resta = (total - p - 1) / max(vel, 1e-9)
        b = "#" * int(25 * hecho) + "-" * (25 - int(25 * hecho))
        fase = "climatologia" if registrado else "equilibrio  "
        sys.stdout.write(f"\r{fase} año {ano:3d} [{b}] {100 * hecho:5.1f} %  {vel:5.1f} pasos/s  "
                         f"quedan {resta / 60:5.1f} min del año | total {(ahora - t0) / 3600:5.2f} h ")
        sys.stdout.flush()

    def al_acabar_ano(ano, info):
        ok, v = info.get("equilibrio", (False, {}))
        s = info.get("serie_equilibrio")
        if not np.isfinite(info["T_atm"]).all():
            raise FloatingPointError(f"temperaturas no finitas al acabar el año {ano}")
        if not s:          # el modelo de partida del arranque caliente (sin nucleo)
            informar(f"Año {ano:3d} (modelo de partida) | cambio maximo de una celda {info['cambio']:.4f} K")
            return
        linea = f"Año {ano:3d} (equilibrio)"
        if s:
            linea += (f" | N {s['N'][-1]:+.3f} W/m2 | aire 2 m {s['T2m'][-1] - 273.15:6.2f} C"
                      + (f" | hielo {100 * s['hielo'][-1]:.2f} % oceano" if s["hielo"][-1] is not None else "")
                      + (f" | agua suelo {100 * s['agua_suelo'][-1]:.1f} %" if s["agua_suelo"][-1] is not None else ""))
        if v:
            linea += " | ventana: " + ", ".join(f"{k} {x[0]:+.4f}{'' if x[2] else ' (NO)'}" for k, x in v.items())
        informar(linea + (" -> EQUILIBRIO" if ok else ""))

    F.AL_PASO = al_paso
    F.AL_ACABAR_ANO = al_acabar_ano
    informar(f"=== M3N v3.1-pre14 con nucleo dinamico | {'modo Tierra' if args.tierra else 'P3N, mapa ' + nombre_mapa} | "
             f"arranque {args.arranque} | {len(orbita)} pasos por año | carpeta {carpeta}")
    try:
        prueba = {} if not args.dias else dict(          # prueba: pocos años y sin exigir equilibrio
            max_anos_equilibrio=2, min_anos=2, max_anos=2, opciones_dos_capas={"max_anos": 2},
            exigir_equilibrio_dos_capas=False)
        r = C.simular_clima(argumentos, I, carpeta, informar=informar, arranque=args.arranque, **prueba)
    except KeyboardInterrupt:
        informar("Interrumpido. Para seguir, lanza la MISMA orden: continua desde el ultimo año guardado.")
        sys.exit(1)
    except (AssertionError, FloatingPointError, IndexError, ValueError) as e:      # valores no validos
        informar(f"ERROR: la simulacion se ha vuelto inestable ({e}). Pasale a Claude {f_reg}.")
        sys.exit(2)
    with open(os.path.join(carpeta, "resultado.pkl"), "wb") as f:
        pickle.dump(r, f)
    escribir_resumen(r, carpeta, time.time() - t0, len(orbita))
    informar(f"Terminado. Resumen en {os.path.join(carpeta, 'resumen.txt')}")


def escribir_resumen(r, carpeta, segundos, pasos_ano):
    c = r["climatologia"]
    eq = c["equilibrio"]
    peso = np.cos(np.radians(LATITUDES_GRADOS))
    T = c["bandas_aire2m_C"].mean(axis=0)
    Pz = c["bandas_precipitacion"].mean(axis=0)
    L = ["RESUMEN de clima_dinamico.py (M3N v3.1-pre12, nucleo dinamico I16)",
         f"modo: {'Tierra' if args.tierra else 'P3N'} | arranque: {args.arranque} | pasos por año: {pasos_ano}"
         + (f" | PRUEBA con años de {args.dias} dias" if args.dias else ""),
         f"tiempo de esta ejecucion: {segundos / 3600:.2f} h",
         "",
         f"EQUILIBRIO: convergido {eq['convergido']} tras {eq['anos']} años"]
    for k, (val, umbral, cumple) in eq["valores"].items():
        L.append(f"  {k:11s} {val:+.5f} (umbral {umbral}) {'ok' if cumple else 'NO'}")
    s = eq["serie"]
    L.append("  serie anual N (W/m2):    " + " ".join(f"{x:+.3f}" for x in s["N"]))
    L.append("  serie anual aire 2 m (C): " + " ".join(f"{x - 273.15:.3f}" for x in s["T2m"]))
    L += ["", f"CLIMATOLOGIA: {c['anos_promediados']} años, cumple el criterio: {c['cumple_criterio']}",
          f"  error maximo aire 2 m {c['error_T_max_K']:.3f} K | precipitacion {c['error_P_max_mm_dia']:.3f} mm/dia "
          f"({100 * c['error_P_rel_max']:.1f} %) | bandas sin cumplir: T {c['bandas_T_sin_cumplir']}, P {c['bandas_P_sin_cumplir']}",
          f"  aire 2 m medio global {float((T * peso).sum() / peso.sum()):.3f} C | precipitacion media global "
          f"{float((Pz * peso).sum() / peso.sum()):.3f} mm/dia", "",
          "POR AÑO REGISTRADO (energia y agua):"]
    for i, d in enumerate(c["por_ano"]):
        e = d["energia"]
        L.append(f"  {i + 1:2d}: cierre energia {e.get('cierre_relativo', float('nan')):+.2e} | error dinamica "
                 f"{e.get('error_dinamica_W_m2', float('nan')):+.4f} W/m2 | resto {e.get('cierre_sin_dinamica', float('nan')):+.1e} | "
                 f"agua atm {d['cierre_agua_atmosfera_kg_m2']:+.1e} tierra {d['cierre_agua_tierra_kg_m2']:+.1e}")
    L += ["", "BANDAS (latitud, aire 2 m C, precipitacion mm/dia):"]
    for lat, t, p in zip(LATITUDES_GRADOS, T, Pz):
        L.append(f"  {lat:+6.1f} {t:8.2f} {p:8.3f}")
    with open(os.path.join(carpeta, "resumen.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

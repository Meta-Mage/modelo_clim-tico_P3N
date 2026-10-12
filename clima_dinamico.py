# clima_dinamico.py -- v3.10.0: simulacion LARGA de M3N con el nucleo dinamico (I16), de principio a fin:
# arranque caliente -> equilibrio -> climatologia de N años (fase6_clima.py; DISENO_FASE_6.3.md §6.10).
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python clima_dinamico.py --tierra        # modo Tierra con el mapa de la Tierra (validacion)
#     python clima_dinamico.py                 # P3N con el mapa activo de C3N
# Opciones:
#     --arranque v31    arranque caliente desde la v3.1 SIN nucleo en vez del modelo de 2 capas (alternativa,
#                       solo si la opcion por defecto se vuelve inestable)
#     --dias N          SOLO PARA PROBAR que todo funciona: "años" de N dias (no da un clima)
#     --anos N          3.15.0: climatologia EXPLORATORIA de N años registrados (N >= 12; por defecto 30, la
#                       normal de la OMM, que es la que vale para el clima oficial). Va a su propia carpeta.
#     --carpeta RUTA    donde va todo (por defecto outputs/clima_dinamico/<tierra|p3n>_<arranque>)
#     --oceano-difusion F, --oceano-capa-mezcla F
#                       4.0.0 (Fase 6.4, DISENO_FASE_6.4.md §4.2): PRUEBAS DE SENSIBILIDAD del oceano, con la
#                       difusion del oceano (D_oc) o la profundidad de la capa de mezcla multiplicadas por F (0,5 o
#                       2). Solo con --anos (son climatologias exploratorias, no el clima oficial). Van a su propia
#                       carpeta, y el modelo de partida del arranque usa la misma capa de mezcla (su difusion, no).
#
#     --continuar-desde CARPETA
#                       3.21.0: una climatologia NUEVA que empieza en el estado final de otra ya terminada (en
#                       equilibrio y sin deriva) del mismo mapa y la misma configuracion, sin repetir el arranque ni
#                       el equilibrio. Se copian su punto de control y su arranque; queda escrito en origen.json, en
#                       el resumen y en la exportacion. Para rehacer la climatologia con diagnosticos nuevos (las
#                       medias zonales por nivel de la 3.21.0) sin repetir los años de equilibrio.
#
# - Barra de progreso: fase, año, % del año, velocidad y tiempo restante del año.
# - Una linea por año en pantalla y en <carpeta>/registro.txt.
# - PUNTO DE CONTROL: si se corta (Ctrl+C, apagon...), la MISMA orden continua donde iba (el resultado es el
#   mismo bit a bit). Lo mas que se pierde es el año en curso.
# - Al terminar: <carpeta>/resumen.txt (lo que hay que pasarle a Claude) y <carpeta>/resultado.pkl.
# - 3.14.0: TRAZABILIDAD. Cada vez que empieza o continua, apunta en <carpeta>/procedencia.json la version y la
#   huella del codigo; si el codigo ha cambiado desde la vez anterior, avisa. El resumen lo recoge.

import argparse
import os
import sys

ap = argparse.ArgumentParser(description="M3N con nucleo dinamico (I16): equilibrio y climatologia")
ap.add_argument("--tierra", action="store_true", help="modo Tierra (mapa y parametros terrestres)")
ap.add_argument("--arranque", choices=("dos_capas", "v31"), default="dos_capas")
ap.add_argument("--dias", type=int, default=None, help="solo para probar: años de N dias")
ap.add_argument("--anos", type=int, default=None, help="climatologia exploratoria de N años (N >= 12)")
ap.add_argument("--carpeta", default=None)
ap.add_argument("--oceano-difusion", type=float, default=1.0, help="sensibilidad: factor de D_oc (con --anos)")
ap.add_argument("--oceano-capa-mezcla", type=float, default=1.0, help="sensibilidad: factor de la capa de mezcla (con --anos)")
ap.add_argument("--continuar-desde", default=None, help="3.21.0: carpeta de una climatologia terminada de la que partir")
args = ap.parse_args()
SENSIBILIDAD = {"oceano_difusion": args.oceano_difusion, "oceano_capa_mezcla": args.oceano_capa_mezcla}
ES_SENSIBILIDAD = any(v != 1.0 for v in SENSIBILIDAD.values())
if any(not (v > 0) for v in SENSIBILIDAD.values()):
    sys.exit("--oceano-difusion y --oceano-capa-mezcla tienen que ser mayores que 0")
if ES_SENSIBILIDAD and not (args.anos or args.dias):
    sys.exit("una prueba de sensibilidad del oceano es una climatologia exploratoria: añade --anos 12 "
             "(DISENO_FASE_6.4.md §4.2)")
if args.anos is not None and args.anos < 12 and not args.dias:
    sys.exit("--anos: una climatologia exploratoria necesita al menos 12 años (WMO-No. 1203 §5.2.3)")
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

VERSION = "3.21.0"   # version de M3N (SemVer, NOMENCLATURA.md): la de README.md y la primera de CHANGELOG.md


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
                                           f"{nombre}_{args.arranque}" + (f"_prueba{args.dias}d" if args.dias else "")
                                           + (f"_exploratoria{args.anos}" if args.anos and args.anos != C.ANOS_CLIMATOLOGIA
                                              and not args.dias else "")
                                           + (f"_oceanoD{args.oceano_difusion:g}" if args.oceano_difusion != 1.0 else "")
                                           + (f"_capa{args.oceano_capa_mezcla:g}" if args.oceano_capa_mezcla != 1.0 else ""))
    os.makedirs(carpeta, exist_ok=True)
    f_reg = os.path.join(carpeta, "registro.txt")

    def informar(texto):
        sys.stdout.write("\r" + " " * 110 + "\r" + texto + "\n")
        sys.stdout.flush()
        with open(f_reg, "a", encoding="utf-8") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + texto + "\n")

    inercias = INERCIA_POR_TIPO
    if args.oceano_capa_mezcla != 1.0:
        # la capacidad de la capa de mezcla es INERCIA_POR_TIPO[AGUA] * raiz(86400 s / pi) (fase2_inercia_multicapa):
        # multiplicar la inercia del agua por F multiplica por F su profundidad (48,7 m por defecto)
        from fase1_geografia import AGUA
        inercias = dict(INERCIA_POR_TIPO)
        inercias[AGUA] = INERCIA_POR_TIPO[AGUA] * args.oceano_capa_mezcla
    argumentos = (orbita, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, inercias, P.PROFUNDIDAD_OPTICA, 0.55)
    opciones_nucleo = None
    if args.oceano_difusion != 1.0:
        opciones_nucleo = {"d_oceano": F.D_OCEANO_V30 * args.oceano_difusion}
    # no se puede continuar una simulacion con otro oceano: una carpeta con simulacion y sin sensibilidad.json es
    # del modelo tal cual (factores 1 y 1)
    import json
    f_sens = os.path.join(carpeta, "sensibilidad.json")
    guardada = None
    if os.path.exists(f_sens):
        with open(f_sens, encoding="utf-8") as f:
            guardada = json.load(f)
    elif any(os.path.exists(os.path.join(carpeta, n)) for n in ("estado_i16.pkl", "climatologia.pkl", "resultado.pkl")):
        guardada = {k: 1.0 for k in SENSIBILIDAD}
    if guardada is not None and guardada != SENSIBILIDAD:
        sys.exit(f"{carpeta} es de otra configuracion del oceano ({guardada}); usa otra carpeta")
    if guardada is None and ES_SENSIBILIDAD:
        with open(f_sens, "w", encoding="utf-8") as f:
            json.dump(SENSIBILIDAD, f)
    origen = preparar_continuacion(args.continuar_desde, carpeta, informar) if args.continuar_desde else \
        leer_origen(carpeta)
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
    informar(f"=== M3N {VERSION} con nucleo dinamico | {'modo Tierra' if args.tierra else 'P3N, mapa ' + nombre_mapa} | "
             f"arranque {args.arranque} | {len(orbita)} pasos por año | carpeta {carpeta}"
             + (f" | SENSIBILIDAD del oceano: D_oc x{args.oceano_difusion:g}, capa de mezcla x{args.oceano_capa_mezcla:g}"
                if ES_SENSIBILIDAD else ""))
    try:
        prueba = {} if not args.dias else dict(          # prueba: pocos años y sin exigir equilibrio
            max_anos_equilibrio=2, anos=2, opciones_dos_capas={"max_anos": 2},
            exigir_equilibrio_dos_capas=False)
        if not args.dias:
            prueba = {"anos": args.anos or C.ANOS_CLIMATOLOGIA}
        r = C.simular_clima(argumentos, I, carpeta, informar=informar, arranque=args.arranque, version=VERSION,
                            opciones_nucleo=opciones_nucleo, **prueba)
    except KeyboardInterrupt:
        informar("Interrumpido. Para seguir, lanza la MISMA orden: continua desde el ultimo año guardado.")
        sys.exit(1)
    except (AssertionError, FloatingPointError, IndexError, ValueError) as e:      # valores no validos
        informar(f"ERROR: la simulacion se ha vuelto inestable ({e}). Pasale a Claude {f_reg}.")
        sys.exit(2)
    from exportar_clima import huella_mapa          # 3.16.0: con que mapa se hizo (lo comprueba la exportacion)
    r["mapa"] = {"nombre": nombre, "huella": huella_mapa(tipo, alt)} if args.tierra else \
        {"nombre": nombre_mapa, "huella": huella_mapa(tipo, alt)}
    r["sensibilidad"] = dict(SENSIBILIDAD)          # 4.0.0: factores del oceano (1 y 1 = el modelo tal cual)
    r["origen"] = origen                            # 3.21.0: climatologia continuada desde otra (o None)
    with open(os.path.join(carpeta, "resultado.pkl"), "wb") as f:
        pickle.dump(r, f)
    escribir_resumen(r, carpeta, time.time() - t0, len(orbita))
    informar(f"Terminado. Resumen en {os.path.join(carpeta, 'resumen.txt')}")


def leer_origen(carpeta):
    """3.21.0: el origen.json de una climatologia continuada desde otra (o None)."""
    import json
    f = os.path.join(carpeta, "origen.json")
    if not os.path.exists(f):
        return None
    with open(f, encoding="utf-8") as fh:
        return json.load(fh)


def preparar_continuacion(origen, carpeta, informar):
    """3.21.0: prepara <carpeta> para una climatologia nueva que parte del estado final de <origen>. Solo la primera
    vez (si la carpeta ya tiene su simulacion, se continua esa y se comprueba que el origen es el mismo)."""
    import hashlib
    import json
    import shutil
    origen = os.path.abspath(origen)
    if os.path.abspath(carpeta) == origen:
        sys.exit("--continuar-desde tiene que ser OTRA carpeta")
    ya = leer_origen(carpeta)
    if ya is not None:
        if ya["carpeta"] != origen:
            sys.exit(f"{carpeta} ya continua desde {ya['carpeta']}, no desde {origen}")
        return ya
    if any(os.path.exists(os.path.join(carpeta, n)) for n in ("estado_i16.pkl", "climatologia.pkl", "resultado.pkl")):
        sys.exit(f"{carpeta} ya tiene una simulacion propia: usa una carpeta nueva para --continuar-desde")
    for n in ("estado_i16.pkl", "climatologia.pkl", "arranque_dos_capas.pkl", "resultado.pkl"):
        if not os.path.exists(os.path.join(origen, n)):
            sys.exit(f"{origen} no tiene {n}: no es una climatologia terminada")
    with open(os.path.join(origen, "climatologia.pkl"), "rb") as fh:
        cl = pickle.load(fh)
    with open(os.path.join(origen, "estado_i16.pkl"), "rb") as fh:
        ck = pickle.load(fh)
    if not ck["registrado"] or ck["ano"] != cl["ano_sim"]:
        sys.exit(f"el punto de control de {origen} no es el del final de su climatologia (¿esta a medias?)")
    if not cl["equilibrio"]["convergido"] and not args.dias:          # con --dias (solo pruebas) no se exige
        sys.exit(f"{origen} no llego al equilibrio: no sirve de partida")
    with open(os.path.join(origen, "resultado.pkl"), "rb") as fh:
        ro = pickle.load(fh)
    if (ro.get("sensibilidad") or {"oceano_difusion": 1.0, "oceano_capa_mezcla": 1.0}) != SENSIBILIDAD:
        sys.exit(f"{origen} es de otra configuracion del oceano")
    d = ro.get("climatologia", {}).get("deriva") or {}
    if d.get("hay_deriva"):
        sys.exit(f"{origen} tiene deriva: no es un estado de equilibrio")
    os.makedirs(carpeta, exist_ok=True)
    for n in ("arranque_dos_capas.pkl", "estado_i16.pkl"):
        shutil.copyfile(os.path.join(origen, n), os.path.join(carpeta, n))
    # el comienzo del año siguiente: fase6_clima lo lee como el "estado previo" (el del final de la otra climatologia)
    shutil.copyfile(os.path.join(origen, "estado_i16.pkl"), os.path.join(carpeta, "estado_i16.pkl.previo"))
    with open(os.path.join(origen, "estado_i16.pkl"), "rb") as fh:
        huella = hashlib.sha256(fh.read()).hexdigest()
    info = {"carpeta": origen, "ano_simulado": int(ck["ano"]), "anos_equilibrio_origen": int(cl["equilibrio"]["anos"]),
            "anos_climatologia_origen": int(cl["acum"].n), "sha256_estado": huella,
            "procedencia_origen": ro.get("climatologia", {}).get("procedencia", [])}
    with open(os.path.join(carpeta, "origen.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=1, ensure_ascii=False)
    informar(f"Climatologia NUEVA desde el estado final de {origen} (año simulado {info['ano_simulado']}: "
             f"{info['anos_equilibrio_origen']} de equilibrio y {info['anos_climatologia_origen']} de climatologia)")
    return info


def escribir_resumen(r, carpeta, segundos, pasos_ano):
    c = r["climatologia"]
    eq = c["equilibrio"]
    peso = np.cos(np.radians(LATITUDES_GRADOS))
    T = c["bandas_aire2m_C"].mean(axis=0)
    Pz = c["bandas_precipitacion"].mean(axis=0)
    L = [f"RESUMEN de clima_dinamico.py (M3N {VERSION}, nucleo dinamico I16)",
         f"modo: {'Tierra' if args.tierra else 'P3N'} | arranque: {args.arranque} | pasos por año: {pasos_ano}"
         + (f" | PRUEBA con años de {args.dias} dias" if args.dias else "")
         + (f" | SENSIBILIDAD del oceano: D_oc x{args.oceano_difusion:g}, capa de mezcla x{args.oceano_capa_mezcla:g}"
            " (no vale como clima oficial)" if ES_SENSIBILIDAD else ""),
         f"tiempo de esta ejecucion: {segundos / 3600:.2f} h",
         "Validacion frente a la Tierra real (solo modo Tierra): python validar_i16.py",
         "",
         f"EQUILIBRIO: convergido {eq['convergido']} tras {eq['anos']} años"]
    og = r.get("origen")
    if og:                                          # 3.21.0
        L[-1] = (f"EQUILIBRIO: el de {og['carpeta']} ({og['anos_equilibrio_origen']} años de equilibrio y "
                 f"{og['anos_climatologia_origen']} de climatologia antes de esta; año simulado de partida "
                 f"{og['ano_simulado']}). Esta climatologia continua desde su estado final.")
    for k, (val, umbral, cumple) in eq["valores"].items():
        L.append(f"  {k:11s} {val:+.5f} (umbral {umbral}) {'ok' if cumple else 'NO'}")
    s = eq["serie"]
    L.append("  serie anual N (W/m2):    " + " ".join(f"{x:+.3f}" for x in s["N"]))
    L.append("  serie anual aire 2 m (C): " + " ".join(f"{x - 273.15:.3f}" for x in s["T2m"]))
    objetivo = c.get("anos_objetivo", c["anos_promediados"])
    tipo = ("EXPLORATORIA (no vale como clima oficial)" if c.get("exploratoria")
            else "normal climatologica de la OMM (WMO-No. 1203)" if objetivo == 30 else "")
    L += ["", f"CLIMATOLOGIA: {c['anos_promediados']} años registrados de {objetivo} | {tipo}",
          f"  error de la media (sigma/raiz(N_eff), informativo): maximo aire 2 m {c['error_T_max_K']:.3f} K | "
          f"precipitacion {c['error_P_max_mm_dia']:.3f} mm/dia ({100 * c['error_P_rel_max']:.1f} %) | bandas por encima "
          f"de la referencia (0,1 K; 5 % o 0,05 mm/dia): T {c['bandas_T_sin_cumplir']}, P {c['bandas_P_sin_cumplir']}",
          f"  aire 2 m medio global {float((T * peso).sum() / peso.sum()):.3f} C | precipitacion media global "
          f"{float((Pz * peso).sum() / peso.sum()):.3f} mm/dia", ""]
    d = c.get("deriva")
    if d:
        f = lambda t: (f"{t['pendiente']:+.4f} +- {t['error']:.4f} /año (p {t['p']:.3f})"
                       + (" SIGNIFICATIVA" if t["significativa"] else ""))
        L.append(f"DERIVA de los años registrados ({d['metodo']}; FDR {d['alfa_fdr']}):")
        L.append(f"  aire 2 m global: {f(d['global_T'])} K | precipitacion global: {f(d['global_P'])} mm/dia")
        malas = [(la, t) for la, t in zip(LATITUDES_GRADOS, d["bandas_T"]) if t["significativa"]]
        L.append("  bandas con deriva significativa: " + (", ".join(f"{la:+.1f} ({t['pendiente']:+.4f} K/año)"
                                                                 for la, t in malas) if malas else "ninguna"))
        L.append("  -> " + ("AVISO: hay deriva; la climatologia no es de un estado en equilibrio" if d["hay_deriva"]
                            else "sin deriva significativa"))
        L.append("")
    L += ["POR AÑO REGISTRADO (energia y agua):"]
    for i, d in enumerate(c["por_ano"]):
        e = d["energia"]
        L.append(f"  {i + 1:2d}: cierre energia {e.get('cierre_relativo', float('nan')):+.2e} | error dinamica "
                 f"{e.get('error_dinamica_W_m2', float('nan')):+.4f} W/m2 | resto {e.get('cierre_sin_dinamica', float('nan')):+.1e} | "
                 f"agua atm {d['cierre_agua_atmosfera_kg_m2']:+.1e} tierra {d['cierre_agua_tierra_kg_m2']:+.1e}")
    L += ["", "PROCEDENCIA (version y huella del codigo de cada tramo de esta simulacion; procedencia.json):"]
    tramos = c.get("procedencia") or []
    if not tramos:
        L.append("  no consta (simulacion empezada antes de la 3.14.0)")
    for t in tramos:
        L.append(f"  desde el año simulado {t['ano_sim']:3d} ({t['anos_registrados']:2d} registrados): M3N {t['version']}, "
                 f"codigo {t['codigo']} ({t['fecha']})")
    if og:
        L.append(f"  los {og['ano_simulado']} años anteriores son de {og['carpeta']} (su procedencia: "
                 + ", ".join(f"M3N {t['version']} {t['codigo']}" for t in og.get("procedencia_origen", [])) + ")")
    elif tramos and tramos[0]["ano_sim"] > 0:
        L.append(f"  AVISO: los {tramos[0]['ano_sim']} primeros años se hicieron antes de la 3.14.0 (codigo no registrado)")
    z = r.get("zonal")
    L.append("  medias zonales por nivel (corte vertical de H3N): " + (f"si, {z['anos']} años" if z else "no"))
    if len({t["codigo"] for t in tramos}) > 1:
        L.append("  AVISO: la simulacion mezcla versiones del codigo")
    L += ["", "BANDAS (latitud, aire 2 m C, precipitacion mm/dia):"]
    for lat, t, p in zip(LATITUDES_GRADOS, T, Pz):
        L.append(f"  {lat:+6.1f} {t:8.2f} {p:8.3f}")
    with open(os.path.join(carpeta, "resumen.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

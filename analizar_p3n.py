# analizar_p3n.py -- 4.0.0 (Fase 6.4, DISENO_FASE_6.4.md §4.7, criterio aprobado por Carlos el 11/10/2026).
# DESCRIPCION del clima de P3N con el nucleo dinamico (I16). P3N no tiene observaciones: este analisis no aprueba
# ni suspende nada por ser distinto de la Tierra. Hace cuatro cosas:
#   A. los bloqueantes (los mismos de la 6.3, DISENO_FASE_6.3.md §6.22), leidos del resumen.txt de la simulacion
#      con las mismas funciones que herramientas/criterio_cierre_6_3.py (A7, exportacion y H3N: a mano);
#   B. la descripcion de la circulacion SIN suponer su forma: todas las celulas de la funcion de corriente (no solo
#      dos buscadas en ventanas terrestres), todos los nucleos de chorro, todas las bandas de viento en superficie
#      y todos los maximos de lluvia, con sus latitudes; temperatura, hielo, transporte de calor y perfil vertical.
#      Las estaciones de P3N las manda la DISTANCIA a la estrella (§3): se usan el año orbital (doceavos desde el
#      periastro) y dos temporadas de un cuarto de año centradas en el periastro y en el apoastro;
#   B'. la orientacion teorica: la anchura de la celula de Hadley del modo Tierra multiplicada por el factor de
#      Held y Hou (1980) para el radio, la gravedad y la rotacion de P3N (solo parametros, §3), si se da la
#      validacion del modo Tierra de la misma version (--referencia-tierra);
#   C. la comparacion con el clima oficial de 2 capas (el json que lee H3N), si es del mismo mapa;
#   D. la sensibilidad al oceano (§4.2): diferencias con otras simulaciones de P3N (--comparar).
#
# Uso (en la carpeta de M3N, con el venv activado), DESPUES de "python clima_dinamico.py":
#     python analizar_p3n.py
# Opciones:
#     --carpeta RUTA            la de clima_dinamico (por defecto outputs/clima_dinamico/p3n_dos_capas)
#     --sin-ano-extra           sin la parte de un año mas (radiacion, hielo, perfil, vientos)
#     --referencia-tierra RUTA  carpeta de validar_i16.py del modo Tierra (campos.npz), para B'
#     --clima-oficial RUTA      json del clima de 2 capas (por defecto el que lee H3N)
#     --comparar RUTA [RUTA..]  carpetas de otras simulaciones de P3N (sensibilidad al oceano)
#     --dias N                  SOLO PARA PROBAR, con una carpeta de "clima_dinamico.py --dias N"
# Resultado: informe.txt y campos.npz en outputs/analisis_p3n/. Diagnostico: no cambia nada del modelo.

import argparse
import os
import sys

ap = argparse.ArgumentParser(description="Analisis del clima de P3N con el nucleo dinamico (Fase 6.4)")
ap.add_argument("--carpeta", default=None)
ap.add_argument("--sin-ano-extra", action="store_true")
ap.add_argument("--referencia-tierra", default=None)
ap.add_argument("--clima-oficial", default=None)
ap.add_argument("--comparar", nargs="+", default=[])
ap.add_argument("--dias", type=int, default=None, help="solo para probar, con una carpeta de --dias N")
ap.add_argument("--salida", default=os.path.join("outputs", "analisis_p3n"))
if __name__ == "__main__":
    args = ap.parse_args()
    if os.environ.get("M3N_MODO", "").strip().lower() == "tierra":
        sys.exit("analizar_p3n.py es para P3N: quita M3N_MODO=tierra del entorno")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
else:                                             # importado (pruebas): solo las funciones
    args = ap.parse_args([])

import json
import math
import pickle
import time

import numpy as np

import parametros as P
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS

LAT = LATITUDES_GRADOS
NIVEL_INTERIOR = (0.15, 0.95)        # interior de la columna para el extremo de cada celula (como validar_i16)
NIVEL_TROPOSFERA = (0.1, 0.5)        # nucleos de chorro de la troposfera (como validar_i16, §6.16)
NIVEL_ESTRATOSFERA = 0.0101
FRACCION_DEBIL = 0.05                # una celula o banda con |extremo| < 5 % del mayor se marca "debil" (no se oculta)
TEMPORADAS = ("periastro", "apoastro")
SALIDA = args.salida


# ------------------------------------------------------------------------------------ utilidades
def peso_filas():
    d = math.radians(180.0 / FILAS) / 2
    return np.sin(np.radians(LAT) + d) - np.sin(np.radians(LAT) - d)


def media_filas(z, mascara=None):
    w = peso_filas() * (1.0 if mascara is None else mascara)
    return float((np.asarray(z) * w).sum() / w.sum())


def calendario(orbita, pasos_dia, n_dias):
    """El año orbital de P3N: dia del periastro y del apoastro (maximo y minimo de la luz media diaria en el
    tope, que solo depende de la distancia), fase de cada dia (0 = periastro, en fraccion de año), doceavo
    (0..11) y temporada ("periastro": el cuarto de año centrado en el periastro; "apoastro": el centrado en el
    apoastro; None el resto). Con 1,7 grados de inclinacion el ciclo anual lo manda la distancia (§3)."""
    toa = np.array([o[0] for o in orbita], dtype=float)
    n_orb = -(-len(toa) // pasos_dia)                 # dias de la orbita (el ultimo, incompleto si lo es)
    diaria = np.array([toa[d * pasos_dia:(d + 1) * pasos_dia].mean() for d in range(n_orb)])
    d_peri, d_apo = int(np.argmax(diaria)), int(np.argmin(diaria))
    fase = ((np.arange(n_dias) - d_peri) / (len(toa) / pasos_dia)) % 1.0
    doceavo = np.minimum((fase * 12).astype(int), 11)
    temporada = np.where((fase < 0.125) | (fase >= 0.875), "periastro",
                         np.where((fase >= 0.375) & (fase < 0.625), "apoastro", ""))
    return {"d_peri": d_peri, "d_apo": d_apo, "fase": fase, "doceavo": doceavo, "temporada": temporada,
            "toa_peri": float(diaria[d_peri]), "toa_apo": float(diaria[d_apo])}


def media_zonal_por(diario, etiqueta, valores):
    """diario (dias, F, C) -> {valor: media zonal (F,)} de los dias con esa etiqueta."""
    z = np.asarray(diario, dtype=float).mean(axis=2)
    return {v: z[etiqueta == v].mean(axis=0) for v in valores if (etiqueta == v).any()}


def tramos_de_signo(perfil, lats):
    """Tramos consecutivos de un mismo signo de un perfil por latitudes (lats de norte a sur). Devuelve una lista
    de (signo, i0, i1, lat_norte, lat_sur), con los bordes interpolados en el cruce por cero (el borde del primer
    y el ultimo tramo es la ultima latitud)."""
    s = np.sign(perfil)
    # un cero exacto no es un tramo: toma el signo del vecino (primero el de su norte; si no hay, el de su sur)
    for i in range(len(s)):
        if s[i] == 0 and i > 0:
            s[i] = s[i - 1]
    for i in range(len(s) - 2, -1, -1):
        if s[i] == 0:
            s[i] = s[i + 1]
    out, i0 = [], 0
    for i in range(1, len(perfil) + 1):
        if i == len(perfil) or s[i] != s[i0]:
            def cruce(a, b):
                return float(lats[a] + (lats[b] - lats[a]) * perfil[a] / (perfil[a] - perfil[b]))
            norte = float(lats[0]) if i0 == 0 else cruce(i0 - 1, i0)
            sur = float(lats[-1]) if i == len(perfil) else cruce(i - 1, i)
            out.append((int(s[i0]), i0, i - 1, norte, sur))
            i0 = i
    return out


def celulas(psi, lat_c, nivel, k5):
    """TODAS las celulas de la funcion de corriente (10^10 kg/s), sin ventanas: los tramos de signo de psi en el
    nivel k5 (~sigma 0,5), cada uno con su extremo en el interior de la columna (sigma 0,15-0,95) dentro de sus
    caras, la latitud del extremo y los bordes (cruces por cero en k5, como TropD, Adam et al. 2018).
    psi > 0: el aire sube por el lado sur de la celula y baja por el norte (como la de Hadley del norte)."""
    nivel = np.asarray(nivel)
    inter = (nivel > NIVEL_INTERIOR[0]) & (nivel < NIVEL_INTERIOR[1])
    out = []
    for signo, i0, i1, ln, ls in tramos_de_signo(psi[k5], lat_c):
        sub = psi[np.ix_(inter, np.arange(i0, i1 + 1))]
        k, j = np.unravel_index(np.argmax(signo * sub), sub.shape)
        out.append({"signo": signo, "extremo": float(sub[k, j]), "lat": float(lat_c[i0 + j]),
                    "norte": ln, "sur": ls})
    mayor = max((abs(c["extremo"]) for c in out), default=0.0)
    for c in out:
        c["debil"] = abs(c["extremo"]) < FRACCION_DEBIL * mayor
    return out


def nucleos_chorro(u, lat, nivel):
    """TODOS los nucleos de chorro de la troposfera (criterio de validar_i16, §6.16: maximo local en latitud y
    altura en sigma 0,1-0,5 y con menos viento en el nivel de encima), del oeste (u > 0) o del este (u < 0, se
    buscan como maximos de -u). Lista de (m/s, lat, sigma), de norte a sur."""
    nivel = np.asarray(nivel)
    u = np.asarray(u, dtype=float)
    nk, ni = u.shape
    out = []
    for s in (1.0, -1.0):
        w = s * u
        for k in range(1, nk - 1):
            if not (NIVEL_TROPOSFERA[0] <= nivel[k] <= NIVEL_TROPOSFERA[1]):
                continue
            for i in range(ni):
                c = w[k, i]
                if c <= 0 or not (c > w[k - 1, i] and c >= w[k + 1, i]):
                    continue
                # empate con la fila del norte: el nucleo se cuenta una vez (en la del norte)
                if (i > 0 and w[k, i - 1] >= c) or (i < ni - 1 and w[k, i + 1] > c):
                    continue
                out.append((float(u[k, i]), float(lat[i]), float(nivel[k])))
    return sorted(out, key=lambda t: -t[1])


def bandas_viento(us, lat):
    """Bandas de viento zonal en superficie (capa baja): tramos de signo con su extremo. Lista de dicts."""
    out = []
    for signo, i0, i1, ln, ls in tramos_de_signo(us, lat):
        tramo = us[i0:i1 + 1]
        j = int(np.argmax(signo * tramo))
        out.append({"signo": signo, "extremo": float(tramo[j]), "lat": float(lat[i0 + j]), "norte": ln, "sur": ls})
    mayor = max((abs(b["extremo"]) for b in out), default=0.0)
    for b in out:
        b["debil"] = abs(b["extremo"]) < FRACCION_DEBIL * mayor
    return out


def extremos_locales(z, lat, maximos=True):
    """Maximos (o minimos) locales estrictos de un perfil por filas, comparando con las filas vecinas (en los
    bordes, con la unica vecina). Lista de (valor, lat), de norte a sur."""
    z = np.asarray(z, dtype=float)
    s = 1.0 if maximos else -1.0
    out = []
    for i in range(len(z)):
        vec = [z[j] for j in (i - 1, i + 1) if 0 <= j < len(z)]
        if all(s * z[i] > s * v for v in vec):
            out.append((float(z[i]), float(lat[i])))
    return out


def factor_held_hou():
    """Anchura de la celula de Hadley ∝ (g H Δθ / (Ω² a² θ0))^½ (Held y Hou 1980, J. Atmos. Sci. 37, 515, como
    proporcionalidad; Guendelman y Kaspi 2019, GRL, doi:10.1029/2018GL080752). Con el mismo contraste y la misma
    altura (§3), el factor de P3N frente a la Tierra solo depende de g, Ω y a."""
    g_t, a_t = 9.80665, 6.371e6                    # los del modo Tierra (parametros.py)
    om_t = 2 * math.pi * (1 / 86400 + 1 / (365.256 * 86400))
    return math.sqrt((P.P3N_GRAVEDAD / g_t) / ((P.OMEGA_SIDERAL / om_t) ** 2 * (P.P3N_RADIO / a_t) ** 2))


# ------------------------------------------------------------------------------------ A. bloqueantes
def seccion_a(L, carpeta):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "herramientas"))
    import criterio_cierre_6_3 as CC
    f = os.path.join(carpeta, "resumen.txt")
    L.append("A. BLOQUEANTES (los de la 6.3, DISENO_FASE_6.3.md §6.22; los mismos calculos que criterio_cierre_6_3.py)")
    if not os.path.exists(f):
        L.append(f"   falta {f}: no se pueden comprobar")
        L.append("")
        return None
    t = CC.Tabla()
    with open(f, encoding="utf-8") as fh:
        try:
            CC.seccion_a(t, fh.read())
        except (ValueError, IndexError, AttributeError) as e:      # resumen incompleto (p. ej., de --dias)
            t.criterio("A", "resumen incompleto", False, f"no se pueden leer todos los bloqueantes ({e})",
                       explicable=False)
    t.criterio("A7", "analisis completo (exportacion y H3N: a mano)", None, "este informe; exportar_clima.py y H3N a mano",
               explicable=False)
    for clave, texto, estado, detalle in t.filas:
        L.append(f"   {clave:3s} {estado:6s} {texto}: {detalle}")
    L.append("")
    return all(e != "❌" for _, _, e, _ in t.filas)


# ------------------------------------------------------------------------------------ B. parte 1
def parte_1(rc, tipo, cal, L, campos):
    from fase1_geografia import TIERRA
    import fase31_agua as AG
    from validar_i16 import transporte
    clim = rc["climatologia"]
    es_tierra = np.asarray(tipo) == TIERRA
    L.append(f"B. CLIMATOLOGIA DE P3N: {clim['anos_promediados']} años (equilibrio tras {clim['equilibrio']['anos']} "
             f"años: {clim['equilibrio']['convergido']})")
    L.append(f"   Año orbital: periastro el dia {cal['d_peri']} (luz en el tope {cal['toa_peri']:.1f} W/m2), apoastro el"
             f" dia {cal['d_apo']} ({cal['toa_apo']:.1f} W/m2). Temporadas: un cuarto de año centrado en cada uno.")
    L.append("")
    T_dia = np.asarray(rc["reg_media"], dtype=float)
    T_zon = T_dia.mean(axis=(0, 2))
    T_doc = media_zonal_por(T_dia, cal["doceavo"], range(12))
    T_doc = {m: T_doc.get(m, np.full(FILAS, np.nan)) for m in range(12)}      # (un doceavo sin dias: NaN)
    T_tmp = media_zonal_por(T_dia, cal["temporada"], TEMPORADAS)
    campos.update(T_zonal=T_zon, T_doceavos=np.array([T_doc[m] for m in range(12)]))
    norte, sur = (LAT > 0).astype(float), (LAT < 0).astype(float)
    T_cel = T_dia.mean(axis=0)
    area = peso_filas()[:, None] * np.ones((1, COLUMNAS))
    L.append("1. AIRE A 2 m (C)")
    L.append(f"   Media global {media_filas(T_zon):6.2f} | hemisferio norte {media_filas(T_zon, norte):6.2f} | "
             f"sur {media_filas(T_zon, sur):6.2f}")
    L.append(f"   Tierra {float((T_cel * area)[es_tierra].sum() / area[es_tierra].sum()):6.2f} | mar "
             f"{float((T_cel * area)[~es_tierra].sum() / area[~es_tierra].sum()):6.2f}")
    ecu = (np.abs(LAT) < 5).astype(float)
    L.append(f"   Ecuador (5 S-5 N) menos 87,5 N {media_filas(T_zon, ecu) - T_zon[0]:6.2f} | menos 87,5 S "
             f"{media_filas(T_zon, ecu) - T_zon[-1]:6.2f}")
    if all(k in T_tmp for k in TEMPORADAS):
        L.append(f"   Temporada del periastro menos la del apoastro, global: "
                 f"{media_filas(T_tmp['periastro']) - media_filas(T_tmp['apoastro']):+6.2f}")
    L.append("   Por filas: lat | anual | ciclo anual (doceavo mas calido - mas frio) | periastro - apoastro")
    for i in range(FILAS):
        serie = np.array([T_doc[m][i] for m in range(12)])
        dif = (T_tmp["periastro"][i] - T_tmp["apoastro"][i]) if all(k in T_tmp for k in TEMPORADAS) else float("nan")
        amp = np.nanmax(serie) - np.nanmin(serie) if np.isfinite(serie).sum() > 1 else float("nan")
        L.append(f"     {LAT[i]:+5.1f} | {T_zon[i]:7.2f} | {amp:6.1f} | {dif:+6.1f}")
    L.append("")

    ag = rc["agua"]
    Pz = np.asarray(ag["precipitacion"]).mean(axis=1)
    Ez = np.asarray(ag["evaporacion"]).mean(axis=1)
    campos.update(P_zonal=Pz, E_zonal=Ez)
    L.append("2. PRECIPITACION (mm/dia)")
    L.append(f"   Media global {media_filas(Pz):6.3f} | evaporacion {media_filas(Ez):6.3f} | agua precipitable "
             f"{media_filas(np.asarray(ag['agua_precipitable']).mean(axis=1)):5.1f} kg/m2")
    L.append(f"   En tierra {float((np.asarray(ag['precipitacion']) * area)[es_tierra].sum() / area[es_tierra].sum()):6.3f}"
             f" | en el mar {float((np.asarray(ag['precipitacion']) * area)[~es_tierra].sum() / area[~es_tierra].sum()):6.3f}")
    diario = ag.get("diario_precipitacion")
    P_tmp = media_zonal_por(diario, cal["temporada"], TEMPORADAS) if diario is not None else {}
    L.append("   Maximos de lluvia de la media zonal (todos, de norte a sur; mm/dia a latitud):")
    for nom, z in (("año", Pz),) + tuple((t, P_tmp[t]) for t in TEMPORADAS if t in P_tmp):
        L.append(f"     {nom:9s} " + "; ".join(f"{v:5.2f} a {la:+5.1f}" for v, la in extremos_locales(z, LAT)))
    L.append("   Minimos (zonas secas de la media zonal, año): "
             + "; ".join(f"{v:5.2f} a {la:+5.1f}" for v, la in extremos_locales(Pz, LAT, maximos=False)))
    L.append("   Por filas: lat | precipitacion | evaporacion | P - E")
    for i in range(FILAS):
        L.append(f"     {LAT[i]:+5.1f} | {Pz[i]:5.2f} | {Ez[i]:5.2f} | {Pz[i] - Ez[i]:+5.2f}")
    L.append("")

    fl = rc["flujos"]
    bordes, tot, resto_toa = transporte(-np.asarray(fl["toa_neto"], dtype=float))
    _, oc, resto_oc = transporte(fl["conv_oceano"])
    atm = tot - oc
    _, lat_t, _ = transporte(AG.L_V * (np.asarray(ag["precipitacion"]) - np.asarray(ag["evaporacion"])) / 86400.0)
    campos.update(bordes=bordes, transporte_total=tot, transporte_oceano=oc, transporte_atmosfera=atm,
                  transporte_latente=lat_t)
    n, s = bordes > 0, bordes < 0
    L.append("3. TRANSPORTE DE CALOR HACIA LOS POLOS (PW)   norte            sur")
    for nom, x in (("Total", tot), ("Atmosfera", atm), ("Oceano", oc), ("Latente", lat_t)):
        L.append(f"   {nom + ', maximo':32s} {x[n].max():5.2f} a {bordes[n][np.argmax(x[n])]:+3.0f}     "
                 f"{-x[s].min():5.2f} a {bordes[s][np.argmin(x[s])]:+3.0f}")
    i0 = int(np.argmin(np.abs(bordes)))
    L.append(f"   A traves del ecuador (>0 = hacia el norte): atmosfera {atm[i0]:+5.2f} PW, oceano {oc[i0]:+5.2f} PW")
    L.append(f"   (desequilibrio global quitado antes de integrar: tope {-resto_toa:+.3f}, oceano {resto_oc:+.3f} W/m2)")
    L.append("")
    return T_zon, Pz


# ------------------------------------------------------------------------------------ B. parte 2
def parte_2(r, vientos, L, campos, ref_tierra=None):
    import validar_v30 as V30
    from validar_i16 import funcion_corriente, perfil_banda
    L.append("B (cont.). UN AÑO MAS desde el final de la climatologia")
    e = r["energia"]
    L.append(f"   cierre de energia {e.get('cierre_relativo', float('nan')):+.1e}; error de la dinamica "
             f"{e.get('error_dinamica_W_m2', float('nan')):+.3f} W/m2")
    L.append("")
    f = V30.informe_flujos(r["flujos"])
    absorbida = f["sw_suelo"] + f["sw_atm"]
    L.append("4. BALANCE DE RADIACION (W/m2; sin nubes hasta la Fase 5.2)")
    L.append(f"   Luz que llega {f['entrada']:6.1f} | reflejada {f['reflejada']:6.1f} (albedo {f['reflejada'] / f['entrada']:.3f})"
             f" | absorbida {absorbida:6.1f} (suelo {f['sw_suelo']:5.1f}, aire {f['sw_atm']:5.1f})")
    L.append(f"   Infrarrojo al espacio {f['olr']:6.1f} | neto en el tope {absorbida - f['olr']:+6.2f} | infrarrojo al suelo"
             f" {f['dlr']:6.1f}")
    L.append("")
    h = V30.hielo(r)
    L.append("5. HIELO MARINO, extension (millones de km2 de P3N; celda entera si hay hielo)")
    L.append(f"   Norte: maximo {h['N_max']:5.2f}, minimo {h['N_min']:5.2f} | Sur: maximo {h['S_max']:5.2f}, minimo {h['S_min']:5.2f}")
    L.append("")
    L.append("6. PERFIL VERTICAL (media anual por bandas)")
    L.append("   Banda      | tropopausa OMM | interpolada | T tropop. interp. | mas frio (500-50 hPa) | gradiente 0-6 km")
    for b in V30.perfiles_por_bandas(r["flujos"], P.P3N_GRAVEDAD):
        if b is None:
            continue
        banda, Tz, zz, pz = b
        d = perfil_banda(Tz, zz, pz)
        g06 = float((Tz[0] - np.interp(6000, zz, Tz)) / ((6000 - zz[0]) / 1000))
        L.append(f"   {banda:10s} | {d['tropopausa_km']:5.1f} km       | {d['tropopausa_interp_km']:5.1f} km    | "
                 f"{d['T_tropopausa_interp']:6.1f} C          | {d['T_frio']:6.1f} C a {d['frio_km']:4.1f} km   | {g06:4.2f} K/km")
    L.append("")
    if vientos is None or "anual" not in vientos:
        return
    sig = vientos["sigma"]
    sm = 0.5 * (sig[1:] + sig[:-1])
    lat_c = vientos["lat_caras"]
    k5 = int(np.argmin(np.abs(sig - 0.5)))
    L.append("7. VIENTOS DEL NUCLEO (sin suponer su forma: todo lo que hay, de norte a sur)")
    L.append("   Celulas (funcion de corriente de masa, 10^10 kg/s; >0: sube por su lado sur y baja por el norte;")
    L.append("   extremo en sigma 0,15-0,95; bordes = cruces por cero a sigma 0,5; 'debil' = menos del 5 % de la mayor):")
    resumen_cel = {}
    for clave in ("anual",) + TEMPORADAS:
        if clave not in vientos:
            continue
        psi = funcion_corriente(vientos[clave]["vps"], sig, vientos["cos_caras"], P.P3N_RADIO, P.P3N_GRAVEDAD) / 1e10
        campos["psi_" + clave] = psi
        campos["u_" + clave] = vientos[clave]["u"]
        cs = celulas(psi, lat_c, sig, k5)
        resumen_cel[clave] = cs
        L.append(f"     {clave}: {len(cs)} celulas ({sum(not c['debil'] for c in cs)} no debiles)")
        for c in cs:
            L.append(f"       {c['extremo']:+7.2f} a {c['lat']:+5.1f} | de {c['norte']:+5.1f} a {c['sur']:+5.1f}"
                     + (" (debil)" if c["debil"] else ""))
    u = vientos["anual"]["u"]
    L.append("   Nucleos de chorro de la troposfera (maximos locales en latitud y altura, sigma 0,1-0,5; >0 del oeste):")
    for clave in ("anual",) + TEMPORADAS:
        if clave not in vientos:
            continue
        nc = nucleos_chorro(vientos[clave]["u"], LAT, sm)
        L.append(f"     {clave:9s} " + ("; ".join(f"{v:+5.1f} m/s a {la:+5.1f}, {1000 * s:4.0f} hPa" for v, la, s in nc)
                                       if nc else "ninguno"))
    L.append("   Estratosfera alta (sigma <= 0,01), viento maximo del oeste por hemisferio (media anual):")
    sel = sm <= NIVEL_ESTRATOSFERA
    for hemi, m in (("norte", LAT > 0), ("sur", LAT < 0)):
        sub = np.where(sel[:, None] & m[None, :], u, -np.inf)
        k, i = np.unravel_index(np.argmax(sub), sub.shape)
        L.append(f"     {hemi:5s} {u[k, i]:+5.1f} m/s a {LAT[i]:+5.1f}, {1000 * sm[k]:5.1f} hPa")
    L.append("   Bandas de viento en la capa baja (media anual; <0 del este, >0 del oeste; extremo y bordes):")
    for b in bandas_viento(u[-1], LAT):
        L.append(f"     {'oeste' if b['signo'] > 0 else 'este ':5s} {b['extremo']:+5.1f} m/s a {b['lat']:+5.1f} | de "
                 f"{b['norte']:+5.1f} a {b['sur']:+5.1f}" + (" (debil)" if b["debil"] else ""))
    L.append("   Viento zonal medio anual por filas (m/s): lat | sigma ~0,25 | ~0,5 | capa baja")
    k25, k50 = int(np.argmin(np.abs(sm - 0.25))), int(np.argmin(np.abs(sm - 0.5)))
    for i in range(0, FILAS, 2):
        L.append(f"     {LAT[i]:+5.1f} | {u[k25, i]:+6.1f} | {u[k50, i]:+6.1f} | {u[-1, i]:+6.1f}")
    L.append("")

    # B'. orientacion teorica
    fh = factor_held_hou()
    L.append("B'. ORIENTACION TEORICA (no es un requisito: §4.7)")
    L.append("   Factor de Held y Hou (1980) para la anchura de la celula de Hadley, P3N/Tierra, con el mismo contraste")
    L.append(f"   de temperatura y la misma altura (solo g, Ω y a): {fh:.3f}")
    if ref_tierra is not None and "anual" in resumen_cel:
        psi_t = ref_tierra["psi_anual"]
        cs_t = celulas(psi_t, lat_c, sig, k5) if psi_t.shape == campos["psi_anual"].shape else None
        if cs_t:
            L.append("   Celulas que tocan el ecuador (las tropicales), borde hacia el polo: modo Tierra x factor | P3N")
            for hemi in ("norte", "sur"):
                borde = (lambda c: c["norte"]) if hemi == "norte" else (lambda c: c["sur"])
                toca = (lambda c: c["sur"] < 5 and c["norte"] > 5) if hemi == "norte" else \
                    (lambda c: c["norte"] > -5 and c["sur"] < -5)
                bt = [borde(c) for c in cs_t if not c["debil"] and toca(c)]
                bp = [borde(c) for c in resumen_cel["anual"] if not c["debil"] and toca(c)]
                if bt and bp:
                    L.append(f"     {hemi:5s} {bt[0]:+5.1f} x {fh:.3f} = {bt[0] * fh:+5.1f} | P3N {bp[0]:+5.1f}")
                else:
                    L.append(f"     {hemi:5s} sin celula que cruce de 5 S a 5 N hacia este lado en "
                             f"{'el modo Tierra' if not bt else 'P3N'}")
        else:
            L.append("   (la referencia del modo Tierra no tiene la misma rejilla vertical: no se compara)")
    elif ref_tierra is None:
        L.append("   (sin --referencia-tierra no se calcula la anchura esperada)")
    L.append("   Con mas rotacion, las celulas se estrechan, aparecen mas chorros y crece el contraste ecuador-polo")
    L.append("   (Kaspi y Showman 2015, Astrophys. J. 804, 60); P3N gira 1,21 veces mas rapido pero es mas pequeño.")
    L.append("")


# ------------------------------------------------------------------------------------ C y D
def seccion_c(L, T_zon, ruta_oficial, huella):
    L.append("C. COMPARACION CON EL CLIMA OFICIAL DE 2 CAPAS (el json que lee H3N)")
    if not ruta_oficial or not os.path.exists(ruta_oficial):
        L.append(f"   no esta {ruta_oficial}: no se compara")
        L.append("")
        return
    with open(ruta_oficial, encoding="utf-8") as fh:
        d = json.load(fh)
    if d.get("mapa", {}).get("huella") != huella:
        L.append(f"   el clima oficial es de otro mapa ({d.get('mapa', {}).get('nombre')}): no se compara")
        L.append("")
        return
    from exportar_clima import leer_serie
    T_of = leer_serie(d["diario"]["aire2m_media"])
    T_of = T_of.reshape(T_of.shape[0], FILAS, COLUMNAS).mean(axis=(0, 2))
    L.append(f"   M3N {d.get('m3n', {}).get('version')} ({d.get('simulacion', {}).get('nivel_modelo', '2 capas')})")
    L.append(f"   Aire a 2 m global: oficial {media_filas(T_of):6.2f} C | con nucleo {media_filas(T_zon):6.2f} C | "
             f"diferencia {media_filas(T_zon) - media_filas(T_of):+5.2f}")
    L.append("   Por filas: lat | oficial | con nucleo | diferencia")
    for i in range(FILAS):
        L.append(f"     {LAT[i]:+5.1f} | {T_of[i]:7.2f} | {T_zon[i]:7.2f} | {T_zon[i] - T_of[i]:+6.2f}")
    L.append("")


def seccion_d(L, rc, carpetas, tipo):
    from fase1_geografia import TIERRA
    if not carpetas:
        return
    L.append("D. SENSIBILIDAD (§4.2): otra simulacion menos esta (aire a 2 m en C, precipitacion en mm/dia)")
    T0 = np.asarray(rc["reg_media"], dtype=float).mean(axis=0)
    P0 = np.asarray(rc["agua"]["precipitacion"], dtype=float)
    es_tierra = np.asarray(tipo) == TIERRA
    for c in carpetas:
        with open(os.path.join(c, "resultado.pkl"), "rb") as fh:
            r = pickle.load(fh)
        if r.get("mapa") != rc.get("mapa"):
            L.append(f"   {c}: es de otro mapa, no se compara")
            continue
        dT = np.asarray(r["reg_media"], dtype=float).mean(axis=0) - T0
        dP = np.asarray(r["agua"]["precipitacion"], dtype=float) - P0
        s = r.get("sensibilidad") or {}
        L.append(f"   {c}  ({r['climatologia']['anos_promediados']} años; D_oc x{s.get('oceano_difusion', 1.0):g}, capa de"
                 f" mezcla x{s.get('oceano_capa_mezcla', 1.0):g})")
        L.append(f"     global: dT {media_filas(dT.mean(axis=1)):+5.2f} | dP {media_filas(dP.mean(axis=1)):+6.3f}"
                 f" | en tierra: dT entre {dT[es_tierra].min():+5.2f} y {dT[es_tierra].max():+5.2f}, |dT| medio "
                 f"{np.abs(dT[es_tierra]).mean():4.2f}")
        L.append("     por filas (cada 10 grados): " + "; ".join(f"{LAT[i]:+5.1f} {dT[i].mean():+5.2f}"
                                                            for i in range(0, FILAS, 2)))
    L.append("")


# ------------------------------------------------------------------------------------ main
def main():
    from puente_c3n import cargar_mapa_activo_de_c3n
    from exportar_clima import huella_mapa, RUTA_SALIDA
    import validar_i16 as V
    carpeta = args.carpeta or os.path.join("outputs", "clima_dinamico",
                                           "p3n_dos_capas" + (f"_prueba{args.dias}d" if args.dias else ""))
    f_res = os.path.join(carpeta, "resultado.pkl")
    if P.MODO_TIERRA:
        sys.exit("analizar_p3n.py es para P3N")
    if not os.path.exists(f_res):
        sys.exit(f"Falta {f_res}: hay que terminar antes 'python clima_dinamico.py'.")
    os.makedirs(SALIDA, exist_ok=True)
    with open(f_res, "rb") as fh:
        rc = pickle.load(fh)
    tipo, alt, nombre = cargar_mapa_activo_de_c3n()
    huella = huella_mapa(tipo, alt)
    if rc.get("mapa") is None or rc["mapa"]["huella"] != huella:
        sys.exit(f"La simulacion de {carpeta} no es del mapa activo ({nombre}): activa en C3N el mapa con el que se hizo.")
    from cache_simulacion import precalcular_orbita_cacheada, huella_codigo, MODULOS_I16
    from temperatura import PASO_TIEMPO
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    pasos_dia = round(P.ROTACION_PERIODO / PASO_TIEMPO)
    if args.dias:
        orbita = orbita[:args.dias * pasos_dia]
    n_dias = np.asarray(rc["reg_media"]).shape[0]
    cal = calendario(orbita, pasos_dia, n_dias)
    L = [f"ANALISIS DEL CLIMA DE P3N CON EL NUCLEO DINAMICO (analizar_p3n.py) -- {time.strftime('%Y-%m-%d %H:%M')}",
         f"carpeta: {carpeta} | mapa {nombre} ({huella}) | codigo de este analisis: {huella_codigo(MODULOS_I16)}",
         "P3N no tiene observaciones: este informe describe y explica; no aprueba ni suspende por diferencias con la"
         " Tierra (DISENO_FASE_6.4.md §4.7).", ""]
    campos = {}
    seccion_a(L, carpeta)
    T_zon, _ = parte_1(rc, tipo, cal, L, campos)
    sensibilidad = any(v != 1.0 for v in (rc.get("sensibilidad") or {}).values())
    if sensibilidad and not args.sin_ano_extra:
        L.append("(Esta carpeta es una prueba de sensibilidad del oceano: el año extra se hace solo con la simulacion"
                 " principal. Aqui, solo la climatologia.)")
        L.append("")
    if not args.sin_ano_extra and not sensibilidad:
        argumentos = (orbita, tipo, alt, P.EMISIVIDAD) + argumentos_resto()
        I = interruptores_i16()
        V.SALIDA = SALIDA
        temporada = cal["temporada"]
        print("Un año mas desde una copia del ultimo punto de control (la simulacion larga no se toca)...", flush=True)
        r, vientos = V.ano_extra(carpeta, argumentos, I,
                                 temporada_del_dia=lambda d: (temporada[d % len(temporada)] or None),
                                 temporadas=TEMPORADAS)
        ref = None
        if args.referencia_tierra:
            f = os.path.join(args.referencia_tierra, "campos.npz")
            ref = dict(np.load(f)) if os.path.exists(f) else None
            if ref is None:
                L.append(f"(no esta {f}: sin referencia del modo Tierra)")
        parte_2(r, vientos, L, campos, ref)
    seccion_c(L, T_zon, args.clima_oficial or RUTA_SALIDA, huella)
    seccion_d(L, rc, args.comparar, tipo)
    texto = "\n".join(L)
    print(texto)
    with open(os.path.join(SALIDA, "informe.txt"), "w", encoding="utf-8") as fh:
        fh.write(texto + "\n")
    np.savez_compressed(os.path.join(SALIDA, "campos.npz"), **campos)
    print(f"\nGuardado en {os.path.join(SALIDA, 'informe.txt')} (y campos.npz): pasale informe.txt a Claude.")


def argumentos_resto():
    """Lo que clima_dinamico.py pasa despues de la emisividad (mismos argumentos que la simulacion larga)."""
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
    return (ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55)


def interruptores_i16():
    import fase2b_atmosfera as F
    from validar_i16 import INTERRUPTORES_I16
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in INTERRUPTORES_I16:
        I[k] = True
    return I


if __name__ == "__main__":
    main()

# validar_v30.py -- v3.0: validacion de la atmosfera de N capas (DISENO_FASE_2.3.md, seccion 10).
#
# Uso (terminal B, venv de M3N), DESPUES de calibrar D con calibrar_v30.py:
#     M3N_MODO=tierra python validar_v30.py     # la Tierra frente a la realidad
#     python validar_v30.py                     # P3N: v2.4.3 frente a v3.0
# (opcional: --d-atm X --d-oc Y para probar otros coeficientes)
#
# --d-atm es el valor de la TIERRA (como sale de la calibracion); en P3N se
# multiplica aqui por parametros.FACTOR_ROTACION_D. Sin --d-atm se usa
# fase2b_atmosfera.D_ATMOSFERA_V30 (y sin --d-oc, D_OCEANO_V30): los calibrados.
#
# Las simulaciones van con el hielo marino encendido (validacion completa);
# en modo Tierra converge en unos 40 años.
#
# Que mide (todo con el año final, medias anuales):
#   - balance global: luz absorbida en el suelo y en el aire, luz
#     reflejada, infrarrojo hacia el suelo (DLR) y hacia el espacio (OLR);
#   - temperatura del aire a 2 m por bandas de latitud;
#   - perfil vertical medio por bandas: altura de la TROPOPAUSA (definicion
#     de la OMM, 1957: el nivel mas bajo en el que el gradiente baja a
#     2 K/km o menos y no vuelve a superarlos de media en los 2 km de
#     encima; buscada entre 550 y 75 hPa como Reichler et al. 2003),
#     temperatura en la tropopausa y a 25 km, gradiente medio
#     entre 0 y 6 km;
#   - hielo marino: extension maxima y minima de cada hemisferio;
#   - transporte de calor hacia los polos (atmosfera y oceano).
# Los valores observados de referencia (la Tierra) se imprimen al lado
# solo donde hay una fuente verificada; ver DISENO_FASE_2.3.md.

import os
import sys
import json
import time
import argparse

import numpy as np
import parametros as P
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, GRADOS_POR_FILA, GRADOS_POR_COL

CARPETA = os.path.join("outputs", "validacion_v30")
BANDAS = [(90, 60), (60, 30), (30, 10), (10, -10), (-10, -30), (-30, -60), (-60, -90)]

# Referencias de la Tierra (solo las verificadas):
# Wild et al. (2019), cielo despejado, W/m2 (DISENO_FASE_2.2.md, tabla de la luz): llega al suelo 247 y el
# suelo absorbe 214 (v3.1: antes se comparaba por error lo absorbido con lo que llega), atmosfera 73,
# reflejada 53 (incluida la del suelo); DLR 314; OLR 267 (estos dos, con la atmosfera REAL, calentada
# tambien por las nubes: el modelo sin nubes debe dar un OLR mayor, igual a lo absorbido).
# NSIDC, extension del hielo marino, media 1981-2010 (DISENO_FASE_3.md, 6.1): Artico max 15,6 / min ~6;
# Antartico max ~18-19 / min ~3 millones de km2.
REF_WILD = {"sw_suelo": 214, "sw_atm": 73, "reflejada": 53, "dlr": 314, "olr": 267}
REF_HIELO = {"N_max": 15.6, "N_min": 6, "S_max": 18.5, "S_min": 3}


def area_celdas():
    lat = np.radians(LATITUDES_GRADOS)
    dphi, dlam = np.radians(GRADOS_POR_FILA), np.radians(GRADOS_POR_COL)
    a = P.P3N_RADIO ** 2 * dlam * (np.sin(lat + dphi / 2) - np.sin(lat - dphi / 2))
    return a[:, None] * np.ones((1, COLUMNAS))


def media_area(x, mascara=None):
    w = area_celdas() if mascara is None else area_celdas() * mascara
    return float((x * w).sum() / w.sum())


def alturas(T, p, ps, gravedad):
    """Altura (m) del centro de cada capa sobre el suelo, integrando la hidrostatica con la T de
    cada capa (T y p de arriba abajo, en la primera dimension; ps, presion en el suelo)."""
    from fase30_multicapa import R_AIRE
    n = T.shape[0]
    z = np.zeros_like(T)
    z[-1] = R_AIRE * T[-1] / gravedad * np.log(ps / p[-1])
    for k in range(n - 2, -1, -1):
        Tm = 0.5 * (T[k] + T[k + 1])
        z[k] = z[k + 1] + R_AIRE * Tm / gravedad * np.log(p[k + 1] / p[k])
    return z


# Rango de busqueda de Reichler, Dameris y Sausen (2003, GRL 30, 2042) ✅: entre 550 y 75 hPa (~5-18 km),
# para no confundir con la tropopausa las inversiones junto al suelo (noche polar).
P_MAX_TROPOPAUSA = 5.5e4   # Pa
P_MIN_TROPOPAUSA = 7.5e3   # Pa


def tropopausa_omm(Tz, zz, pz):
    """Tropopausa de la OMM sobre un perfil (de abajo arriba, z en m, p en Pa). Devuelve (z, T) o (nan, nan)."""
    gam = -np.diff(Tz) / np.diff(zz) * 1000          # K/km entre niveles consecutivos
    for i in range(len(gam)):
        if pz[i] > P_MAX_TROPOPAUSA:
            continue
        if pz[i] < P_MIN_TROPOPAUSA:
            break
        if gam[i] <= 2.0:
            z0, T0 = zz[i], Tz[i]
            dentro = (zz > z0) & (zz <= z0 + 2000)
            if not dentro.any():
                return z0, T0
            if np.all((T0 - Tz[dentro]) / ((zz[dentro] - z0) / 1000) <= 2.0):
                return z0, T0
    return np.nan, np.nan


def perfil_por_bandas(fl, gravedad):
    """Media zonal (por area) del perfil medio anual en cada banda: tropopausa, T a 25 km, gradiente 0-6 km."""
    T = fl["T_atm"]; p = fl["p_capas"]; ps_c = fl["p_superficie"]
    filas = []
    for (n, s) in BANDAS:
        m = (LATITUDES_GRADOS < n) & (LATITUDES_GRADOS > s)
        # media de las celdas de la banda con el suelo por debajo de ~500 m (presion en superficie > 940 hPa):
        # el perfil sobre las montañas empieza en otra base
        Tb, pb, sb, wb = [], [], [], []
        for i in np.where(m)[0]:
            for j in range(COLUMNAS):
                if ps_c[i, j] > 9.4e4:
                    Tb.append(T[:, i, j]); pb.append(p[:, i, j]); sb.append(ps_c[i, j]); wb.append(np.cos(np.radians(LATITUDES_GRADOS[i])))
        if not Tb:
            filas.append(None); continue
        wb = np.array(wb)
        Tm = (np.array(Tb) * wb[:, None]).sum(0) / wb.sum()
        pm = (np.array(pb) * wb[:, None]).sum(0) / wb.sum()
        z = alturas(Tm, pm, float((np.array(sb) * wb).sum() / wb.sum()), gravedad)
        Tz, zz, pz = Tm[::-1], z[::-1], pm[::-1]       # de abajo arriba
        zt, Tt = tropopausa_omm(Tz, zz, pz)
        T25 = float(np.interp(25000, zz, Tz)) if zz[-1] > 25000 else np.nan
        # v3.1: desde el centro de la capa mas baja (np.interp no extrapola por debajo de el)
        g06 = float((Tz[0] - np.interp(6000, zz, Tz)) / ((6000 - zz[0]) / 1000))
        filas.append({"banda": f"{n}..{s}", "tropopausa_km": zt / 1000, "T_tropopausa": Tt - 273.15,
                      "T_25km": T25 - 273.15, "gradiente_0_6km": g06})
    return filas


def hielo(r):
    h = r["hielo_espesor"]                          # (dias, F, C)
    a = area_celdas() / 1e12                          # millones de km2
    norte = LATITUDES_GRADOS[:, None] > 0
    N = ((h > 0) * a * norte).sum(axis=(1, 2)); S = ((h > 0) * a * ~norte).sum(axis=(1, 2))
    return {"N_max": float(N.max()), "N_min": float(N.min()), "S_max": float(S.max()), "S_min": float(S.min())}


def simular(tipo, alt, multi, d_atm, d_oc, max_anos=60, v31=False, vapor_radiativo=False):
    import fase2b_atmosfera as F
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
    from cache_simulacion import precalcular_orbita_cacheada
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    I = dict(F.INTERRUPTORES_FASE2B)
    I["atmosfera_multicapa"] = multi
    if v31:
        for k in ("ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "albedo_espectral"):
            I[k] = True
        I["vapor_radiativo"] = vapor_radiativo
    return F.simular_fase2b(orbita, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
                            P.PROFUNDIDAD_OPTICA, 0.55, interruptores=I,
                            d_atmosfera=d_atm, d_oceano=d_oc, max_anos=max_anos)


def bandas_aire(r, tipo):
    from fase1_geografia import TIERRA
    media = r["reg_media"].mean(axis=0)
    out = []
    for (n, s) in BANDAS:
        m = ((LATITUDES_GRADOS < n) & (LATITUDES_GRADOS > s))[:, None] * np.ones((1, COLUMNAS), bool)
        t = m & (tipo == TIERRA); a = m & (tipo != TIERRA)
        out.append((f"{n}..{s}", media_area(media, m), media_area(media, t) if t.any() else np.nan,
                    media_area(media, a) if a.any() else np.nan))
    return out


def informe_flujos(fl):
    sw_s, sw_a = media_area(fl["sw_suelo"]), media_area(fl["sw_atm"])
    entrada = P.S3N_LUMINOSIDAD / (16 * np.pi * P.SEMIEJE_MAYOR ** 2 * np.sqrt(1 - P.ORBITA_EXCENTRICIDAD ** 2))
    return {"entrada": entrada, "sw_suelo": sw_s, "sw_atm": sw_a, "reflejada": entrada - sw_s - sw_a,
            "dlr": media_area(fl["dlr"]), "olr": media_area(fl["olr_celda"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d-atm", type=float, default=None, help="D de la atmosfera, valor de la Tierra")
    ap.add_argument("--d-oc", type=float, default=None)
    ap.add_argument("--max-anos", type=int, default=60)
    ap.add_argument("--v31", action="store_true", help="v3.1: con el ciclo del agua (I11-I14)")
    ap.add_argument("--vapor-radiativo", action="store_true", help="v3.1: con el prototipo de vapor radiativo (I15)")
    a = ap.parse_args()
    import fase2b_atmosfera as F
    import fase30_multicapa as M30
    d_tierra = F.D_ATMOSFERA_V30 if a.d_atm is None else a.d_atm
    d_atm = d_tierra * P.FACTOR_ROTACION_D
    d_oc = F.D_OCEANO_V30 if a.d_oc is None else a.d_oc
    os.makedirs(CARPETA, exist_ok=True)
    if P.MODO_TIERRA:
        from modo_tierra import mapa_tierra
        tipo, alt = mapa_tierra(); nombre = "Tierra"
    else:
        from puente_c3n import cargar_mapa_activo_de_c3n
        tipo, alt, nombre = cargar_mapa_activo_de_c3n()
    print(f"Validacion v3.0 -- {nombre} -- D atm {d_tierra} (Tierra) x {P.FACTOR_ROTACION_D:.3f} = {d_atm:.3f}, D oc {d_oc}", flush=True)

    t = time.time()
    r3 = simular(tipo, alt, True, d_atm, d_oc, a.max_anos, v31=a.v31 or a.vapor_radiativo,
                 vapor_radiativo=a.vapor_radiativo)
    print(f"v3.0: {r3['anos']} años, {time.time() - t:.0f} s", flush=True)
    r2 = None
    if not P.MODO_TIERRA:
        t = time.time()
        r2 = simular(tipo, alt, False, None, None, a.max_anos)      # v2.4.3 tal cual (sus D de siempre)
        print(f"v2.4.3: {r2['anos']} años, {time.time() - t:.0f} s", flush=True)

    fl = r3["flujos"]
    flu = informe_flujos(fl)
    perf = perfil_por_bandas(fl, P.P3N_GRAVEDAD)
    hi3 = hielo(r3)
    b3 = bandas_aire(r3, tipo)
    bordes, atm = M30.transporte_meridional(fl["conv_atmosfera"], P.P3N_RADIO)
    _, oc = M30.transporte_meridional(fl["conv_oceano"], P.P3N_RADIO)

    L = []
    L.append(f"VALIDACION {'v3.1' if a.v31 else 'v3.0'} -- {nombre} -- {time.strftime('%Y-%m-%d')}")
    L.append(f"D atm = {d_tierra} x {P.FACTOR_ROTACION_D:.3f} = {d_atm:.3f}; D oc = {d_oc}; N = {M30.N_CAPAS_ATM} capas; "
             f"{r3['anos']} años; desequilibrio {r3['energia']['diferencia_relativa']:.1e}")
    L.append("")
    L.append("BALANCE GLOBAL (W/m2)        modelo" + ("   Tierra sin nubes (Wild 2019)" if P.MODO_TIERRA else ""))
    for k, nom in (("entrada", "Luz que llega (media)"), ("sw_suelo", "Luz absorbida, suelo"), ("sw_atm", "Luz absorbida, aire"),
                   ("reflejada", "Luz reflejada"), ("dlr", "Infrarrojo al suelo (DLR)"), ("olr", "Infrarrojo al espacio (OLR)")):
        ref = f"{REF_WILD[k]:8.0f}" if P.MODO_TIERRA and k in REF_WILD else ""
        L.append(f"  {nom:28s} {flu[k]:7.1f} {ref}")
    if P.MODO_TIERRA:
        L.append("  (DLR y OLR observados son con la atmosfera real, calentada tambien por las nubes; sin nubes, el OLR")
        L.append("   del modelo debe igualar lo absorbido. Ver DISENO_FASE_2.3.md.)")
    L.append("")
    L.append(f"AIRE A 2 m, media anual (C): global {media_area(r3['reg_media'].mean(0)):.2f}"
             + (f" (v2.4.3: {media_area(r2['reg_media'].mean(0)):.2f})" if r2 else ""))
    L.append("  Banda      |   v3.0: todo  tierra   agua" + ("  |  v2.4.3: todo  tierra   agua" if r2 else ""))
    b2 = bandas_aire(r2, tipo) if r2 else None
    for i, (bn, x, t_, w) in enumerate(b3):
        extra = f"  |        {b2[i][1]:6.1f} {b2[i][2]:6.1f} {b2[i][3]:6.1f}" if r2 else ""
        L.append(f"  {bn:10s} |        {x:6.1f} {t_:6.1f} {w:6.1f}{extra}")
    L.append("")
    L.append("PERFIL VERTICAL (medias anuales, celdas con el suelo por debajo de ~500 m)")
    L.append("  Banda      | tropopausa (OMM) | T tropopausa | T a 25 km | gradiente 0-6 km")
    for f in perf:
        if f is None:
            continue
        L.append(f"  {f['banda']:10s} |   {f['tropopausa_km']:5.1f} km       |  {f['T_tropopausa']:6.1f} C   | "
                 f"{f['T_25km']:6.1f} C  |  {f['gradiente_0_6km']:4.2f} K/km")
    L.append("")
    L.append("HIELO MARINO (millones de km2; celda entera si hay hielo)" + ("        observado (NSIDC 1981-2010)" if P.MODO_TIERRA else ""))
    hi2 = hielo(r2) if r2 else None
    for k, nom in (("N_max", "Norte, maximo"), ("N_min", "Norte, minimo"), ("S_max", "Sur, maximo"), ("S_min", "Sur, minimo")):
        extra = f"   {REF_HIELO[k]:5.1f}" if P.MODO_TIERRA else (f"   (v2.4.3: {hi2[k]:5.1f})" if hi2 else "")
        L.append(f"  {nom:15s} {hi3[k]:6.1f}{extra}")
    L.append("")
    tot = atm + oc
    n, s = bordes > 0, bordes < 0
    L.append(f"TRANSPORTE HACIA LOS POLOS: atmosfera max {atm[n].max():.2f} PW ({bordes[n][np.argmax(atm[n])]:+.0f}) / "
             f"{-atm[s].min():.2f} PW ({bordes[s][np.argmin(atm[s])]:+.0f}); oceano max {oc[n].max():.2f} / {-oc[s].min():.2f} PW")
    if P.MODO_TIERRA:
        L.append("  observado (Trenberth y Caron 2001): atmosfera 5,0 PW a 43 N, parecido hacia 40 S")
    ag = r3.get("agua")
    if ag:
        # v3.1: ciclo del agua. Referencias verificadas: GPCP v2.3 (Adler et al. 2018) 2,69 mm/dia (tierra
        # 2,24, oceano 2,90); agua precipitable ~24,9 kg/m2 (Trenberth y Smith 2005: 1,27e16 kg de vapor);
        # transporte latente ~ la mitad del atmosferico en latitudes medias (Hwang y Frierson 2010).
        from fase1_geografia import TIERRA
        bt, la = M30.transporte_meridional(ag["transporte_latente_conv"], P.P3N_RADIO)
        i40n, i40s = int(np.argmin(np.abs(bordes - 40))), int(np.argmin(np.abs(bordes + 40)))
        L.append("")
        L.append("CICLO DEL AGUA (v3.1)" + ("            observado" if P.MODO_TIERRA else ""))
        for nom, x, ref in (("Precipitacion global (mm/dia)", media_area(ag["precipitacion"]), "2,69 (GPCP)"),
                            ("  sobre tierra", media_area(ag["precipitacion"], tipo == TIERRA), "2,24"),
                            ("  sobre el oceano", media_area(ag["precipitacion"], tipo != TIERRA), "2,90"),
                            ("Evaporacion global (mm/dia)", media_area(ag["evaporacion"]), ""),
                            ("Agua precipitable (kg/m2)", media_area(ag["agua_precipitable"]), "~24,9"),
                            ("Nieve: fraccion de la precipitacion", media_area(ag["nieve"]) / max(media_area(ag["precipitacion"]), 1e-12), ""),
                            ("Cubo medio en tierra (kg/m2)", media_area(ag["cubo_medio"], tipo == TIERRA), ""),
                            ("Parte latente del transporte a 40 N / 40 S", f"{la[i40n] / atm[i40n]:.2f} / {la[i40s] / atm[i40s]:.2f}", "~0,5")):
            xs = f"{x:8.3f}" if isinstance(x, float) else f"{x:>8s}"
            L.append(f"  {nom:42s} {xs}   {ref if P.MODO_TIERRA else ''}")
        L.append(f"  Cierres: agua de la atmosfera {ag['cierre_agua_atmosfera_kg_m2']:.1e}, de la tierra "
                 f"{ag['cierre_agua_tierra_kg_m2']:.1e} kg/m2; energia {r3['energia'].get('cierre_relativo') or 0:.1e}")
        L.append("  Por bandas: precipitacion / evaporacion (mm/dia)")
        for (bn, b0, b1) in [(f"{n}..{s_}", n, s_) for (n, s_) in BANDAS]:
            m = ((LATITUDES_GRADOS < b0) & (LATITUDES_GRADOS > b1))[:, None] * np.ones((1, COLUMNAS), bool)
            L.append(f"    {bn:10s}  {media_area(ag['precipitacion'], m):6.2f} / {media_area(ag['evaporacion'], m):6.2f}")
    texto = "\n".join(L)
    print(); print(texto)
    base = os.path.join(CARPETA, f"validacion_{'v31vr' if a.vapor_radiativo else 'v31' if a.v31 else 'v30'}_{nombre}_{time.strftime('%Y-%m-%d_%H%M')}")
    open(base + ".txt", "w", encoding="utf-8").write(texto + "\n")
    np.savez_compressed(base + ".npz", **{k: v for k, v in fl.items()},
                        reg_media_v30=r3["reg_media"], **({"reg_media_v243": r2["reg_media"]} if r2 else {}))
    print(f"\nGuardado en {base}.txt (y .npz con los campos)")


if __name__ == "__main__":
    main()

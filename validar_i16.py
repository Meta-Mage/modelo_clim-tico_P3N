# validar_i16.py -- v3.1-pre16: VALIDACION del modo Tierra de M3N con el nucleo dinamico (I16) frente a la
# Tierra real (DISENO_FASE6_3.md §6.14). Diagnostico: no cambia nada del modelo ni de la simulacion larga.
#
# Uso (en la carpeta de M3N, con el venv activado), DESPUES de "python clima_dinamico.py --tierra":
#     python validar_i16.py
# Opciones:
#     --carpeta RUTA     la de clima_dinamico (por defecto outputs/clima_dinamico/tierra_dos_capas)
#     --sin-ano-extra    solo la parte 1 (sin el año extra)
#     --dias N           SOLO PARA PROBAR, con una carpeta de "clima_dinamico.py --tierra --dias N"
#
# Que hace:
#   PARTE 1 (al momento) -- con la climatologia de la simulacion larga (resultado.pkl, media de N años):
#     - aire a 2 m: media global y de cada hemisferio, por bandas de latitud y ciclo anual, frente a la
#       climatologia CRU 1961-1990 (Jones et al. 1999), que esta en la misma rejilla de 5 grados;
#     - precipitacion: global, por bandas, la banda de lluvias ecuatorial (ITCZ) en enero, julio y el año,
#       frente a GPCP (Huffman et al. 2009; Adler et al. 2018); evaporacion y agua precipitable;
#     - transporte de calor hacia los polos, total (desde el balance en el tope), del oceano y de la atmosfera
#       (el resto), y la parte latente, frente a Trenberth y Caron (2001).
#   PARTE 2 (un año mas, ~15 min en el PC de Carlos) -- desde una COPIA del ultimo punto de control (la
#     simulacion larga no se toca), con todos los diagnosticos de ese año:
#     - balance de radiacion en el tope y en la superficie frente a CERES EBAF (Loeb et al. 2018), con cielo
#       despejado y con nubes, y Wild et al. (2019): M3N NO TIENE NUBES (Fase 5b), asi que lo comparable es el
#       cielo despejado; la diferencia con el cielo real es el efecto de las nubes;
#     - hielo marino (extension maxima y minima) frente a NSIDC;
#     - perfil vertical por bandas (tropopausa, gradiente 0-6 km);
#     - VIENTOS del nucleo: corrientes en chorro, vientos en superficie (alisios y del oeste) y celulas de
#       Hadley (funcion de corriente de masa), anual, diciembre-febrero y junio-agosto. Todavia sin
#       referencia observada verificada: se dan los valores del modelo.
# Los valores observados estan en referencias_tierra.py, cada uno con su fuente. Resultado: informe.txt y
# campos.npz en outputs/validacion_i16/.

import argparse
import os
import sys

ap = argparse.ArgumentParser(description="Validacion del modo Tierra de M3N con I16 frente a la Tierra real")
ap.add_argument("--carpeta", default=None)
ap.add_argument("--sin-ano-extra", action="store_true")
ap.add_argument("--dias", type=int, default=None, help="solo para probar, con una carpeta de --dias N")
ap.add_argument("--salida", default=os.path.join("outputs", "validacion_i16"))
if __name__ == "__main__":
    args = ap.parse_args()
    os.environ["M3N_MODO"] = "tierra"             # antes de importar parametros: siempre modo Tierra
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
else:                                             # importado (pruebas): solo las funciones, sin tocar el entorno
    args = ap.parse_args([])

import hashlib
import math
import pickle
import shutil
import time

import numpy as np

import parametros as P
import referencias_tierra as R
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS

LAT = LATITUDES_GRADOS
SALIDA = args.salida
INTERRUPTORES_I16 = ("atmosfera_multicapa", "ciclo_agua", "conveccion_humeda", "suelo_termico_agua",
                     "albedo_espectral", "nucleo_dinamico")          # los mismos que clima_dinamico.py
NOMBRES_MES = ("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")


# ---------------------------------------------------------------------------------------- utilidades
def peso_filas():
    """Area relativa de cada fila (proporcional a sin(lat + 2,5) - sin(lat - 2,5))."""
    d = math.radians(180.0 / FILAS) / 2
    return np.sin(np.radians(LAT) + d) - np.sin(np.radians(LAT) - d)


def media_filas(z, mascara=None):
    """Media por area de un perfil por filas (36,), opcionalmente solo en las filas de la mascara."""
    w = peso_filas() * (1.0 if mascara is None else mascara)
    return float((np.asarray(z) * w).sum() / w.sum())


def mes_de_cada_dia(n_dias):
    """Mes (0 = enero) de cada dia de la orbita del modo Tierra (dia 0 = solsticio de diciembre)."""
    fin = np.cumsum(R.DIAS_MES)
    ano = fin[-1]
    cal = (np.arange(n_dias) + 0.5 - R.DIA_1_ENERO) % ano          # dia del año civil (0 = 1 de enero)
    return np.searchsorted(fin, cal, side="right")


def mes_del_dia(dia):
    """Mes (0 = enero) del dia 'dia' de la orbita del modo Tierra."""
    return int(np.searchsorted(np.cumsum(R.DIAS_MES), (dia + 0.5 - R.DIA_1_ENERO) % sum(R.DIAS_MES), side="right"))


def medias_mensuales_zonales(diario):
    """diario: (dias, FILAS, COLUMNAS) -> (12, FILAS), media zonal de cada mes."""
    mes = mes_de_cada_dia(diario.shape[0])
    zonal = diario.mean(axis=2)
    return np.array([zonal[mes == m].mean(axis=0) if (mes == m).any() else np.full(FILAS, np.nan)
                     for m in range(12)])


def transporte(convergencia):
    """Transporte hacia el norte (PW) en los bordes entre filas a partir de una convergencia (W/m2) con
    media global cero (se le quita antes la media; se devuelve tambien lo quitado, W/m2)."""
    import fase30_multicapa as M30
    c = np.asarray(convergencia, dtype=float)
    w = np.cos(np.radians(LAT))[:, None] * np.ones((1, COLUMNAS))       # el area del operador de M30
    resto = float((c * w).sum() / w.sum())
    bordes, t = M30.transporte_meridional(c - resto, P.P3N_RADIO)
    return bordes, t, resto


def funcion_corriente(vps_zonal, sigma_semi, cos_caras, radio, gravedad):
    """Funcion de corriente de masa (kg/s) en las caras entre filas y en las interfaces sigma, integrada
    desde el tope: psi_{k+1/2} = 2 pi a cos(phi) / g * sum_{j<=k} [v p_s]_j dsigma_j.
    vps_zonal: (N capas, FILAS-1) media zonal de v*p_s en las caras. Devuelve (N+1, FILAS-1)."""
    dsig = np.diff(sigma_semi)
    acum = np.concatenate([np.zeros((1, vps_zonal.shape[1])), np.cumsum(vps_zonal * dsig[:, None], axis=0)])
    return 2 * math.pi * radio * cos_caras[None, :] / gravedad * acum


def cruce_hacia_el_polo(perfil, lats, i0, paso):
    """Latitud (interpolada) donde el perfil cambia de signo yendo desde el indice i0 en la direccion paso."""
    s = np.sign(perfil[i0])
    i = i0
    while 0 <= i + paso < len(perfil):
        if np.sign(perfil[i + paso]) != s:
            a, b = perfil[i], perfil[i + paso]
            return float(lats[i] + (lats[i + paso] - lats[i]) * a / (a - b))
        i += paso
    return float("nan")


# ------------------------------------------------------------------------------------ preparar
def argumentos_simulacion():
    """Los mismos argumentos y la misma configuracion que clima_dinamico.py --tierra."""
    import fase2b_atmosfera as F
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
    from cache_simulacion import precalcular_orbita_cacheada
    from modo_tierra import mapa_tierra
    from temperatura import PASO_TIEMPO
    orbita = precalcular_orbita_cacheada(P.S3N_LUMINOSIDAD, P.INCLINACION_AXIAL_RAD, P.SEMIEJE_MAYOR)
    if args.dias:
        orbita = orbita[:args.dias * round(P.ROTACION_PERIODO / PASO_TIEMPO)]
    tipo, alt = mapa_tierra()
    I = dict(F.INTERRUPTORES_FASE2B)
    for k in INTERRUPTORES_I16:
        I[k] = True
    return (orbita, tipo, alt, P.EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, P.PROFUNDIDAD_OPTICA, 0.55), I


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------ parte 2: año extra
def ano_extra(carpeta, argumentos, I):
    """Un año mas desde una COPIA del ultimo punto de control, con todos los diagnosticos y los vientos del
    nucleo. Se guarda en SALIDA/ano_extra.pkl y se reutiliza si el punto de control no ha cambiado."""
    import fase2b_atmosfera as F
    import fase6_nucleo as N6
    f_sim = os.path.join(carpeta, "estado_i16.pkl")
    f_clima = os.path.join(carpeta, "climatologia.pkl")
    for f in (f_sim, f_clima, os.path.join(carpeta, "arranque_dos_capas.pkl")):
        if not os.path.exists(f):
            sys.exit(f"Falta {f}: hay que terminar antes 'python clima_dinamico.py --tierra'.")
    huella = sha256(f_sim)
    f_guardado = os.path.join(SALIDA, "ano_extra.pkl")
    if os.path.exists(f_guardado):
        with open(f_guardado, "rb") as f:
            d = pickle.load(f)
        if d["huella"] == huella:
            print("Año extra: ya estaba hecho con este punto de control (se reutiliza).", flush=True)
            return d["r"], d["vientos"]
    with open(f_clima, "rb") as f:
        ano_sim = pickle.load(f)["ano_sim"]
    with open(f_sim, "rb") as f:
        ck = pickle.load(f)
    if not ck["registrado"] or ck["ano"] != ano_sim:
        sys.exit("El punto de control no es el del final de la climatologia (¿esta a medias?). "
                 "Termina antes 'python clima_dinamico.py --tierra'.")
    with open(os.path.join(carpeta, "arranque_dos_capas.pkl"), "rb") as f:
        estado_ini = pickle.load(f)["estado"]
    copia = os.path.join(SALIDA, "estado_i16_copia.pkl")
    shutil.copyfile(f_sim, copia)
    for extra in (copia + ".previo",):
        if os.path.exists(extra):
            os.remove(extra)

    # vientos: se leen (sin cambiar nada) en cada subpaso del nucleo
    from temperatura import PASO_TIEMPO
    pasos_dia = round(P.ROTACION_PERIODO / PASO_TIEMPO)
    V = {"n": 0}
    original = N6.NucleoSeco.avanzar

    def avanzar(self, ant, act, forzamiento=None, **kw):
        res = original(self, ant, act, forzamiento, **kw)
        ps, T, u, v = res[1]
        dia = (V["n"] // 2) // pasos_dia
        mes = mes_del_dia(dia)
        temporada = "dic-feb" if mes in (11, 0, 1) else "jun-ago" if mes in (5, 6, 7) else None
        ps_cara = 0.5 * (ps[1:] + ps[:-1])
        campos = {"u": u.mean(axis=2), "vps": (v * ps_cara[None]).mean(axis=2), "ps": ps.mean(axis=1)}
        for clave in ("anual",) + ((temporada,) if temporada else ()):
            d = V.setdefault(clave, {"n": 0, **{k: 0.0 for k in campos}})
            d["n"] += 1
            for k, x in campos.items():
                d[k] = d[k] + x
        if "sigma" not in V:
            V["sigma"] = self.sh.copy()
            V["cos_caras"] = np.cos(self.Rj.phi_f[1:-1])
            V["lat_caras"] = np.degrees(self.Rj.phi_f[1:-1])
        V["n"] += 1
        return res

    t0 = time.time()

    def al_paso(ano, p, total, registrado):
        if p % 200 == 0 or p + 1 == total:
            hecho = (p + 1) / total
            vel = (p + 1) / max(time.time() - t0, 1e-9)
            b = "#" * int(25 * hecho) + "-" * (25 - int(25 * hecho))
            sys.stdout.write(f"\rAño extra [{b}] {100 * hecho:5.1f} %  {vel:5.1f} pasos/s  quedan "
                             f"{(total - p - 1) / max(vel, 1e-9) / 60:5.1f} min ")
            sys.stdout.flush()

    N6.NucleoSeco.avanzar = avanzar
    F.AL_PASO = al_paso
    try:
        r = F.simular_fase2b(*argumentos, max_anos=ck["ano"], interruptores=I, estado_inicial=estado_ini,
                             archivo_estado=copia, guardar_al_terminar=True)
    finally:
        N6.NucleoSeco.avanzar = original
        F.AL_PASO = None
    print(flush=True)
    for extra in (copia, copia + ".previo"):
        if os.path.exists(extra):
            os.remove(extra)
    r = {k: r[k] for k in ("flujos", "energia", "hielo_espesor")}     # lo que usa la parte 2
    vientos = {"sigma": V["sigma"], "cos_caras": V["cos_caras"], "lat_caras": V["lat_caras"]}
    for clave in ("anual", "dic-feb", "jun-ago"):
        if clave in V:
            d = V[clave]
            vientos[clave] = {k: d[k] / d["n"] for k in ("u", "vps", "ps")}
    with open(f_guardado + ".tmp", "wb") as f:
        pickle.dump({"huella": huella, "r": r, "vientos": vientos}, f)
    os.replace(f_guardado + ".tmp", f_guardado)
    return r, vientos


# ------------------------------------------------------------------------------------ informe
def parte_1(rc, tipo, L, campos):
    from fase1_geografia import TIERRA
    import fase31_agua as AG
    clim = rc["climatologia"]
    L.append(f"PARTE 1 -- CLIMATOLOGIA DE LA SIMULACION LARGA: {clim['anos_promediados']} años "
             f"(equilibrio tras {clim['equilibrio']['anos']} años: {clim['equilibrio']['convergido']})")
    L.append("")

    # --- aire a 2 m
    T_dia = np.asarray(rc["reg_media"], dtype=float)              # (dias, F, C), C
    T_zon = T_dia.mean(axis=(0, 2))
    T_mes = medias_mensuales_zonales(T_dia)
    cru, cru_mes = np.array(R.CRU_T_ANUAL), np.array(R.CRU_T_MES)
    norte, sur = (LAT > 0).astype(float), (LAT < 0).astype(float)
    campos.update(T_zonal=T_zon, T_mes=T_mes)
    L.append("1. AIRE A 2 m (C)                     M3N    observado (CRU 1961-1990, Jones et al. 1999)   M3N - obs.")
    for nom, m in (("Media global", None), ("Hemisferio norte", norte), ("Hemisferio sur", sur)):
        a, b = media_filas(T_zon, m), media_filas(cru, m)
        L.append(f"   {nom:34s} {a:6.2f}   {b:6.2f}{'':38s}{a - b:+6.2f}")
    L.append(f"   (otras referencias globales: " + "; ".join(f"{k} {v:.2f}" for k, v in R.REF_T_GLOBAL.items()) + ")")
    a, b = media_filas(T_zon, norte) - media_filas(T_zon, sur), media_filas(cru, norte) - media_filas(cru, sur)
    L.append(f"   Norte menos sur                    {a:+6.2f}   {b:+6.2f}")
    ecu = (np.abs(LAT) < 5).astype(float)
    for nom, fila in (("Ecuador (5 S-5 N) menos 87,5 N", 0), ("Ecuador (5 S-5 N) menos 87,5 S", FILAS - 1)):
        a, b = media_filas(T_zon, ecu) - T_zon[fila], media_filas(cru, ecu) - cru[fila]
        L.append(f"   {nom:34s} {a:6.2f}   {b:6.2f}")
    L.append("   Por bandas: lat | M3N anual | obs. | M3N - obs. | ciclo anual (mes mas calido - mas frio): M3N  obs.")
    for i in range(FILAS):
        amp, amp_o = np.nanmax(T_mes[:, i]) - np.nanmin(T_mes[:, i]), cru_mes[:, i].max() - cru_mes[:, i].min()
        L.append(f"     {LAT[i]:+5.1f} | {T_zon[i]:7.2f} | {cru[i]:7.2f} | {T_zon[i] - cru[i]:+6.2f} |"
                 f"{'':38s}{amp:6.1f} {amp_o:6.1f}")
    L.append("   (CRU: aire a 2 m en tierra y temperatura del agua en el mar; M3N: aire a 2 m en todas partes)")
    L.append("")

    # --- precipitacion
    ag = rc["agua"]
    Pz = np.asarray(ag["precipitacion"]).mean(axis=1)
    Ez = np.asarray(ag["evaporacion"]).mean(axis=1)
    gp, gp_mes = np.array(R.GPCP_P_ANUAL), np.array(R.GPCP_P_MES)
    campos.update(P_zonal=Pz, E_zonal=Ez)
    L.append("2. PRECIPITACION (mm/dia)              M3N    observado")
    L.append(f"   Media global                     {media_filas(Pz):6.3f}   {media_filas(gp):6.3f} (GPCP v2.2 1979-2010); "
             f"{R.REF_P_GLOBAL:.2f} (GPCP v2.3, Adler et al. 2018)")
    L.append(f"   Evaporacion global               {media_filas(Ez):6.3f}   (= precipitacion en equilibrio)")
    L.append(f"   Agua precipitable (kg/m2)        {media_filas(np.asarray(ag['agua_precipitable']).mean(axis=1)):6.1f}"
             f"   {R.REF_AGUA_PRECIPITABLE:.1f} (Trenberth y Smith 2005)")
    tropico = np.abs(LAT) < 25
    idx = np.where(tropico)[0]
    L.append("   Banda de lluvias ecuatorial (ITCZ): latitud del maximo de la media zonal (entre 25 S y 25 N)")
    diario = ag.get("diario_precipitacion")
    P_mes = medias_mensuales_zonales(np.asarray(diario, dtype=float)) if diario is not None else None
    for nom, mod, obs in (("año", Pz, gp),
                          ("enero", None if P_mes is None else P_mes[0], gp_mes[0]),
                          ("julio", None if P_mes is None else P_mes[6], gp_mes[6])):
        im = idx[np.nanargmax(mod[idx])] if mod is not None and np.isfinite(mod[idx]).any() else None
        io = idx[np.argmax(obs[idx])]
        L.append(f"     {nom:6s}  M3N {LAT[im]:+5.1f} ({mod[im]:5.2f} mm/dia)   obs. {LAT[io]:+5.1f} ({obs[io]:5.2f})"
                 if im is not None else f"     {nom:6s}  M3N  --   obs. {LAT[io]:+5.1f} ({obs[io]:5.2f})")
    for nom, m in (("Subtropicos norte (15-40 N): minimo", (LAT > 15) & (LAT < 40)),
                   ("Subtropicos sur (15-40 S): minimo", (LAT < -15) & (LAT > -40))):
        i_m, i_o = np.where(m)[0][np.argmin(Pz[m])], np.where(m)[0][np.argmin(gp[m])]
        L.append(f"   {nom:36s} M3N {Pz[i_m]:5.2f} a {LAT[i_m]:+5.1f}   obs. {gp[i_o]:5.2f} a {LAT[i_o]:+5.1f}")
    L.append("   Por bandas: lat | M3N | obs. (GPCP v2.2) | M3N/obs. | evaporacion M3N")
    for i in range(FILAS):
        L.append(f"     {LAT[i]:+5.1f} | {Pz[i]:5.2f} | {gp[i]:5.2f} | {Pz[i] / gp[i]:5.2f} | {Ez[i]:5.2f}")
    L.append("")

    # --- transporte de calor
    fl = rc["flujos"]
    # en equilibrio, lo que cada celda recibe por transporte (convergencia) es lo que pierde en el tope
    bordes, tot, resto_toa = transporte(-np.asarray(fl["toa_neto"], dtype=float))
    _, oc, resto_oc = transporte(fl["conv_oceano"])
    atm = tot - oc
    conv_lat = AG.L_V * (np.asarray(ag["precipitacion"]) - np.asarray(ag["evaporacion"])) / 86400.0   # W/m2
    _, lat_t, _ = transporte(conv_lat)
    campos.update(bordes=bordes, transporte_total=tot, transporte_oceano=oc, transporte_atmosfera=atm,
                  transporte_latente=lat_t)
    n, s = bordes > 0, bordes < 0
    i35n, i35s = int(np.argmin(np.abs(bordes - 35))), int(np.argmin(np.abs(bordes + 35)))
    i40n, i40s = int(np.argmin(np.abs(bordes - 40))), int(np.argmin(np.abs(bordes + 40)))
    L.append("3. TRANSPORTE DE CALOR HACIA LOS POLOS (PW)   norte            sur")
    L.append(f"   Total, maximo                    {tot[n].max():5.2f} a {bordes[n][np.argmax(tot[n])]:+3.0f}     "
             f"{-tot[s].min():5.2f} a {bordes[s][np.argmin(tot[s])]:+3.0f}")
    L.append(f"   Atmosfera, maximo                {atm[n].max():5.2f} a {bordes[n][np.argmax(atm[n])]:+3.0f}     "
             f"{-atm[s].min():5.2f} a {bordes[s][np.argmin(atm[s])]:+3.0f}     obs. {R.REF_TRANSPORTE['atm_max_PW']:.1f} a "
             f"{R.REF_TRANSPORTE['atm_max_lat']:.0f} N, parecido hacia 40 S (Trenberth y Caron 2001)")
    L.append(f"   Oceano, maximo                   {oc[n].max():5.2f} a {bordes[n][np.argmax(oc[n])]:+3.0f}     "
             f"{-oc[s].min():5.2f} a {bordes[s][np.argmin(oc[s])]:+3.0f}")
    L.append(f"   Parte del oceano a 35 grados     {oc[i35n] / tot[i35n]:5.2f}            {oc[i35s] / tot[i35s]:5.2f}"
             f"            obs. {R.REF_TRANSPORTE['oceano_35N']:.2f} / {R.REF_TRANSPORTE['oceano_35S']:.2f}")
    L.append(f"   Parte latente de la atmosfera a 40   {lat_t[i40n] / atm[i40n]:5.2f}        {lat_t[i40s] / atm[i40s]:5.2f}")
    L.append(f"   (total desde el balance en el tope; oceano desde la convergencia de su difusion; atmosfera = total -"
             f" oceano; latente = L_v (P - E). Desequilibrio global quitado antes de integrar: tope {-resto_toa:+.3f},"
             f" oceano {resto_oc:+.3f} W/m2)")
    L.append("   Por bordes: lat | total | atmosfera | oceano | latente")
    for i in range(len(bordes)):
        if abs(bordes[i]) % 10 == 5 or abs(bordes[i]) % 10 == 0:
            L.append(f"     {bordes[i]:+5.0f} | {tot[i]:+6.2f} | {atm[i]:+6.2f} | {oc[i]:+6.2f} | {lat_t[i]:+6.2f}")
    L.append("")
    L.append("   Para comparar, la v3.1 con la MISMA fisica pero SIN nucleo (transporte por difusion; DISENO_V3.1.md"
             " §10 bis, 05/10/2026, superficie): global ~19,0 C, polos +8,5 C (norte) y -6,2 C (sur), sin hielo marino.")
    L.append("")


def parte_2(r, vientos, L, campos):
    import validar_v30 as V30
    L.append("PARTE 2 -- UN AÑO MAS desde el final de la climatologia (todos los diagnosticos de ese año)")
    e = r["energia"]
    L.append(f"   cierre de energia {e.get('cierre_relativo', float('nan')):+.1e}; error de la dinamica "
             f"{e.get('error_dinamica_W_m2', float('nan')):+.3f} W/m2 (lo devuelve el corrector: "
             f"{e.get('correccion_energia_W_m2', float('nan')):+.3f})")
    L.append("")
    fl = r["flujos"]
    f = V30.informe_flujos(fl)
    c = R.REF_CERES
    absorbida = f["sw_suelo"] + f["sw_atm"]
    neto = absorbida - f["olr"]
    L.append("4. BALANCE DE RADIACION (W/m2)      M3N   | CERES cielo despejado | CERES con nubes (real) | Wild 2019 despejado")
    filas = (("Luz que llega", f["entrada"], c["entrada"], c["entrada"], None),
             ("Luz reflejada", f["reflejada"], c["reflejada_despejado"], c["reflejada_todo_cielo"], None),
             ("Albedo planetario", f["reflejada"] / f["entrada"], c["reflejada_despejado"] / c["entrada"],
              c["reflejada_todo_cielo"] / c["entrada"], None),
             ("Luz absorbida (total)", absorbida, c["entrada"] - c["reflejada_despejado"],
              c["entrada"] - c["reflejada_todo_cielo"], None),
             ("  en el suelo", f["sw_suelo"], None, None, R.REF_WILD["sw_suelo"]),
             ("  en el aire", f["sw_atm"], None, None, R.REF_WILD["sw_atm"]),
             ("Infrarrojo al espacio (OLR)", f["olr"], c["olr_despejado"], c["olr_todo_cielo"], None),
             ("Neto en el tope", neto, c["entrada"] - c["reflejada_despejado"] - c["olr_despejado"],
              c["neto_todo_cielo"], None),
             ("Infrarrojo al suelo (DLR)", f["dlr"], None, None, R.REF_WILD["dlr"]))
    fmt = lambda x: "      --   " if x is None else (f"{x:10.3f} " if abs(x) < 1 else f"{x:10.1f} ")
    for nom, m, d, t, w in filas:
        L.append(f"   {nom:30s} {fmt(m)}|  {fmt(d)}          |  {fmt(t)}           |  {fmt(w)}")
    L.append(f"   Efecto observado de las nubes (CERES, con nubes - despejado): luz {c['reflejada_despejado'] - c['reflejada_todo_cielo']:+.1f},"
             f" infrarrojo {c['olr_despejado'] - c['olr_todo_cielo']:+.1f}, neto "
             f"{c['reflejada_despejado'] - c['reflejada_todo_cielo'] + c['olr_despejado'] - c['olr_todo_cielo']:+.1f} W/m2."
             " M3N no tiene nubes (Fase 5b).")
    L.append("   (CERES EBAF Ed4.0, jul 2005-jun 2015, Loeb et al. 2018, tabla 5. El cielo despejado de CERES no esta en"
             " equilibrio: su neto es lo que las nubes le quitan a la Tierra real)")
    L.append("")
    h = V30.hielo(r)
    L.append("5. HIELO MARINO, extension (millones de km2; celda entera si hay hielo)   M3N   observado (NSIDC 1981-2010)")
    for k, nom in (("N_max", "Norte, maximo"), ("N_min", "Norte, minimo"), ("S_max", "Sur, maximo"), ("S_min", "Sur, minimo")):
        L.append(f"   {nom:16s} {h[k]:6.1f}   {R.REF_HIELO[k]:5.1f}")
    L.append("")
    L.append("6. PERFIL VERTICAL (media anual; celdas con el suelo por debajo de ~500 m)")
    L.append("   Banda      | tropopausa (OMM) | T tropopausa | T a 25 km | gradiente 0-6 km")
    for p in V30.perfil_por_bandas(fl, P.P3N_GRAVEDAD):
        if p is not None:
            L.append(f"   {p['banda']:10s} |   {p['tropopausa_km']:5.1f} km       |  {p['T_tropopausa']:6.1f} C   | "
                     f"{p['T_25km']:6.1f} C  |  {p['gradiente_0_6km']:4.2f} K/km")
    L.append("   (sin referencia observada verificada en el repositorio todavia)")
    L.append("")
    if vientos is None or "anual" not in vientos:
        return
    sig = vientos["sigma"]
    sm = 0.5 * (sig[1:] + sig[:-1])
    lat_c = vientos["lat_caras"]
    L.append("7. VIENTOS DEL NUCLEO (sin referencia observada verificada todavia: valores del modelo)")
    u = vientos["anual"]["u"]                                         # (N, FILAS) media zonal
    for nom, m in (("Chorro del norte", (LAT > 15) & (LAT < 70)), ("Chorro del sur", (LAT < -15) & (LAT > -70))):
        arriba = sm < 0.5
        sub = np.where(arriba[:, None] & m[None, :], u, -np.inf)
        k, i = np.unravel_index(np.argmax(sub), sub.shape)
        L.append(f"   {nom:18s} (media anual, maximo del viento del oeste): {u[k, i]:5.1f} m/s a {LAT[i]:+5.1f},"
                 f" sigma {sm[k]:.2f} (~{1000 * sm[k]:.0f} hPa)")
    us = u[-1]
    for nom, m, f_ in (("Alisios norte (0-30 N): minimo", (LAT > 0) & (LAT < 30), np.argmin),
                       ("Alisios sur (0-30 S): minimo", (LAT < 0) & (LAT > -30), np.argmin),
                       ("Oeste norte (30-70 N): maximo", (LAT > 30) & (LAT < 70), np.argmax),
                       ("Oeste sur (30-70 S): maximo", (LAT < -30) & (LAT > -70), np.argmax)):
        i = np.where(m)[0][f_(us[m])]
        L.append(f"   Capa baja, {nom:32s} {us[i]:+5.1f} m/s a {LAT[i]:+5.1f}")
    L.append("   Celulas de Hadley (funcion de corriente de masa, 10^10 kg/s; >0 = la del norte, <0 = la del sur):")
    k5 = int(np.argmin(np.abs(sig - 0.5)))                            # interfaz mas cercana a sigma 0,5
    for clave in ("anual", "dic-feb", "jun-ago"):
        if clave not in vientos:
            continue
        psi = funcion_corriente(vientos[clave]["vps"], sig, vientos["cos_caras"], P.P3N_RADIO, P.P3N_GRAVEDAD) / 1e10
        campos["psi_" + clave] = psi
        campos["u_" + clave] = vientos[clave]["u"]
        trop = np.abs(lat_c) < 40
        inter = (sig > 0.15) & (sig < 0.95)
        sub = np.where(inter[:, None] & trop[None, :], psi, np.nan)
        kn, jn = np.unravel_index(np.nanargmax(sub), sub.shape)
        ks, js = np.unravel_index(np.nanargmin(sub), sub.shape)
        borde_n = cruce_hacia_el_polo(psi[k5], lat_c, jn, -1)        # lat_c va de norte a sur
        borde_s = cruce_hacia_el_polo(psi[k5], lat_c, js, +1)
        L.append(f"     {clave:8s} norte {psi[kn, jn]:+6.2f} a {lat_c[jn]:+5.1f}, borde (sigma 0,5) {borde_n:+5.1f} | "
                 f"sur {psi[ks, js]:+6.2f} a {lat_c[js]:+5.1f}, borde {borde_s:+5.1f}")
    L.append("   Viento zonal medio anual por bandas (m/s): lat | sigma ~0,25 | ~0,5 | capa baja")
    k25, k50 = int(np.argmin(np.abs(sm - 0.25))), int(np.argmin(np.abs(sm - 0.5)))
    for i in range(0, FILAS, 2):
        L.append(f"     {LAT[i]:+5.1f} | {u[k25, i]:+6.1f} | {u[k50, i]:+6.1f} | {u[-1, i]:+6.1f}")
    L.append("")


def main():
    from modo_tierra import mapa_tierra
    carpeta = args.carpeta or os.path.join("outputs", "clima_dinamico",
                                           "tierra_dos_capas" + (f"_prueba{args.dias}d" if args.dias else ""))
    f_res = os.path.join(carpeta, "resultado.pkl")
    if not P.MODO_TIERRA:
        sys.exit("validar_i16.py es solo para el modo Tierra")
    if not os.path.exists(f_res):
        sys.exit(f"Falta {f_res}: hay que terminar antes 'python clima_dinamico.py --tierra'.")
    os.makedirs(SALIDA, exist_ok=True)
    with open(f_res, "rb") as f:
        rc = pickle.load(f)
    tipo, _ = mapa_tierra()
    L = [f"VALIDACION DEL MODO TIERRA CON EL NUCLEO DINAMICO (validar_i16.py) -- {time.strftime('%Y-%m-%d %H:%M')}",
         f"carpeta: {carpeta}", ""]
    campos = {}
    parte_1(rc, tipo, L, campos)
    if not args.sin_ano_extra:
        argumentos, I = argumentos_simulacion()
        print("Parte 2: un año mas desde una copia del ultimo punto de control (la simulacion larga no se toca)...",
              flush=True)
        r, vientos = ano_extra(carpeta, argumentos, I)
        parte_2(r, vientos, L, campos)
    texto = "\n".join(L)
    print(texto)
    with open(os.path.join(SALIDA, "informe.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    np.savez_compressed(os.path.join(SALIDA, "campos.npz"), **campos)
    print(f"\nGuardado en {os.path.join(SALIDA, 'informe.txt')} (y campos.npz): pasale informe.txt a Claude.")


if __name__ == "__main__":
    main()

# validar_i16.py -- 3.12.0; 3.13.0: referencias de vientos, Hadley y tropopausa, y arreglos (§6.15);
# 3.14.0: el chorro de la troposfera es un nucleo (maximo local en latitud y altura), chorros de invierno y
# huella del codigo en el año extra (§6.16).
# VALIDACION del modo Tierra de M3N con el nucleo dinamico (I16) frente a la Tierra real (DISENO_FASE_6.3.md §6.14).
# Diagnostico: no cambia nada del modelo ni de la simulacion larga.
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
#       (el resto), y la parte latente, frente a Trenberth y Caron (2001); el transporte a traves del ecuador
#       frente a Donohoe et al. (2013).
#   PARTE 2 (un año mas, ~15 min en el PC de Carlos) -- desde una COPIA del ultimo punto de control (la
#     simulacion larga no se toca), con todos los diagnosticos de ese año:
#     - balance de radiacion en el tope y en la superficie frente a CERES EBAF (Loeb et al. 2018), con cielo
#       despejado y con nubes, y Wild et al. (2019): M3N NO TIENE NUBES (Fase 5.2), asi que lo comparable es el
#       cielo despejado; la diferencia con el cielo real es el efecto de las nubes;
#     - hielo marino (extension maxima y minima) frente a NSIDC;
#     - perfil vertical por bandas (tropopausa de la OMM en un nivel e interpolada, punto mas frio, T a 25 km,
#       gradiente 0-6 km) frente a ERA-Interim;
#     - VIENTOS del nucleo: chorros de la troposfera y de la estratosfera, vientos en la capa baja (alisios y
#       del oeste) y celulas de Hadley (funcion de corriente de masa), anual, diciembre-febrero y junio-agosto,
#       frente a ERA-Interim 1979-2016 (3.13.0), medido con las mismas funciones en la misma rejilla.
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


# -------------------------------------------------------- metricas (3.13.0): las mismas para M3N y ERA-Interim
# Cada funcion recibe campos ya en la rejilla de 5 grados de M3N y una coordenada vertical normalizada
# ("nivel": sigma en M3N; p / 1000 hPa en ERA-Interim), para que el modelo y la referencia se midan igual
# (herramientas/referencias_era_interim.py calcula con ellas los valores de referencias_tierra.REF_ERAI).
NIVEL_HADLEY = (0.15, 0.95)           # interior de la columna donde se busca el maximo de cada celula
# Cada celula de Hadley se busca en SU lado (3.13.0; la 3.12.0 buscaba las dos en 40 S-40 N y en diciembre-
# febrero cogio como "celula del sur" la de Ferrel del norte, que tiene su mismo signo): la del norte entre
# 20 S y 40 N (su maximo cruza el ecuador en el invierno del norte), la del sur entre 40 S y 20 N. Asi la
# celula de Ferrel del otro hemisferio (por encima de ~25-30 grados) queda siempre fuera.
VENTANA_HADLEY = {"norte": (-20.0, 40.0), "sur": (-40.0, 20.0)}
NIVEL_TROPOSFERA = (0.1, 0.5)         # chorro de la troposfera: ~100-500 hPa
NIVEL_ESTRATOSFERA = 0.0101           # viento de la estratosfera alta: ~10 hPa y por encima
# 3.14.0 (§6.16): el chorro de la troposfera es un NUCLEO, un maximo local del viento a la vez en la vertical y en
# latitud dentro de la ventana de niveles (como Manney et al. 2011, ACP 11, 6115, §3.1: maximos en 2-D entre 400
# y 100 hPa). Ademas el viento tiene que BAJAR en el nivel de encima: si sigue subiendo, ese punto es la parte de
# abajo del chorro de la estratosfera, no un chorro de la troposfera (el mismo criterio con el que Manney et al.
# 2011, §3.2, separan el chorro del vortice polar). En la 3.13.0 el chorro era el maximo de la ventana sin mas, y
# en M3N, cuyo viento crece con la altura hasta el tope en latitudes medias, salia siempre en el borde de arriba
# (~112 hPa) y a 52,5 grados: era la cola del chorro de la estratosfera. Si no hay nucleo, se da NaN ("sin
# nucleo"). Los umbrales de 40 y 30 m/s de Manney et al. son para campos instantaneos y no se usan aqui (medias
# zonales de una temporada o un año).
TEMPORADAS_INVIERNO = (("dic-feb", "norte"), ("jun-ago", "sur"))   # chorros de invierno de cada hemisferio
VENTANAS_VIENTO_BAJO = (("Alisios norte (0-30 N): minimo", "alisios_n", (0, 30), np.argmin),
                        ("Alisios sur (0-30 S): minimo", "alisios_s", (-30, 0), np.argmin),
                        ("Oeste norte (30-70 N): maximo", "oeste_n", (30, 70), np.argmax),
                        ("Oeste sur (30-70 S): maximo", "oeste_s", (-70, -30), np.argmax))


def celulas_hadley(psi, lat_c, nivel, k5):
    """psi (niveles, caras), caras de norte a sur; nivel: coordenada normalizada de cada nivel de psi; k5: el
    nivel mas cercano a 0,5. Devuelve {"norte": (maximo, lat, borde), "sur": (minimo, lat, borde)}; el borde es
    el cruce por cero de psi en k5 yendo hacia el polo desde el maximo (como TropD, Adam et al. 2018)."""
    nivel = np.asarray(nivel)
    inter = (nivel > NIVEL_HADLEY[0]) & (nivel < NIVEL_HADLEY[1])
    out = {}
    for clave, f, paso in (("norte", np.nanargmax, -1), ("sur", np.nanargmin, +1)):   # lat_c de norte a sur
        a, b = VENTANA_HADLEY[clave]
        m = (lat_c > a) & (lat_c < b)
        sub = np.where(inter[:, None] & m[None, :], psi, np.nan)
        k, j = np.unravel_index(f(sub), sub.shape)
        out[clave] = (float(psi[k, j]), float(lat_c[j]), cruce_hacia_el_polo(psi[k5], lat_c, j, paso))
    return out


def nucleo_troposfera(u, lat, nivel, m):
    """El nucleo del chorro de la troposfera (3.14.0, ver arriba) en las filas de la mascara m: el mayor de los
    puntos de la ventana de niveles en los que u es mayor que en el nivel de encima, no menor que en el de
    debajo y no menor que en las filas vecinas. Niveles de arriba (indice 0) a abajo. Sin nucleo: NaN."""
    nk, ni = u.shape
    mejor = (np.nan, np.nan, np.nan)
    for k in range(1, nk - 1):
        if not (NIVEL_TROPOSFERA[0] <= nivel[k] <= NIVEL_TROPOSFERA[1]):
            continue
        for i in np.where(m)[0]:
            c = u[k, i]
            if not (c > u[k - 1, i] and c >= u[k + 1, i]):
                continue
            if (i > 0 and u[k, i - 1] > c) or (i < ni - 1 and u[k, i + 1] > c):
                continue
            if np.isnan(mejor[0]) or c > mejor[0]:
                mejor = (float(c), float(lat[i]), float(nivel[k]))
    return mejor


def chorros(u, lat, nivel):
    """u (capas, filas) media zonal, de arriba abajo; nivel de cada capa. Entre 15 y 70 grados de cada
    hemisferio: el nucleo del chorro de la troposfera (nivel 0,1-0,5; NaN si no hay) y el maximo del viento del
    oeste en la estratosfera alta (nivel <= 0,01). Devuelve {(capa, hemisferio): (m/s, lat, nivel)}."""
    nivel = np.asarray(nivel)
    u = np.asarray(u, dtype=float)
    out = {}
    hemis = (("norte", (lat > 15) & (lat < 70)), ("sur", (lat < -15) & (lat > -70)))
    for hemi, m in hemis:
        out[("troposfera", hemi)] = nucleo_troposfera(u, lat, nivel, m)
    sel = nivel <= NIVEL_ESTRATOSFERA
    for hemi, m in hemis:
        sub = np.where(sel[:, None] & m[None, :], u, -np.inf)
        k, i = np.unravel_index(np.argmax(sub), sub.shape)
        out[("estratosfera", hemi)] = (float(u[k, i]), float(lat[i]), float(nivel[k]))
    return out


def texto_chorro(t, con_nivel=True):
    """(m/s, lat, nivel) -> texto; 'sin nucleo' si es NaN."""
    v, la, niv = t
    if v is None or not np.isfinite(v):
        return "sin nucleo".ljust(25 if con_nivel else 16)
    return f"{v:5.1f} m/s a {la:+5.1f}, {1000 * niv:4.0f} hPa" if con_nivel else f"{v:5.1f} a {la:+5.1f}"


def vientos_bajos(us, lat):
    """us (filas): viento zonal de la capa baja. Devuelve {clave: (m/s, lat)} (alisios y vientos del oeste)."""
    out = {}
    for _, clave, (a, b), f in VENTANAS_VIENTO_BAJO:
        m = (lat > a) & (lat < b)
        i = np.where(m)[0][f(us[m])]
        out[clave] = (float(us[i]), float(lat[i]))
    return out


def tropopausa_interpolada(Tz, zz, pz):
    """Tropopausa de la OMM con el cruce de 2 K/km interpolado entre los centros de las capas (a la manera
    de Reichler et al. 2003, GRL 30, 2042, aqui en altura): con ~2 km entre niveles cerca de la tropopausa,
    la de la OMM en un nivel (validar_v30.tropopausa_omm) salta de nivel en nivel. Devuelve (z, T)."""
    import validar_v30 as V30
    zt, Tt = V30.tropopausa_omm(Tz, zz, pz)
    if np.isnan(zt):
        return zt, Tt
    i = int(np.where(zz == zt)[0][0])
    gam = -np.diff(Tz) / np.diff(zz) * 1000
    if i == 0 or gam[i - 1] <= 2.0:
        return float(zt), float(Tt)
    zm = 0.5 * (zz[1:] + zz[:-1])
    z = zm[i - 1] + (gam[i - 1] - 2.0) / (gam[i - 1] - gam[i]) * (zm[i] - zm[i - 1])
    return float(z), float(np.interp(z, zz, Tz))


def perfil_banda(Tz, zz, pz):
    """Perfil medio de una banda (de abajo arriba; K, m, Pa): tropopausa de la OMM en un nivel y
    interpolada, punto mas frio entre 500 y 50 hPa ("cold point") y T a 25 km. Temperaturas en C."""
    import validar_v30 as V30
    zt, Tt = V30.tropopausa_omm(Tz, zz, pz)
    zi, Ti = tropopausa_interpolada(Tz, zz, pz)
    sel = np.where((pz >= 5e3) & (pz <= 5e4))[0]
    i = sel[np.argmin(Tz[sel])]
    T25 = float(np.interp(25000, zz, Tz)) if zz[-1] > 25000 else np.nan
    return {"tropopausa_km": zt / 1000, "T_tropopausa": Tt - 273.15, "tropopausa_interp_km": zi / 1000,
            "T_tropopausa_interp": Ti - 273.15, "frio_km": zz[i] / 1000, "frio_hPa": pz[i] / 100,
            "T_frio": Tz[i] - 273.15, "T_25km": T25 - 273.15}


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
    from cache_simulacion import huella_codigo, MODULOS_I16
    codigo = huella_codigo(MODULOS_I16)      # 3.14.0: el año extra se rehace tambien si cambia el codigo
    f_guardado = os.path.join(SALIDA, "ano_extra.pkl")
    if os.path.exists(f_guardado):
        with open(f_guardado, "rb") as f:
            d = pickle.load(f)
        if d["huella"] == huella and d.get("codigo") == codigo:
            print("Año extra: ya estaba hecho con este punto de control y este codigo (se reutiliza).", flush=True)
            return d["r"], d["vientos"]
        print("Año extra: el guardado es de otro punto de control o de otro codigo"
              + (" (o anterior a la 3.14.0, sin huella del codigo)" if "codigo" not in d else "") + ": se rehace.",
              flush=True)
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
        pickle.dump({"huella": huella, "codigo": codigo, "r": r, "vientos": vientos}, f)
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
    i0 = int(np.argmin(np.abs(bordes)))
    L.append(f"   A traves del ecuador (>0 = hacia el norte): atmosfera {atm[i0]:+5.2f} PW, oceano {oc[i0]:+5.2f} PW."
             f"   obs. atmosfera {R.REF_TRANSPORTE_ECUADOR['atm_PW']:+.1f} PW, banda de lluvias en "
             f"{R.REF_TRANSPORTE_ECUADOR['itcz_centroide']:.2f} N (Donohoe et al. 2013)")
    L.append(f"   (total desde el balance en el tope; oceano desde la convergencia de su difusion; atmosfera = total -"
             f" oceano; latente = L_v (P - E). Desequilibrio global quitado antes de integrar: tope {-resto_toa:+.3f},"
             f" oceano {resto_oc:+.3f} W/m2)")
    L.append("   Por bordes: lat | total | atmosfera | oceano | latente")
    for i in range(len(bordes)):
        if abs(bordes[i]) % 10 == 5 or abs(bordes[i]) % 10 == 0:
            L.append(f"     {bordes[i]:+5.0f} | {tot[i]:+6.2f} | {atm[i]:+6.2f} | {oc[i]:+6.2f} | {lat_t[i]:+6.2f}")
    L.append("")
    L.append("   Para comparar, la v3.1 con la MISMA fisica pero SIN nucleo (transporte por difusion; DISENO_FASE_5.1.md"
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
             " M3N no tiene nubes (Fase 5.2).")
    L.append("   (CERES EBAF Ed4.0, jul 2005-jun 2015, Loeb et al. 2018, tabla 5. El cielo despejado de CERES no esta en"
             " equilibrio: su neto es lo que las nubes le quitan a la Tierra real)")
    L.append("")
    h = V30.hielo(r)
    L.append("5. HIELO MARINO, extension (millones de km2; celda entera si hay hielo)   M3N   observado (NSIDC 1981-2010)")
    for k, nom in (("N_max", "Norte, maximo"), ("N_min", "Norte, minimo"), ("S_max", "Sur, maximo"), ("S_min", "Sur, minimo")):
        L.append(f"   {nom:16s} {h[k]:6.1f}   {R.REF_HIELO[k]:5.1f}")
    L.append("")
    L.append("6. PERFIL VERTICAL (media anual; M3N: celdas con el suelo por debajo de ~500 m; obs.: ERA-Interim 1979-2016)")
    L.append("   Banda      | tropopausa OMM: en un nivel | interpolada  | T tropop. interp. | mas frio (500-50 hPa) | T a 25 km | gradiente 0-6 km")
    E = R.REF_ERAI["tropopausa"]
    for b in V30.perfiles_por_bandas(fl, P.P3N_GRAVEDAD):
        if b is None:
            continue
        banda, Tz, zz, pz = b
        d = perfil_banda(Tz, zz, pz)
        g06 = float((Tz[0] - np.interp(6000, zz, Tz)) / ((6000 - zz[0]) / 1000))
        o = E[banda]
        tropico = banda in ("30..10", "10..-10", "-10..-30")
        frio_m = f"{d['T_frio']:6.1f} C a {d['frio_km']:4.1f} km" if tropico else "        --        "
        frio_o = f"{o['T_frio']:6.1f} C a {o['frio_km']:4.1f} km" if tropico else "        --        "
        L.append(f"   {banda:10s} M3N  | {d['tropopausa_km']:5.1f} km                 | {d['tropopausa_interp_km']:5.1f} km     |"
                 f"  {d['T_tropopausa_interp']:6.1f} C         | {frio_m}  | {d['T_25km']:6.1f} C  |  {g06:4.2f} K/km")
        L.append(f"   {'':10s} obs. | {o['tropopausa_km']:5.1f} km                 | {o['tropopausa_interp_km']:5.1f} km     |"
                 f"  {o['T_tropopausa_interp']:6.1f} C         | {frio_o}  | {o['T_25km']:6.1f} C  |     --")
    L.append("   (OMM: primer nivel con gradiente <= 2 K/km que se mantiene 2 km; M3N tiene ~2 km entre niveles cerca de la")
    L.append("   tropopausa, por eso tambien la interpolada (Reichler et al. 2003). El punto mas frio solo en los tropicos.")
    L.append("   Comprobacion independiente de ERA-Interim: Seidel et al. 2001, radiosondas, ecuador: ~16,5 km y ~-81 C)")
    L.append("")
    if vientos is None or "anual" not in vientos:
        return
    sig = vientos["sigma"]
    sm = 0.5 * (sig[1:] + sig[:-1])
    lat_c = vientos["lat_caras"]
    L.append("7. VIENTOS DEL NUCLEO (obs.: ERA-Interim 1979-2016 en la rejilla de 5 grados de M3N, medido igual)")
    u = vientos["anual"]["u"]                                         # (N, FILAS) media zonal
    ch, cho = chorros(u, LAT, sm), R.REF_ERAI["chorros"]
    L.append("   Chorros, media anual (troposfera: nucleo, maximo local en latitud y altura; 'sin nucleo' si el viento")
    L.append("   sigue creciendo hacia la estratosfera. Estratosfera: maximo por encima de ~10 hPa):")
    for (capa, hemi), t in ch.items():
        nom = f"{'Troposfera' if capa == 'troposfera' else 'Estratosfera'} ({hemi})"
        L.append(f"     {nom:22s} M3N {texto_chorro(t)}   obs. {texto_chorro(cho[f'{capa}_{hemi}'])}")
    cht = R.REF_ERAI["chorros_temporada"]
    L.append("   Chorros de invierno de cada hemisferio (dic-feb el norte, jun-ago el sur) y estratosfera de verano:")
    for temp, hemi in TEMPORADAS_INVIERNO:
        if temp not in vientos:
            continue
        verano = "sur" if hemi == "norte" else "norte"
        c = chorros(vientos[temp]["u"], LAT, sm)
        for capa, h, nom in (("troposfera", hemi, f"Troposfera ({hemi}, {temp})"),
                             ("estratosfera", hemi, f"Estratosfera ({hemi}, {temp})"),
                             ("estratosfera", verano, f"Estratosfera ({verano}, verano)")):
            L.append(f"     {nom:30s} M3N {texto_chorro(c[(capa, h)])}   obs. {texto_chorro(cht[temp][f'{capa}_{h}'])}")
    L.append("   (en verano la estratosfera real tiene viento del ESTE, negativo, por el calentamiento del ozono)")
    vb, vbo, vb10 = vientos_bajos(u[-1], LAT), R.REF_ERAI["vientos_bajos_975hPa"], R.REF_ERAI["vientos_bajos_10m"]
    for nom, clave, _, _ in VENTANAS_VIENTO_BAJO:
        L.append(f"   Capa baja, {nom:32s} M3N {vb[clave][0]:+5.1f} m/s a {vb[clave][1]:+5.1f}   "
                 f"obs. 975 hPa {vbo[clave][0]:+5.1f} a {vbo[clave][1]:+5.1f} (10 m: {vb10[clave][0]:+5.1f})")
    L.append("   Celulas de Hadley (funcion de corriente de masa, 10^10 kg/s; >0 = la del norte, <0 = la del sur;")
    L.append("   cada una buscada en su lado; borde = cruce por cero a sigma 0,5 / 500 hPa):")
    k5 = int(np.argmin(np.abs(sig - 0.5)))                            # interfaz mas cercana a sigma 0,5
    for clave in ("anual", "dic-feb", "jun-ago"):
        if clave not in vientos:
            continue
        psi = funcion_corriente(vientos[clave]["vps"], sig, vientos["cos_caras"], P.P3N_RADIO, P.P3N_GRAVEDAD) / 1e10
        campos["psi_" + clave] = psi
        campos["u_" + clave] = vientos[clave]["u"]
        c, co = celulas_hadley(psi, lat_c, sig, k5), R.REF_ERAI["hadley"][clave]
        for cel in ("norte", "sur"):
            L.append(f"     {clave:8s} {cel:5s}  M3N {c[cel][0]:+6.2f} a {c[cel][1]:+5.1f}, borde {c[cel][2]:+5.1f}"
                     f"   obs. {co[cel][0]:+6.2f} a {co[cel][1]:+5.1f}, borde {co[cel][2]:+5.1f}")
    L.append("   Viento zonal medio anual por filas (m/s): lat | sigma ~0,25: M3N  obs. | ~0,5: M3N  obs. | capa baja: M3N  obs.")
    k25, k50 = int(np.argmin(np.abs(sm - 0.25))), int(np.argmin(np.abs(sm - 0.5)))
    ub = R.REF_ERAI["u_bandas"]
    for i in range(0, FILAS, 2):
        L.append(f"     {LAT[i]:+5.1f} | {u[k25, i]:+6.1f} {ub[250][i]:+6.1f} | {u[k50, i]:+6.1f} {ub[500][i]:+6.1f} |"
                 f" {u[-1, i]:+6.1f} {ub[975][i]:+6.1f}")
    L.append("   (obs. de la capa baja: ERA-Interim a 975 hPa; por debajo de ~65 S ese nivel queda bajo el suelo antartico)")
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
    from cache_simulacion import huella_codigo, MODULOS_I16
    tramos = rc.get("climatologia", {}).get("procedencia") or []
    L.append(f"codigo de esta validacion: {huella_codigo(MODULOS_I16)} | simulacion: "
             + (", ".join(f"M3N {t['version']} ({t['codigo']}) desde el año {t['ano_sim']}" for t in tramos)
                if tramos else "procedencia no consta (anterior a la 3.14.0)"))
    L.append("")
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

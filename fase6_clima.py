# fase6_clima.py -- v3.10.0: clima de M3N con el nucleo dinamico (I16): arranque caliente, equilibrio y
# climatologia de N años. Diseno: DISENO_FASE_6.3.md §6.3 (decisiones de Carlos del 06/10) y §6.10.
#
#   1. Equilibrio del modelo de 2 capas (v2.4.3) con el mismo mapa y los mismos parametros: su oceano, su
#      suelo y su hielo son el ARRANQUE CALIENTE (el aire arranca en reposo).
#   2. Modelo con I16 desde ese estado hasta el equilibrio (fase6_equilibrio.py: medias globales anuales en
#      una ventana de 5 años).
#   3. Años REGISTRADOS, uno tras otro, hasta que el error de la media de la climatologia cumple el criterio
#      en todas las bandas de latitud (o se llega al maximo de años): media de N años por dia del año,
#      extremos absolutos, desviacion entre años, medias horarias.
#
# Todo con PUNTO DE CONTROL en una carpeta: si se corta (Ctrl+C, apagon...), la misma orden continua donde iba
# y el resultado es el mismo bit a bit.
#
# Criterio del error de la media (§6.10.3; detalles aprobados por Carlos el 07/10):
#   - banda = cada fila de la rejilla (media zonal de la media anual);
#   - error = sigma / sqrt(N_eff), N_eff = N (1 - r1) / (1 + r1), r1 = autocorrelacion de un año al siguiente
#     (acotada a >= 0): la aproximacion para una serie AR(1) de Wilks (2011, Statistical Methods in the
#     Atmospheric Sciences, 3.a ed.) ✅; Zwiers y von Storch (1995, J. Climate) dan una version mas exacta;
#   - aire a 2 m: error < 0,1 K; precipitacion: error < 5 % o < 0,05 mm/dia (bandas casi secas);
#   - minimo 5 años registrados, maximo 30.

import hashlib
import json
import os
import time
import pickle
import shutil

import numpy as np

import fase2b_atmosfera as F

MIN_ANOS = 5
MAX_ANOS = 30
ERROR_T = 0.1             # K
ERROR_P_REL = 0.05
ERROR_P_ABS = 0.05        # mm/dia
MAX_ANOS_EQUILIBRIO = 50

# registros diarios del resultado de simular_fase2b: (clave, como se combinan los extremos)
DIARIOS = (("reg_media", None), ("reg_min", "min"), ("reg_max", "max"), ("suelo_media", None),
           ("suelo_min", "min"), ("suelo_max", "max"), ("cl_media", None), ("tr_media", None),
           ("hielo_espesor", None))
HORARIOS = ("horario_aire2m", "horario_superficie")


def autocorrelacion_1(y):
    """Autocorrelacion de retardo 1 (estimador habitual, con la varianza de toda la serie). y: (N, ...)."""
    y = np.asarray(y, dtype=float)
    a = y - y.mean(axis=0)
    den = (a * a).sum(axis=0)
    num = (a[1:] * a[:-1]).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)


def error_media(y):
    """Error de la media de N años (por banda): sigma / sqrt(N_eff), N_eff = N (1 - r1)/(1 + r1), r1 >= 0."""
    y = np.asarray(y, dtype=float)
    n = y.shape[0]
    if n < 2:
        return np.full(y.shape[1:], np.inf)
    sigma = y.std(axis=0, ddof=1)
    r1 = np.maximum(autocorrelacion_1(y), 0.0)
    n_eff = n * (1 - r1) / (1 + r1)
    return sigma / np.sqrt(np.maximum(n_eff, 1.0))


def criterio_climatologia(bandas_T, bandas_P, min_anos=MIN_ANOS):
    """bandas_T, bandas_P: (N años, filas). Devuelve (cumple, info)."""
    n = len(bandas_T)
    eT = error_media(bandas_T)
    eP = error_media(bandas_P)
    Pm = np.mean(bandas_P, axis=0)
    okT = eT < ERROR_T
    okP = (eP < ERROR_P_REL * Pm) | (eP < ERROR_P_ABS)
    info = {"anos": n, "error_T_max_K": float(np.max(eT)), "error_P_max_mm_dia": float(np.max(eP)),
            "error_P_rel_max": float(np.max(np.where(Pm > 0, eP / np.where(Pm > 0, Pm, 1.0), 0.0))),
            "bandas_T_sin_cumplir": int((~okT).sum()), "bandas_P_sin_cumplir": int((~okP).sum())}
    return bool(n >= min_anos and okT.all() and okP.all()), info


class Acumulador:
    """Sumas, minimos, maximos y sumas de cuadrados de los años registrados (para el punto de control)."""

    def __init__(self):
        self.n = 0
        self.suma, self.suma2, self.ext = {}, {}, {}
        self.ref = {}              # el primer año: las sumas de cuadrados van respecto a el (sin cancelacion)
        self.bandas_T, self.bandas_P = [], []
        self.por_ano = []           # diagnosticos pequeños de cada año (energia, agua)

    def anadir(self, r):
        self.n += 1
        for k, modo in DIARIOS:
            x = np.asarray(r[k], dtype=float)
            self._sumar(k, x)
            if modo is not None:
                f = np.minimum if modo == "min" else np.maximum
                self.ext[k] = x.copy() if k not in self.ext else f(self.ext[k], x)
        self._cuadrado("reg_media", np.asarray(r["reg_media"], float))
        for k in HORARIOS:
            if r.get(k) is not None:
                self._sumar(k, np.asarray(r[k], dtype=float))
        a = r["agua"]
        for k in ("precipitacion", "nieve", "evaporacion", "escorrentia", "agua_precipitable",
                  "diario_precipitacion", "diario_nieve", "diario_evaporacion"):
            self._sumar("agua_" + k, np.asarray(a[k], dtype=float))
        self._cuadrado("agua_precipitacion", np.asarray(a["precipitacion"], float))
        for k in ("toa_neto", "T_atm", "conv_oceano", "conv_atmosfera"):
            self._sumar("flujos_" + k, np.asarray(r["flujos"][k], dtype=float))
        self.bandas_T.append(np.asarray(r["reg_media"], float).mean(axis=(0, 2)))       # C, media zonal anual
        self.bandas_P.append(np.asarray(a["precipitacion"], float).mean(axis=1))        # mm/dia
        self.por_ano.append({"energia": {k: v for k, v in r["energia"].items() if np.ndim(v) == 0},
                             "cierre_agua_atmosfera_kg_m2": a["cierre_agua_atmosfera_kg_m2"],
                             "cierre_agua_tierra_kg_m2": a["cierre_agua_tierra_kg_m2"]})

    def _cuadrado(self, k, x):
        if k not in self.ref:
            self.ref[k] = x.copy()
        d = x - self.ref[k]
        self.suma2[k] = (d * d) if k not in self.suma2 else self.suma2[k] + d * d
        self.suma2[k + "_d"] = d.copy() if k + "_d" not in self.suma2 else self.suma2[k + "_d"] + d

    def _sumar(self, k, x):
        self.suma[k] = x.copy() if k not in self.suma else self.suma[k] + x


def _guardar(obj, ruta):
    tmp = f"{ruta}.{os.getpid()}.tmp"
    with open(tmp, "wb") as f:
        pickle.dump(obj, f)
    os.replace(tmp, ruta)


def _huella(*cosas):
    return hashlib.sha256(pickle.dumps(cosas)).hexdigest()


def arranque_caliente(argumentos, carpeta, informar=print, opciones=None, exigir_equilibrio=True, nombre="dos_capas"):
    """Equilibrio del modelo de 2 capas (interruptores por defecto, con aceleracion) con los mismos
    argumentos: su estado final, guardado en la carpeta (y comprobado con una huella al reutilizarlo).
    Con opciones={"interruptores": ...} y otro nombre, el de otra configuracion (p. ej., la v3.1 sin nucleo)."""
    ruta = os.path.join(carpeta, f"arranque_{nombre}.pkl")
    huella = _huella(*[np.asarray(a).tobytes() if isinstance(a, np.ndarray) else repr(a) for a in argumentos[1:]],
                     len(argumentos[0]), float(argumentos[0][0][0]), float(argumentos[0][-1][0]),
                     sorted((opciones or {}).items()))
    if os.path.exists(ruta):
        with open(ruta, "rb") as f:
            d = pickle.load(f)
        if d["huella"] == huella:
            return d["estado"]
        raise ValueError(f"{ruta} es de otra simulacion (otro mapa o parametros); borralo o usa otra carpeta")
    informar(f"1/3 Equilibrio del modelo de partida ({nombre}) para el arranque caliente...")
    r = F.simular_fase2b(*argumentos, **(opciones or {}))
    informar(f"    {r['anos']} años, convergido: {r['convergido']}")
    if exigir_equilibrio and not r["convergido"]:
        raise RuntimeError(f"el modelo de partida ({nombre}) no ha llegado al equilibrio: no sirve de arranque caliente")
    _guardar({"huella": huella, "estado": r["estado_final"]}, ruta)
    return r["estado_final"]


def leer_procedencia(carpeta):
    """La lista de tramos de procedencia.json de una simulacion larga ([] si no hay: anterior a la 3.14.0)."""
    ruta = os.path.join(carpeta, "procedencia.json")
    if not os.path.exists(ruta):
        return []
    with open(ruta, encoding="utf-8") as f:
        return json.load(f)


def registrar_procedencia(carpeta, version, codigo, ano_sim, anos_registrados, informar=print):
    """3.14.0 (DISENO_FASE_6.3.md §6.16): TRAZABILIDAD de una simulacion larga. Cada vez que empieza o continua,
    añade a <carpeta>/procedencia.json un tramo con la fecha, la version de M3N, la huella del codigo
    (cache_simulacion.MODULOS_I16), el año simulado y los años registrados en ese momento. Si el codigo no es el
    del tramo anterior, AVISA (no para: continuar con otro codigo puede ser lo que se quiere, pero tiene que
    constar). Devuelve la lista de tramos. No cambia nada de la simulacion."""
    tramos = leer_procedencia(carpeta)
    if tramos and tramos[-1]["codigo"] != codigo:
        informar(f"    AVISO: el codigo de M3N ha cambiado desde el ultimo tramo de esta simulacion (version "
                 f"{tramos[-1]['version']}, huella {tramos[-1]['codigo']}; ahora {version}, {codigo}). El resultado "
                 "mezclara versiones; queda registrado en procedencia.json y en el resumen.")
    elif not tramos and (ano_sim or anos_registrados):
        informar("    AVISO: esta simulacion empezo antes de la 3.14.0: no consta con que codigo se hicieron los "
                 f"{ano_sim} años anteriores.")
    tramos.append({"fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "version": version, "codigo": codigo,
                   "ano_sim": int(ano_sim), "anos_registrados": int(anos_registrados)})
    _escribir_json(tramos, os.path.join(carpeta, "procedencia.json"))
    return tramos


def _escribir_json(obj, ruta):
    with open(ruta + ".tmp", "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False)
    os.replace(ruta + ".tmp", ruta)


def simular_clima(argumentos, interruptores, carpeta, min_anos=MIN_ANOS, max_anos=MAX_ANOS,
                  max_anos_equilibrio=MAX_ANOS_EQUILIBRIO, informar=print, opciones_dos_capas=None,
                  exigir_equilibrio_dos_capas=True, estado_inicial=None, arranque="dos_capas", version=None):
    """argumentos: los 8 posicionales de simular_fase2b (orbita, tipo, altitud, emisividad, albedos,
    inercias, profundidad optica, D). interruptores: con I16 encendido. Devuelve un dict como el de
    simular_fase2b con las MEDIAS de los años registrados, mas 'extremos', 'desviacion_entre_anos' y
    'climatologia' (años promediados, errores, criterio, equilibrio, diagnosticos de cada año y, desde la
    3.14.0, 'procedencia': los tramos de procedencia.json). version: la de M3N, para la procedencia."""
    if not interruptores.get("nucleo_dinamico", False):
        raise ValueError("simular_clima es para el modelo con nucleo dinamico (I16)")
    os.makedirs(carpeta, exist_ok=True)
    # estado_inicial (dict con T_col y HIELO): para arrancar desde otro estado en vez del modelo de 2 capas
    # arranque="dos_capas" (decidido, §6.3): oceano, suelo y hielo del equilibrio del modelo de 2 capas.
    # arranque="v31" (alternativa 🔶, §6.10): los de la v3.1 con los mismos interruptores pero SIN nucleo, por si
    # la superficie del modelo de 2 capas resulta demasiado lejos del equilibrio del modelo con nucleo.
    if estado_inicial is not None:
        estado_ini = estado_inicial
    elif arranque == "dos_capas":
        estado_ini = arranque_caliente(argumentos, carpeta, informar, opciones_dos_capas, exigir_equilibrio_dos_capas)
    elif arranque == "v31":
        sin_nucleo = dict(interruptores, nucleo_dinamico=False)
        estado_ini = arranque_caliente(argumentos, carpeta, informar, dict(opciones_dos_capas or {}, interruptores=sin_nucleo),
                                       exigir_equilibrio_dos_capas, nombre="v31")
    else:
        raise ValueError(f"arranque desconocido: {arranque}")
    f_sim = os.path.join(carpeta, "estado_i16.pkl")
    f_previo = f_sim + ".previo"
    f_acum = os.path.join(carpeta, "climatologia.pkl")
    kw = dict(interruptores=interruptores, estado_inicial=estado_ini, archivo_estado=f_sim, guardar_al_terminar=True)

    def leer(ruta):
        with open(ruta, "rb") as f:
            return pickle.load(f)

    acum, equilibrio, ano_esperado = Acumulador(), None, None
    if os.path.exists(f_acum):
        d = leer(f_acum)
        acum, equilibrio, ano_esperado = d["acum"], d["equilibrio"], d["ano_sim"]
        informar(f"-> continua: {acum.n} años registrados")
        cumple, info = criterio_climatologia(acum.bandas_T, acum.bandas_P, min_anos)
        if cumple or acum.n >= max_anos:               # ya estaba terminada: no se simula nada mas
            return combinar(acum, equilibrio, info, cumple, leer_procedencia(carpeta))
    from cache_simulacion import huella_codigo, MODULOS_I16
    ck0 = leer(f_sim) if os.path.exists(f_sim) else None
    tramos = registrar_procedencia(carpeta, version, huella_codigo(MODULOS_I16), ck0["ano"] if ck0 else 0, acum.n,
                                   informar)

    while True:
        # ¿desde donde sigue el simulador? (un corte entre el final de un año registrado y el guardado de la
        # climatologia se deshace volviendo al punto de control previo: ese año se repite, identico)
        ck = leer(f_sim) if os.path.exists(f_sim) else None
        if ck is not None and ck["registrado"] and ck["ano"] != ano_esperado:
            previo = leer(f_previo)
            shutil.copyfile(f_previo, f_sim)
            ck = previo
        if ano_esperado is not None:
            max_anos_llamada = ano_esperado                       # el siguiente año registrado
        elif ck is None:
            max_anos_llamada = max_anos_equilibrio                # empezar
        elif ck["convergido"] or ck["registrado"]:
            max_anos_llamada = ck["ano"]                          # ya en equilibrio: directo al año registrado
        else:
            max_anos_llamada = max_anos_equilibrio                # seguir buscando el equilibrio
        if equilibrio is None:
            informar("2/3 Modelo con nucleo dinamico hasta el equilibrio (+ 1.er año registrado)...")
        r = F.simular_fase2b(*argumentos, max_anos=max_anos_llamada, **kw)
        ck = leer(f_sim)
        if equilibrio is None:
            # el estado del comienzo del año registrado (en el .previo) dice si se llego al equilibrio
            antes = leer(f_previo)
            ok, valores = F.EQ.evaluar(antes["serie_equilibrio"])
            equilibrio = {"convergido": bool(antes["convergido"]), "anos": int(antes["ano"]),
                          "serie": antes["serie_equilibrio"], "valores": valores}
            if not antes["convergido"]:
                informar(f"    AVISO: sin equilibrio tras {antes['ano']} años; la climatologia sera la de un estado sin equilibrar")
        ano_esperado = ck["ano"]
        acum.anadir(r)
        del r
        cumple, info = criterio_climatologia(acum.bandas_T, acum.bandas_P, min_anos)
        _guardar({"acum": acum, "equilibrio": equilibrio, "ano_sim": ano_esperado}, f_acum)
        if os.path.exists(f_previo):
            os.remove(f_previo)
        informar(f"3/3 Año registrado {acum.n}: error aire 2 m {info['error_T_max_K']:.3f} K, "
                 f"precipitacion {info['error_P_max_mm_dia']:.3f} mm/dia ({100 * info['error_P_rel_max']:.1f} %)"
                 + (" -> CUMPLE" if cumple else ""))
        if cumple or acum.n >= max_anos:
            break
    return combinar(acum, equilibrio, info, cumple, tramos)


def combinar(acum, equilibrio, info, cumple, procedencia=None):
    n = acum.n
    m = {k: v / n for k, v in acum.suma.items()}
    r = {k: m[k] for k, _ in DIARIOS}
    for k in HORARIOS:
        r[k] = m.get(k)
    r["agua"] = {k[5:]: v for k, v in m.items() if k.startswith("agua_")}
    r["flujos"] = {k[7:]: v for k, v in m.items() if k.startswith("flujos_")}
    r["extremos"] = {k: v for k, v in acum.ext.items()}
    # varianza entre años (n-1), con las desviaciones respecto al primer año: sum(d^2) - (sum d)^2 / n
    var = lambda k: np.maximum(acum.suma2[k] - acum.suma2[k + "_d"] ** 2 / n, 0.0) / max(n - 1, 1)
    r["desviacion_entre_anos"] = {"reg_media": np.sqrt(var("reg_media")),
                                  "precipitacion": np.sqrt(var("agua_precipitacion"))}
    r["climatologia"] = {"anos_promediados": n, "cumple_criterio": cumple, **info,
                         "equilibrio": equilibrio, "por_ano": acum.por_ano,
                         "bandas_aire2m_C": np.array(acum.bandas_T), "bandas_precipitacion": np.array(acum.bandas_P),
                         "procedencia": list(procedencia or [])}
    return r

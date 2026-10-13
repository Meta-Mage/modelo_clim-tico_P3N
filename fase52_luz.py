# fase52_luz.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §3.2, §3.3 y §7): LECTORES Y PRIMEROS CALCULOS
# de los datos de la luz que descarga herramientas/descargar_luz.py (en el PC de Carlos). No se usa todavia en el
# modelo.
#   - El espectro de S3N con los modelos PHOENIX (Husser et al. 2013): lector propio de FITS (sin astropy: solo hace
#     falta el primer bloque de datos, un vector), interpolacion en la temperatura de la estrella y el reparto de su
#     luz entre las bandas de la luz de M3N.
#   - Las secciones eficaces del ozono de Serdyuchenko et al. 2014.
#   - La absorcion inducida por colisiones (CIA) en el formato de HITRAN (Richard et al. 2012, J. Quant. Spectrosc.
#     Radiat. Transf. 113, 1276 ⚠️: cabecera de 100 caracteres por bloque y despues "numero de onda, coeficiente").

import numpy as np

C2_UM_K = 14387.769          # h c / k en µm K (CODATA 2018)


# ---------------------------------------------------------------------------------------------------- FITS minimo
def leer_fits_1d(ruta):
    """El vector del primer bloque (HDU primario) de un FITS: cabecera en bloques de 2880 bytes de tarjetas de 80
    caracteres hasta END; datos en big-endian con BITPIX (-32, -64, 16 o 32), con BSCALE y BZERO si los hay."""
    with open(ruta, "rb") as f:
        tarjetas = {}
        fin = False
        while not fin:
            bloque = f.read(2880)
            if len(bloque) < 2880:
                raise ValueError(f"{ruta}: cabecera FITS incompleta")
            for i in range(0, 2880, 80):
                t = bloque[i:i + 80].decode("ascii", errors="replace")
                clave = t[:8].strip()
                if clave == "END":
                    fin = True
                    break
                if t[8:10] == "= ":
                    valor = t[10:].split("/")[0].strip()
                    tarjetas[clave] = valor
        bitpix = int(tarjetas["BITPIX"])
        naxis = int(tarjetas["NAXIS"])
        n = 1
        for k in range(1, naxis + 1):
            n *= int(tarjetas[f"NAXIS{k}"])
        tipo = {-32: ">f4", -64: ">f8", 16: ">i2", 32: ">i4"}[bitpix]
        datos = np.frombuffer(f.read(n * abs(bitpix) // 8), dtype=tipo).astype(float)
    if len(datos) != n:
        raise ValueError(f"{ruta}: faltan datos")
    escala = float(tarjetas.get("BSCALE", 1.0))
    cero = float(tarjetas.get("BZERO", 0.0))
    return datos * escala + cero


# ---------------------------------------------------------------------------------------------------- la estrella
def espectro_phoenix(ruta_onda, modelos, T):
    """(longitud de onda en µm, flujo) a la temperatura T interpolando el logaritmo del flujo linealmente en la
    temperatura entre los dos modelos que la rodean. modelos: {Teff: ruta del FITS}. La longitud de onda de PHOENIX
    viene en Å (vacio) y el flujo en erg/s/cm^2/cm; aqui solo importa la FORMA (se normaliza despues)."""
    lam = leer_fits_1d(ruta_onda) * 1e-4
    ts = sorted(modelos)
    if not ts[0] <= T <= ts[-1]:
        raise ValueError(f"T = {T} fuera de los modelos ({ts[0]}-{ts[-1]} K)")
    i = max(j for j in range(len(ts) - 1) if ts[j] <= T) if T < ts[-1] else len(ts) - 2
    t0, t1 = ts[i], ts[i + 1]
    f0, f1 = leer_fits_1d(modelos[t0]), leer_fits_1d(modelos[t1])
    w = (T - t0) / (t1 - t0)
    flujo = np.exp((1 - w) * np.log(np.maximum(f0, 1e-300)) + w * np.log(np.maximum(f1, 1e-300)))
    return lam, flujo


def fracciones_bandas(lam_um, flujo_lam, bandas_cm1, T_cola):
    """Fraccion de la luz de la estrella en cada banda (cm^-1). El flujo por longitud de onda se integra en
    longitud de onda (trapecios); lo que queda mas alla del ultimo punto del modelo (5,5 µm en PHOENIX) se completa
    con la cola de cuerpo negro a T_cola, escalada para coincidir con el modelo en su ultimo tramo. Devuelve
    (fracciones por banda, fraccion total cubierta por las bandas)."""
    lam = np.asarray(lam_um, dtype=float)
    F = np.asarray(flujo_lam, dtype=float)
    orden = np.argsort(lam)
    lam, F = lam[orden], F[orden]
    acumulado = np.concatenate([[0.0], np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(lam))])
    # cola: Planck en lambda escalado al flujo medio del ultimo 1 % del modelo
    def planck(l):
        return l ** -5 / np.expm1(C2_UM_K / (l * T_cola))
    sel = lam > lam[-1] * 0.99
    escala = np.mean(F[sel]) / np.mean(planck(lam[sel]))
    l_cola = np.geomspace(lam[-1], lam[-1] * 1e3, 20001)
    p = planck(l_cola) * escala
    cola = float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(l_cola)))
    total = acumulado[-1] + cola

    def integral_hasta(l):          # ∫ desde el principio hasta l (µm), con la cola si hace falta
        if l <= lam[-1]:
            return float(np.interp(l, lam, acumulado))
        lc = l_cola[l_cola <= l]
        pc = planck(lc) * escala
        return acumulado[-1] + float(np.sum(0.5 * (pc[1:] + pc[:-1]) * np.diff(lc)))
    fr = []
    for a, b in bandas_cm1:
        l_a, l_b = 1e4 / b, 1e4 / a
        fr.append((integral_hasta(l_b) - integral_hasta(l_a)) / total)
    fr = np.array(fr)
    return fr, float(fr.sum())


# ---------------------------------------------------------------------------------------------------- el ozono
def leer_serdyuchenko(ruta):
    """(longitud de onda en nm (vacio), temperaturas en K (293 ... 193), secciones (temperaturas, longitudes) en
    cm^2/molecula). Las lineas que no empiezan por un numero (cabeceras) se saltan."""
    filas = []
    with open(ruta, encoding="latin-1") as f:
        for linea in f:
            p = linea.split()
            if not p:
                continue
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            if len(v) == 12:
                filas.append(v)
    a = np.array(filas)
    if a.ndim != 2 or len(a) < 10:
        raise ValueError(f"{ruta}: no tiene las 12 columnas esperadas")
    return a[:, 0], np.arange(293.0, 192.0, -10.0), a[:, 1:].T


def seccion_ozono(lam_nm, T, datos):
    """Seccion eficaz del ozono (cm^2) a la temperatura T, interpolando linealmente en T (fuera, la del extremo) y en
    longitud de onda."""
    l, Ts, s = datos
    T = float(np.clip(T, Ts.min(), Ts.max()))
    i = int(np.clip(np.searchsorted(-Ts, -T) - 1, 0, len(Ts) - 2))       # Ts decreciente
    w = (Ts[i] - T) / (Ts[i] - Ts[i + 1])
    fila = (1 - w) * s[i] + w * s[i + 1]
    return np.interp(lam_nm, l, fila, left=np.nan, right=np.nan)


# ---------------------------------------------------------------------------------------------------- colisiones
def leer_cia(ruta):
    """Bloques de un archivo CIA de HITRAN: lista de dicts {par, T, nu (cm^-1), k (cm^5/molecula^2)}. Cabecera de
    cada bloque: simbolo (20 caracteres), numero de onda minimo y maximo, numero de puntos, temperatura, ..."""
    bloques = []
    with open(ruta, encoding="latin-1") as f:
        lineas = f.read().splitlines()
    i = 0
    while i < len(lineas):
        cab = lineas[i]
        if not cab.strip():
            i += 1
            continue
        par = cab[:20].strip()
        campos = cab[20:].split()
        n = int(campos[2])
        T = float(campos[3])
        datos = np.array([[float(x) for x in lineas[i + 1 + j].split()[:2]] for j in range(n)])
        bloques.append({"par": par, "T": T, "nu": datos[:, 0], "k": datos[:, 1]})
        i += 1 + n
    return bloques


def k_cia(bloques, nu, T):
    """Coeficiente CIA (cm^5/molecula^2) en nu a la temperatura T: interpolacion lineal en T entre los bloques que
    la rodean y que cubren nu (los bloques pueden cubrir rangos distintos); 0 donde no hay datos."""
    nu = np.asarray(nu, dtype=float)
    out = np.zeros_like(nu)
    por_rango = {}
    for b in bloques:
        por_rango.setdefault((round(b["nu"][0], 3), round(b["nu"][-1], 3)), []).append(b)
    rangos = sorted(por_rango)
    for (a0, a1), (b0, b1) in zip(rangos[:-1], rangos[1:]):
        if b0 < a1:              # dos medidas del mismo tramo: no se suman a ciegas (habra que elegir una)
            raise ValueError(f"bloques CIA solapados: {a0}-{a1} y {b0}-{b1} cm^-1")
    for grupo in por_rango.values():
        grupo = sorted(grupo, key=lambda b: b["T"])
        Ts = np.array([b["T"] for b in grupo])
        dentro = (nu >= grupo[0]["nu"][0]) & (nu <= grupo[0]["nu"][-1])
        if not dentro.any():
            continue
        Tc = float(np.clip(T, Ts[0], Ts[-1]))
        j = int(np.clip(np.searchsorted(Ts, Tc) - 1, 0, max(len(Ts) - 2, 0)))
        if len(grupo) == 1:
            v = np.interp(nu[dentro], grupo[0]["nu"], grupo[0]["k"])
        else:
            w = (Tc - Ts[j]) / (Ts[j + 1] - Ts[j])
            v = ((1 - w) * np.interp(nu[dentro], grupo[j]["nu"], grupo[j]["k"])
                 + w * np.interp(nu[dentro], grupo[j + 1]["nu"], grupo[j + 1]["k"]))
        out[dentro] += np.maximum(v, 0.0)
    return out

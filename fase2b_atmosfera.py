# fase2b_atmosfera.py -- Fase 2b: atmosfera con cuerpo.
#
# Diseño completo, fuentes y decisiones: DISENO_FASE2B.md. Resumen:
#   - Dos capas de aire con temperatura y capacidad propias por celda:
#     capa limite (CL, 1000-900 hPa) y troposfera (TR, 900-0 hPa).
#   - Reparto de la luz segun el balance terrestre sin nubes (Wild et
#     al. 2019): tau = 0.18; de lo atenuado, 73 % lo absorbe la
#     atmosfera y 27 % vuelve al espacio.
#   - Infrarrojo de dos capas grises (Pierrehumbert 2010, cap. 4).
#   - Calor sensible suelo <-> CL, formula bulk con estabilidad de
#     Louis (1979), resuelto de forma implicita.
#   - Ajuste convectivo CL <-> TR (Manabe y Wetherald 1967), 6.5 K/km.
#   - Difusion horizontal reubicada: parte atmosferica en la TR, parte
#     oceanica en la capa de mezcla del oceano (Trenberth y Caron 2001).
#   - Altitud dentro de la fisica (aire de referencia mas frio con la
#     altura) y temperatura del aire a 2 m como temperatura de referencia.
#
# SIETE INTERRUPTORES (INTERRUPTORES_FASE2B). Con todos apagados, el
# modelo reproduce la Fase 2 (fase2_combinado.py, v2.2) -- es la
# prueba V0 del diseño. Restriccion: sin capacidad de la atmosfera
# (I2 apagado) la atmosfera es "instantanea" (modelo de una capa
# clasico) y solo admite I1; I3, I4 e I5 necesitan I2.

import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from parametros import CONSTANTE_SB, ROTACION_PERIODO, P3N_GRAVEDAD
from temperatura import PASO_TIEMPO as PASO_TIEMPO_POR_DEFECTO
from rejilla import FILAS, COLUMNAS, LATITUDES_GRADOS, LONGITUDES_GRADOS
from fase1_geografia import TIERRA, AGUA, correccion_altitud, GRADIENTE_TERMICO
from fase2_difusion import construir_matriz_difusion, DPHI_RAD, DLAMBDA_RAD
from fase2_inercia_multicapa import (
    construir_columna_rejilla, preparar_conduccion_implicita,
    K_DIFUSIVIDAD_TIERRA, N_CAPAS_DEFECTO,
)

# ================================================================
# INTERRUPTORES
# ================================================================
INTERRUPTORES_FASE2B = {
    "luz_absorbida": True,          # I1: reparto de la luz (Wild 2019) en vez de tau=0.3 perdido
    "capacidad_atmosfera": True,    # I2: CL y TR con temperatura y capacidad propias
    "calor_sensible": True,         # I3: calor sensible suelo <-> CL (Louis 1979)
    "capa_limite_radiativa": True,  # I4: la CL absorbe/emite infrarrojo (eps_b > 0)
    "ajuste_convectivo": True,      # I5: ajuste convectivo CL <-> TR
    "difusion_reubicada": True,     # I6: difusion en TR + oceano, no en la piel del suelo
    "altitud_en_fisica": True,      # I7: altitud dentro de la fisica, no a posteriori
}
INTERRUPTORES_APAGADOS = {k: False for k in INTERRUPTORES_FASE2B}

# ================================================================
# PARAMETROS (ver DISENO_FASE2B.md, seccion 3)
# ================================================================
R_AIRE = 287.05          # J/kg/K
CP_AIRE = 1004.0         # J/kg/K
PRESION_SUPERFICIE = 1.0e5   # Pa (1 bar, decision de Carlos 01/10/2026)
PRESION_TECHO_CL = 0.9e5     # Pa: la CL va de 1000 a 900 hPa

CAPACIDAD_CL = (PRESION_SUPERFICIE - PRESION_TECHO_CL) / P3N_GRAVEDAD * CP_AIRE   # ~1.11e6 J/m2/K
CAPACIDAD_TR = PRESION_TECHO_CL / P3N_GRAVEDAD * CP_AIRE                          # ~9.98e6 J/m2/K

ALTURA_ESCALA = R_AIRE * 255.0 / P3N_GRAVEDAD                  # ~8.1 km
Z_CENTRO_CL = ALTURA_ESCALA * math.log(1000 / 950)             # ~0.41 km
Z_CENTRO_TR = ALTURA_ESCALA * math.log(1000 / 550)             # ~4.83 km
GRADIENTE_ADIABATICO = P3N_GRAVEDAD / CP_AIRE                   # K/m, ~9.0 K/km
GRADIENTE_CRITICO = 6.5e-3                                       # K/m (Manabe y Wetherald 1967)
DT_CRITICA_CL_TR = GRADIENTE_CRITICO * (Z_CENTRO_TR - Z_CENTRO_CL)   # ~28.7 K

# Luz (Wild et al. 2019, sin nubes)
TAU_NUEVO = 0.18                     # transmision media 72.6 % con la masa de aire de M3N
FRACCION_ATENUADO_ABSORBIDA = 0.73   # del atenuado: 73 % atmosfera, 27 % espacio
FACTOR_DIFUSIVIDAD = 1.66            # luz reflejada por el suelo, difusa al subir

# Infrarrojo: CALIBRADO (02/10/2026) para reproducir el balance terrestre
# sin nubes de Wild et al. (2019) -- infrarrojo hacia el suelo 314 W/m2,
# hacia el espacio 267 W/m2 -- con las temperaturas de la Atmosfera
# Estandar (US 1976) promediadas en masa en cada capa: suelo 289 K
# (emision 398 W/m2, Wild), CL 285.3 K, TR 246.3 K. Absorcion total del
# infrarrojo: 1-(1-0.740)(1-0.661) = 0.912 (antes 0.77 en una sola capa,
# que estaba calibrada para la Tierra CON nubes).
EMISIVIDAD_CL = 0.740
EMISIVIDAD_TR = 0.661

# Calor sensible
VIENTO = 5.0                     # m/s, provisional hasta la Fase 6
KARMAN = 0.4
Z_REF = 10.0                     # m, altura de referencia del aire
Z0M = {TIERRA: 0.01, AGUA: 2e-4}       # suelo desnudo (Carlos 02/10); oceano ~Charnock a 5 m/s
Z0H = {TIERRA: 0.001, AGUA: 2e-5}      # z0h = z0m/10
C_HN_OCEANO = 1.1e-3             # Large y Pond 1982; Smith 1988
LOUIS_B, LOUIS_C, LOUIS_D = 5.0, 5.0, 5.0

DEPURAR = False

# Difusion reubicada (I6): fraccion oceanica de D. CALIBRADA (02/10/2026)
# para que el oceano lleve ~15 % del transporte hacia los polos a 35
# grados (Trenberth y Caron 2001: 22 % en el norte, 8 % en el sur),
# diagnosticando el transporte de cada componente en el mapa prueba1.
FRACCION_DIFUSION_OCEANO = 0.13


def coeficientes_neutros(tipo_superficie):
    """C_N (momento) y C_HN (calor) neutros por celda."""
    z0m = np.where(tipo_superficie == TIERRA, Z0M[TIERRA], Z0M[AGUA])
    z0h = np.where(tipo_superficie == TIERRA, Z0H[TIERRA], Z0H[AGUA])
    c_n = (KARMAN / np.log(Z_REF / z0m)) ** 2
    c_hn_tierra = KARMAN ** 2 / (np.log(Z_REF / z0m) * np.log(Z_REF / z0h))
    c_hn = np.where(tipo_superficie == TIERRA, c_hn_tierra, C_HN_OCEANO)
    return c_n, c_hn, z0m, z0h


def factor_estabilidad_louis(ri, c_n, z0m):
    """Funcion de estabilidad para el calor de Louis (1979), b=c=d=5."""
    inestable = ri < 0
    ri_abs = np.abs(ri)
    f_inestable = 1 - 3 * LOUIS_B * ri / (1 + 3 * LOUIS_B * LOUIS_C * c_n * np.sqrt(Z_REF / z0m * ri_abs))
    f_estable = 1 / (1 + 3 * LOUIS_B * ri * np.sqrt(1 + LOUIS_D * ri_abs))
    return np.where(inestable, f_inestable, f_estable)


def construir_matriz_difusion_enmascarada(D, mascara):
    """
    Como construir_matriz_difusion() (fase2_difusion.py), pero el flujo
    entre dos celdas vecinas solo existe si AMBAS estan en la mascara
    (p. ej. solo oceano: no hay transporte oceanico a traves de la
    costa). Con D uniforme el operador conserva la energia igual que el
    original (flujo simetrico ponderado por area). Sin la mascara (todo
    True) coincide con construir_matriz_difusion(D*unos).
    """
    phi = np.radians(LATITUDES_GRADOS)
    cos_phi = np.cos(phi)
    cos_n = np.cos(phi + DPHI_RAD / 2)
    cos_s = np.cos(phi - DPHI_RAD / 2)
    filas_idx, cols_idx, vals = [], [], []

    def idx(i, j):
        return i * COLUMNAS + j

    for i in range(FILAS):
        coef_zonal = D / (DLAMBDA_RAD ** 2 * cos_phi[i] ** 2)
        coef_norte = D * cos_n[i] / (DPHI_RAD ** 2 * cos_phi[i])
        coef_sur = D * cos_s[i] / (DPHI_RAD ** 2 * cos_phi[i])
        for j in range(COLUMNAS):
            if not mascara[i, j]:
                continue
            p = idx(i, j)
            diag = 0.0
            vecinos = [(i, (j + 1) % COLUMNAS, coef_zonal), (i, (j - 1) % COLUMNAS, coef_zonal)]
            if i > 0:
                vecinos.append((i - 1, j, coef_norte))
            if i < FILAS - 1:
                vecinos.append((i + 1, j, coef_sur))
            for i2, j2, coef in vecinos:
                if mascara[i2, j2]:
                    filas_idx.append(p); cols_idx.append(idx(i2, j2)); vals.append(coef)
                    diag -= coef
            filas_idx.append(p); cols_idx.append(p); vals.append(diag)
    N = FILAS * COLUMNAS
    return sp.csr_matrix((vals, (filas_idx, cols_idx)), shape=(N, N))


def preparar_luz(albedo_grid, tau, reparto_nuevo):
    """
    Devuelve luz(toa, decl, ang_h_lon0) -> (absorbida_suelo, absorbida_atmosfera).
    Misma geometria y masa de aire que preparar_irradiancia_absorbida_
    rejilla() (fase1_geografia.py). reparto_nuevo=False: esquema antiguo
    (lo atenuado se pierde; atmosfera = 0).
    """
    lat_rad = np.radians(LATITUDES_GRADOS).reshape(-1, 1)
    sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
    lon_rad = np.radians(LONGITUDES_GRADOS)
    uno_menos_albedo = 1 - albedo_grid
    t_subida = math.exp(-tau * FACTOR_DIFUSIVIDAD)

    def luz(toa, declinacion, ang_h_lon0):
        ang_h = (ang_h_lon0 + lon_rad).reshape(1, -1)
        cos_cenital = sin_lat * np.sin(declinacion) + cos_lat * np.cos(declinacion) * np.cos(ang_h)
        cos_cenital = np.clip(cos_cenital, -1.0, 1.0)
        cenital = np.arccos(cos_cenital)
        altura_solar = 90 - np.degrees(cenital)
        dia = altura_solar > 0
        suelo = np.zeros(cenital.shape)
        atmos = np.zeros(cenital.shape)
        if dia.any():
            h = altura_solar[dia]
            masa = 1 / np.sin(np.radians(h + 244 / (165 + 47 * h ** 1.1)))
            transmitancia = np.exp(-masa * tau)
            inst = toa * np.cos(cenital[dia])
            llega = inst * transmitancia
            suelo[dia] = llega * uno_menos_albedo[dia]
            if reparto_nuevo:
                reflejada = llega * albedo_grid[dia]
                atmos[dia] = FRACCION_ATENUADO_ABSORBIDA * (inst * (1 - transmitancia) + reflejada * (1 - t_subida))
        return suelo, atmos

    return luz


def simular_fase2b(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
    interruptores=None, K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.015,
    eps_cl=EMISIVIDAD_CL, eps_tr=EMISIVIDAD_TR, fraccion_oceano=FRACCION_DIFUSION_OCEANO,
    acelerar=True,
):
    """
    Simulacion de la Fase 2b. Devuelve un dict con:
      'anos', 'T_final' (C), 'reg_min', 'reg_media', 'reg_max' (temperatura
      de referencia: aire a 2 m si hay atmosfera con cuerpo; si no, suelo
      como en la Fase 2), 'suelo_min/media/max', 'cl_media', 'tr_media'
      (registros diarios, C), 'energia' (balance en lo alto de la
      atmosfera, año final).
    """
    I = dict(INTERRUPTORES_FASE2B if interruptores is None else interruptores)
    atm = I["capacidad_atmosfera"]
    if not atm and (I["calor_sensible"] or I["capa_limite_radiativa"] or I["ajuste_convectivo"]):
        raise ValueError("calor_sensible, capa_limite_radiativa y ajuste_convectivo necesitan capacidad_atmosfera")
    if I["difusion_reubicada"] and not atm:
        raise ValueError("difusion_reubicada necesita capacidad_atmosfera (la difusion atmosferica actua sobre la TR)")

    sigma = CONSTANTE_SB
    es_tierra = tipo_superficie == TIERRA
    albedo_grid = np.where(es_tierra, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    tau = TAU_NUEVO if I["luz_absorbida"] else profundidad_optica
    luz = preparar_luz(albedo_grid, tau, I["luz_absorbida"])
    eps_b = eps_cl if I["capa_limite_radiativa"] else 0.0
    eps_t = eps_tr if I["capa_limite_radiativa"] else emisividad

    capacidades_reales, conductancias = construir_columna_rejilla(tipo_superficie, inercia_por_tipo, K_difusividad, n_capas)

    # ---- difusion horizontal ----
    D_grid = np.full((FILAS, COLUMNAS), D)
    L_total = construir_matriz_difusion(D_grid)
    if I["difusion_reubicada"]:
        L_tr = construir_matriz_difusion(np.full((FILAS, COLUMNAS), D * (1 - fraccion_oceano)))
        fact_tr = splu((sp.diags(np.full(FILAS * COLUMNAS, CAPACIDAD_TR / paso_tiempo)) - L_tr).tocsc())
        L_oc = construir_matriz_difusion_enmascarada(D * fraccion_oceano, ~es_tierra)

    def preparar_superficie(capacidades):
        """Todo lo que depende de la capacidad de la superficie (cambia
        durante el arranque rapido, ver mas abajo)."""
        C0 = capacidades[..., 0]
        resolver = preparar_conduccion_implicita(capacidades, conductancias, paso_tiempo)
        L_sup = L_oc if I["difusion_reubicada"] else L_total
        fact = splu((sp.diags(C0.flatten() / paso_tiempo) - L_sup).tocsc())
        return C0, resolver, fact

    # ---- calor sensible ----
    c_n, c_hn, z0m, z0h = coeficientes_neutros(tipo_superficie)
    desfase_theta = GRADIENTE_ADIABATICO * Z_CENTRO_CL - (
        GRADIENTE_TERMICO * altitud_metros / 1000 if I["altitud_en_fisica"] else 0.0)
    fraccion_2m = np.log(2.0 / z0h) / np.log(Z_REF / z0h)

    # ---- estado inicial: balance medio anual con difusion (Newton) ----
    suma_s = np.zeros((FILAS, COLUMNAS)); suma_a = np.zeros((FILAS, COLUMNAS))
    for toa, decl, ang in datos_orbita:
        s, a = luz(toa, decl, ang)
        suma_s += s; suma_a += a
    n_pasos = len(datos_orbita)
    eps_eff = 1 - (1 - eps_b) * (1 - eps_t)
    forzante = ((suma_s + suma_a / 2) / n_pasos).flatten()
    k = (1 - eps_eff / 2) * sigma
    Ts = np.maximum((forzante / k) ** 0.25, 150.0)
    for _ in range(30):
        res = forzante - k * Ts ** 4 + L_total @ Ts
        corr = splu((sp.diags(-4 * k * Ts ** 3) + L_total).tocsc()).solve(-res)
        Ts = Ts + corr
        if np.max(np.abs(corr)) < 1e-6:
            break
    Ts = Ts.reshape(FILAS, COLUMNAS)
    T_col = np.repeat(Ts[..., None], n_capas, axis=2)
    T_tr = (Ts ** 4 / 2) ** 0.25
    T_cl = Ts - desfase_theta

    def paso(T_col, T_cl, T_tr, s_abs, a_abs, acumular=None):
        C0, resolver_conduccion, fact_sup = SUPERFICIE
        Ts = T_col[..., 0]
        F0 = sigma * Ts ** 4
        if atm:
            Eb = sigma * T_cl ** 4
            Et = sigma * T_tr ** 4
            sube_cl = (1 - eps_b) * F0 + eps_b * Eb
            baja_tr = eps_t * Et
            dlr = (1 - eps_b) * baja_tr + eps_b * Eb
            olr = (1 - eps_t) * sube_cl + eps_t * Et
            neto_s = s_abs + dlr - F0
            neto_b = eps_b * (F0 + baja_tr) - 2 * eps_b * Eb
            neto_t = eps_t * sube_cl - 2 * eps_t * Et + a_abs
            T_cl = T_cl + (paso_tiempo / CAPACIDAD_CL) * neto_b
            T_tr = T_tr + (paso_tiempo / CAPACIDAD_TR) * neto_t
        else:
            # atmosfera instantanea (una capa, sin capacidad): eps*F0 + a = 2*eps*sigma*Tt^4
            dlr = (eps_t * F0 + a_abs) / 2
            olr = (1 - eps_t) * F0 + dlr
            if I["luz_absorbida"]:
                neto_s = s_abs + dlr - F0
            else:
                neto_s = s_abs - (1 - emisividad / 2) * CONSTANTE_SB * Ts ** 4   # forma exacta de la Fase 2
        if acumular is not None:
            acumular["abs"] += np.sum((s_abs + a_abs) * PESO)
            acumular["olr"] += np.sum(olr * PESO)
        T_col = T_col.copy()
        T_col[..., 0] = Ts + (paso_tiempo / C0) * neto_s

        # difusion horizontal (implicita)
        if I["difusion_reubicada"]:
            T_tr = fact_tr.solve((CAPACIDAD_TR / paso_tiempo) * T_tr.flatten()).reshape(FILAS, COLUMNAS)
            T_col[..., 0] = fact_sup.solve((C0.flatten() / paso_tiempo) * T_col[..., 0].flatten()).reshape(FILAS, COLUMNAS)
        else:
            T_col[..., 0] = fact_sup.solve((C0.flatten() / paso_tiempo) * T_col[..., 0].flatten()).reshape(FILAS, COLUMNAS)

        # calor sensible suelo <-> CL (implicito, cerrado, conserva energia)
        if I["calor_sensible"]:
            Ts = T_col[..., 0]
            theta_aire = T_cl + desfase_theta
            dif = Ts - theta_aire
            theta_media = 0.5 * (Ts + theta_aire)
            ri = P3N_GRAVEDAD * Z_REF * (theta_aire - Ts) / (theta_media * VIENTO ** 2)
            c_h = c_hn * factor_estabilidad_louis(ri, c_n, z0m)
            rho = PRESION_SUPERFICIE / (R_AIRE * theta_aire)
            g = rho * CP_AIRE * c_h * VIENTO
            a_s = C0 / paso_tiempo
            a_b = CAPACIDAD_CL / paso_tiempo
            dif_nueva = dif / (1 + g * (1 / a_s + 1 / a_b))
            flujo = g * dif_nueva
            T_col[..., 0] = Ts - flujo / a_s
            T_cl = T_cl + flujo / a_b

        # ajuste convectivo CL <-> TR (instantaneo, conserva energia)
        if I["ajuste_convectivo"]:
            exceso = (T_cl - T_tr) > DT_CRITICA_CL_TR
            if exceso.any():
                energia = CAPACIDAD_CL * T_cl + CAPACIDAD_TR * T_tr
                tr_aj = (energia - CAPACIDAD_CL * DT_CRITICA_CL_TR) / (CAPACIDAD_CL + CAPACIDAD_TR)
                T_tr = np.where(exceso, tr_aj, T_tr)
                T_cl = np.where(exceso, tr_aj + DT_CRITICA_CL_TR, T_cl)

        T_col = resolver_conduccion(T_col)
        return T_col, T_cl, T_tr

    PESO = np.cos(np.radians(LATITUDES_GRADOS)).reshape(-1, 1) * np.ones((1, COLUMNAS))

    SUPERFICIE = preparar_superficie(capacidades_reales)

    # ACELERACION DE LA CONVERGENCIA (solo con atmosfera con cuerpo).
    # Tras los primeros años, lo que queda por converger es casi siempre
    # un unico "modo lento" (oceano polar, que de noche casi no
    # intercambia calor con el aire estable): cada año cambia un
    # porcentaje casi constante r del año anterior. Cuando r se mantiene
    # estable tres años seguidos (variacion < 0.03), se salta de una vez
    # lo que falta, delta*r/(1-r) -- suma de la serie geometrica -- en
    # todo el estado (suelo, CL, TR). Despues se siguen simulando años
    # normales y la convergencia se decide con el MISMO criterio de
    # siempre: el salto solo acorta el camino. Validado contra la
    # simulacion sin aceleracion.
    historial_r = []
    cambio_anterior = None
    delta_anterior = None
    saltos = 0

    anos = max_anos
    for ano in range(max_anos):
        Ts0, Tt0 = T_col[..., 0].copy(), T_tr.copy()
        estado0 = (T_col.copy(), T_cl.copy(), T_tr.copy())
        for toa, decl, ang in datos_orbita:
            s, a = luz(toa, decl, ang)
            T_col, T_cl, T_tr = paso(T_col, T_cl, T_tr, s, a)
        cambio = max(np.max(np.abs(T_col[..., 0] - Ts0)), np.max(np.abs(T_tr - Tt0)) if atm else 0.0)
        if atm and acelerar and cambio >= tolerancia_convergencia:
            if cambio_anterior is not None:
                historial_r.append(cambio / cambio_anterior)
            cambio_anterior = cambio
            if len(historial_r) >= 3 and max(historial_r[-3:]) - min(historial_r[-3:]) < 0.03 and historial_r[-1] < 0.9:
                r = historial_r[-1]
                factor = r / (1 - r)
                T_col = T_col + factor * (T_col - estado0[0])
                T_cl = T_cl + factor * (T_cl - estado0[1])
                T_tr = T_tr + factor * (T_tr - estado0[2])
                saltos += 1
                historial_r = []
                cambio_anterior = None
                if DEPURAR:
                    print(f"  -> salto con r={r:.3f} (factor {factor:.2f})", flush=True)
        if DEPURAR:
            ds = np.abs(T_col[..., 0] - Ts0); dt_ = np.abs(T_tr - Tt0)
            i = np.unravel_index(np.argmax(ds), ds.shape); j = np.unravel_index(np.argmax(dt_), dt_.shape)
            print(f"  año {ano+1}: suelo {ds.max():.4f} en {i} ({'tierra' if es_tierra[i] else 'agua'}), TR {dt_.max():.4f} en {j}", flush=True)
        if cambio < tolerancia_convergencia:
            anos = ano + 1
            break

    # ---- año final con registro ----
    pasos_dia = round(ROTACION_PERIODO / paso_tiempo)
    acum = {"abs": 0.0, "olr": 0.0}
    reg = {k: [] for k in ("min", "media", "max", "s_min", "s_media", "s_max", "cl", "tr")}

    def temp_referencia(T_col, T_cl):
        Ts = T_col[..., 0]
        if I["calor_sensible"]:
            return Ts + fraccion_2m * (T_cl + desfase_theta - Ts)   # aire a 2 m
        return Ts

    def a_celsius(T):
        T = T - 273.15
        return T if I["altitud_en_fisica"] else correccion_altitud(T, altitud_metros)

    for p, (toa, decl, ang) in enumerate(datos_orbita):
        if p % pasos_dia == 0:
            if p > 0:
                for clave, valor in (("min", mn), ("max", mx), ("s_min", smn), ("s_max", smx)):
                    reg[clave].append(a_celsius(valor))
                reg["media"].append(a_celsius(suma / cnt)); reg["s_media"].append(a_celsius(ssuma / cnt))
                reg["cl"].append(cl_suma / cnt - 273.15); reg["tr"].append(tr_suma / cnt - 273.15)
            mn = np.full((FILAS, COLUMNAS), np.inf); mx = -mn; suma = np.zeros((FILAS, COLUMNAS))
            smn = mn.copy(); smx = mx.copy(); ssuma = suma.copy(); cl_suma = suma.copy(); tr_suma = suma.copy(); cnt = 0
        s, a = luz(toa, decl, ang)
        T_col, T_cl, T_tr = paso(T_col, T_cl, T_tr, s, a, acum)
        tref = temp_referencia(T_col, T_cl)
        mn = np.minimum(mn, tref); mx = np.maximum(mx, tref); suma = suma + tref
        Ts = T_col[..., 0]
        smn = np.minimum(smn, Ts); smx = np.maximum(smx, Ts); ssuma = ssuma + Ts
        cl_suma = cl_suma + T_cl; tr_suma = tr_suma + T_tr; cnt += 1
    if cnt == pasos_dia:
        for clave, valor in (("min", mn), ("max", mx), ("s_min", smn), ("s_max", smx)):
            reg[clave].append(a_celsius(valor))
        reg["media"].append(a_celsius(suma / cnt)); reg["s_media"].append(a_celsius(ssuma / cnt))
        reg["cl"].append(cl_suma / cnt - 273.15); reg["tr"].append(tr_suma / cnt - 273.15)

    return {
        "anos": anos, "saltos": saltos,
        "T_final": a_celsius(temp_referencia(T_col, T_cl)),
        "reg_min": np.array(reg["min"]), "reg_media": np.array(reg["media"]), "reg_max": np.array(reg["max"]),
        "suelo_min": np.array(reg["s_min"]), "suelo_media": np.array(reg["s_media"]), "suelo_max": np.array(reg["s_max"]),
        "cl_media": np.array(reg["cl"]), "tr_media": np.array(reg["tr"]),
        "energia": {"absorbido": acum["abs"], "olr": acum["olr"],
                    "diferencia_relativa": abs(acum["abs"] - acum["olr"]) / acum["abs"]},
    }

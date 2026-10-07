# fase31_agua.py -- v3.1 (Fase 5a): ciclo del agua sobre la atmosfera de N capas.
#
# Diseño, fuentes y decisiones: DISENO_FASE5A.md (decisiones de Carlos) y
# DISENO_V3.1.md (implementacion). Este modulo contiene la fisica "de
# columna" del agua, vectorizada sobre todas las celdas; la usa
# simular_fase2b() (fase2b_atmosfera.py) con los interruptores I11-I15.
#
# Contenido:
#   - saturacion del vapor (Ambaum 2020, sobre agua y sobre hielo);
#   - conveccion humeda: Betts-Miller simplificado de Frierson (2007),
#     traduccion fiel del codigo de Isca (qe_moist_convection.F90);
#   - condensacion de gran escala (Isca, lscale_cond.F90, sin reevaporacion);
#   - fraccion de lluvia/nieve (rampa del CLM5);
#   - albedo de la nieve y del hielo segun el espectro de la estrella;
#   - propiedades termicas del suelo segun su agua (CLM5: Farouki 1981).
#
# Convenios: capas con indice 0 = la mas alta (como fase30_multicapa); q es
# la humedad ESPECIFICA (kg de vapor por kg de aire humedo); unidades SI.

import math
import numpy as np

try:
    from numba import njit, prange
    HAY_NUMBA = True
except ImportError:
    HAY_NUMBA = False
    prange = range

from fase30_multicapa import R_AIRE, CP_AIRE, KAPPA, P0

# ================================================================
# CONSTANTES DEL AGUA (a) -- DISENO_FASE5A.md, seccion 5
# ================================================================
# Calores latentes CONSTANTES en el balance de energia (decision 1 de Carlos,
# opcion A). L_F es el MISMO valor que ya usa el hielo marino (fase2b,
# L_FUSION = 3,34e5; Isca HLF = 3,34e5), para que la nieve que cae sobre el
# hielo y el hielo cierren la energia entre si.
L_V = 2.50084e6          # J/kg (MetPy, a T0)
L_F = 3.34e5             # J/kg
L_S = L_V + L_F          # J/kg
R_V = 461.52             # J/kg/K
EPSILON = R_AIRE / R_V   # 0,62197 (masa molar del agua / del aire seco)
# Solo para la SATURACION (formula de Ambaum 2020, ecs. 13 y 17):
T0_SAT = 273.16          # K, punto triple
E0_SAT = 611.2           # Pa (MetPy)
C_PL = 4219.4            # J/kg/K, agua liquida (MetPy)
C_PV = 1860.0            # J/kg/K, vapor (MetPy)
C_PI = 2090.0            # J/kg/K, hielo (MetPy)
T_FUSION = 273.15        # K, agua dulce (nieve, suelo)

# ================================================================
# PARAMETROS DE LA FASE 5a (c salvo que se diga otra cosa)
# ================================================================
# Conveccion humeda (Frierson 2007; Isca qe_moist_convection, valores por defecto ✅):
TAU_CONVECCION = 7200.0  # s
RH_REF_CONVECCION = 0.8  # (c) INCIERTO: 0,6 (control de Frierson 2007b) - 0,8 (Isca) -> prueba de sensibilidad
T_MIN_LCL = 173.0        # K (Isca)
SIGMA_TOPE_CONVECCION = 0.1   # la conveccion no sube por encima de p = 0,1 p_superficie (~16 km en la Tierra)
# Condensacion de gran escala: humedad relativa 1 (Isca, hc = 1 ✅; Manabe et al. 1965 ✅). Con N = 20
# capas la resolucion vertical es la de un modelo de capas finas (como Isca), no la de las dos capas
# gruesas para las que se habia previsto una RH critica < 1 (DISENO_FASE5A 3.7). Sensibilidad: 0,9.
RH_CONDENSACION = 1.0
# Lluvia o nieve: rampa lineal con el aire a 2 m, todo nieve por debajo de 0 C y todo lluvia por
# encima de 2 C (CLM5, nota tecnica ✅). Jennings et al. (2018) ✅: umbral del 50 % observado en
# tierra de 1,0 C de media, coherente con el centro de la rampa.
T_TODO_NIEVE = 273.15
T_TODO_LLUVIA = 275.15
# Cubo de Manabe (1969): capacidad 15 cm y beta = min(1, W / (0,75 W_fc)) (Vallis et al. 2018,
# Isca ✅, secundaria)
W_CAMPO = 150.0          # kg/m2 (= mm)
UMBRAL_BETA = 0.75
# Nieve en tierra:
S_MAX = 1000.0           # kg/m2: tope (CLM4.5 usaba 1000 mm ✅ via DISENO_FASE5A 3.10); lo que sobra
                         # vuelve al oceano como descarga glaciar (agua y energia conservadas)
S_CUBRE = 10.0           # kg/m2 (c) ⚠️: con S = S_CUBRE la nieve cubre la mitad de la celda a efectos
                         # de albedo (fraccion S/(S + S_CUBRE)); unos 4 cm de nieve de 250 kg/m3


# Vapor radiativo (I15, PROTOTIPO, 05/10/2026): d(tau)/d(p/P0) = A + B*q (Byrne y O'Gorman 2013;
# forma verificada en Isca two_stream_gray_rad, esquema 'byrne', con A = 0,8678 y B = 1997,9 para
# SU modelo). Recalibrado para M3N igual que la v3.0 (calibrar_tau_vapor.py): con la temperatura de
# la Atmosfera Estandar 1976 (suelo 289 K), la forma del perfil de vapor de la atmosfera estandar de
# EE. UU. de la AFGL (Anderson et al. 1986) escalada al agua precipitable media de la Tierra
# (24,9 kg/m2, Trenberth y Smith 2005), reproduce los flujos sin nubes de Wild et al. (2019):
# DLR 314 y OLR 267 W/m2. Espesor total de la referencia 2,45 (el mismo que el tau fijo de la v3.0).
A_LW_SECO = 0.5656
B_LW_VAPOR = 773.1

# ================================================================
# SATURACION (a) -- Ambaum (2020); verificado contra MetPy
# ================================================================
def es_agua(T):
    """Presion de vapor de saturacion sobre agua liquida (Pa)."""
    L = L_V - (C_PL - C_PV) * (T - T0_SAT)
    return E0_SAT * (T0_SAT / T) ** ((C_PL - C_PV) / R_V) * np.exp((L_V / T0_SAT - L / T) / R_V)


def es_hielo(T):
    """Presion de vapor de saturacion sobre hielo (Pa)."""
    L = L_S - (C_PI - C_PV) * (T - T0_SAT)
    return E0_SAT * (T0_SAT / T) ** ((C_PI - C_PV) / R_V) * np.exp((L_S / T0_SAT - L / T) / R_V)


def _es_tabla(T, tabla):
    """e_s por interpolacion lineal en una tabla cada 0,005 K (error relativo < 3e-7)."""
    x = np.clip((T - T_TABLA_MIN) / T_TABLA_PASO, 0.0, N_TABLA - 1.000001)
    i = x.astype(np.int64)
    f = x - i
    return tabla[i] * (1.0 - f) + tabla[i + 1] * f


def qs_y_derivada(T, p, hielo=None):
    """Humedad especifica de saturacion q_s = eps*e/(p - (1-eps)*e) y su derivada respecto a T.
    hielo: mascara (mismo forma o difundible) donde se usa la saturacion sobre hielo.
    d(ln e_s)/dT = L(T)/(R_v T^2) exactamente (Clausius-Clapeyron con el L de Kirchhoff)."""
    if hielo is True:                                  # todo sobre hielo (sin calcular lo del agua)
        e = _es_tabla(T, _TABLA_ES_HIELO)
        Lt = L_S - (C_PI - C_PV) * (T - T0_SAT)
    else:
        e = _es_tabla(T, _TABLA_ES)
        Lt = L_V - (C_PL - C_PV) * (T - T0_SAT)
    if hielo is not None and hielo is not True:
        e = np.where(hielo, _es_tabla(T, _TABLA_ES_HIELO), e)
        Lt = np.where(hielo, L_S - (C_PI - C_PV) * (T - T0_SAT), Lt)
    e = np.minimum(e, 0.5 * p)                       # salvaguarda (solo actua por encima de ~80 C)
    den = p - (1 - EPSILON) * e
    qs = EPSILON * e / den
    dqs = qs * p / den * Lt / (R_V * T * T)
    return qs, dqs


# ================================================================
# CONVECCION HUMEDA: BETTS-MILLER SIMPLIFICADO (Frierson 2007)
# ================================================================
# Traduccion de qe_moist_convection.F90 de Isca (SBM_convection_scheme y sus
# subrutinas), columna a columna. Diferencias, todas de detalle:
#   - la temperatura del LCL se resuelve por Newton directamente (Isca usa
#     una tabla precalculada con el mismo Newton);
#   - la saturacion es la de Ambaum sobre agua (Isca: su tabla escomp);
#   - las constantes son las de M3N (R, c_p, L_v) y la g de cada planeta.
# Conserva exactamente la entalpia humeda de la columna (cp*T + L_v*q) y el
# agua (lo que se condensa es la precipitacion).

def _sbm_columnas(T, q, pf, ph, lpf, pk, dlnph, dt, tau, rhbm, g, cp, rd, rv, lv, tmin, sigma_tope):
    n, m = T.shape
    dT = np.zeros((n, m)); dq = np.zeros((n, m)); prec = np.zeros(m)
    kappa = rd / cp
    eps = rd / rv
    for c in prange(m):
        Tp = np.empty(n); rp = np.empty(n); rin = np.empty(n); Tv = np.empty(n)
        dTp = np.empty(n); dqp = np.empty(n); Tref = np.empty(n); qref = np.empty(n)
        ks = n - 1                                     # capa de superficie
        # tope de la conveccion: la capa mas alta con p >= sigma_tope * p_superficie (en Isca, con muchas
        # capas, la flotacion se acaba siempre por debajo de la tropopausa; aqui se acota por seguridad)
        ktope = 0
        while ktope < ks and pf[ktope, c] < sigma_tope * ph[n, c]:
            ktope += 1
        for k in range(n):
            rin[k] = q[k, c] / (1.0 - q[k, c])
            Tv[k] = T[k, c] * (1.0 + q[k, c] * (rv / rd - 1.0))
            Tp[k] = T[k, c]; rp[k] = rin[k]
        # ---- CAPE (CAPE_calculation) ----
        nocape = True; cape = 0.0; kLZB = -1; kLFC = -1; skip = False
        T0 = T[ks, c]; r0 = rin[ks]
        es = _es_agua_escalar(T0)
        rs = rd * es / rv / (pf[ks, c] - es)
        if r0 >= rs:                                   # saturado en superficie
            kLCL = ks
            Tp[ks] = T0 + (r0 - rs) / ((cp / lv) + (lv * rs) / rv / (T0 * T0))
            es = _es_agua_escalar(Tp[ks]); rp[ks] = rd * es / rv / (pf[ks, c] - es)
        else:
            theta0 = T[ks, c] * (P0 / pf[ks, c]) ** kappa
            if r0 <= 0.0:
                skip = True; kLCL = 0
            else:
                valor = math.log(theta0 ** (-1.0 / kappa) * P0 * r0 / (eps + r0))
                TLCL = T0
                for _it in range(100):
                    es = _es_agua_escalar(TLCL)
                    f = valor - math.log(es * TLCL ** (-1.0 / kappa))
                    df = 1.0 / kappa / TLCL - _lv_kirchhoff(TLCL) / rv / (TLCL * TLCL)
                    paso = f / df
                    TLCL -= paso
                    if abs(paso) < 1e-7:
                        break
                pLCL = P0 * (TLCL / theta0) ** (1.0 / kappa)
                if pLCL < pf[0, c]:
                    pLCL = pf[0, c]; TLCL = theta0 * (pLCL / P0) ** kappa
                k = ks
                while k >= 0 and pf[k, c] > pLCL:
                    Tp[k] = theta0 * pk[k, c]
                    es = _es_agua_escalar(Tp[k]); rp[k] = rd * es / rv / (pf[k, c] - es)
                    k -= 1
                kLCL = k
                if kLCL < 0:
                    skip = True
                else:
                    a = kappa * TLCL + (lv / cp) * r0
                    b = lv * lv * r0 / (cp * rv * TLCL * TLCL)
                    Tp[kLCL] = TLCL + a / (1.0 + b) * math.log(pf[kLCL, c] / pLCL) / 2.0
                    if Tp[kLCL] < tmin:
                        skip = True
                    else:
                        es = _es_agua_escalar(Tp[kLCL])
                        rp[kLCL] = rd * es / rv / ((pf[kLCL, c] + pLCL) / 2.0 - es)
                        a = kappa * Tp[kLCL] + (lv / cp) * rp[kLCL]
                        b = lv * lv * rp[kLCL] / (cp * rv * Tp[kLCL] * Tp[kLCL])
                        Tp[kLCL] = TLCL + a / (1.0 + b) * math.log(pf[kLCL, c] / pLCL)
                        if Tp[kLCL] < tmin:
                            skip = True
                        else:
                            es = _es_agua_escalar(Tp[kLCL]); rp[kLCL] = rd * es / rv / (pf[kLCL, c] - es)
                            tvp = Tp[kLCL] * (1.0 + rp[kLCL] / (1.0 + rp[kLCL]) * (rv / rd - 1.0))
                            if not (tvp < Tv[kLCL]):
                                cape += rd * (tvp - Tv[kLCL]) * dlnph[kLCL, c]
                                nocape = False; kLFC = kLCL
        if skip:
            nocape = True
        # ---- CAPE por encima del LCL (CAPE_above_LCL) ----
        if not skip:
            for k in range(kLCL - 1, ktope - 1, -1):
                a = kappa * Tp[k + 1] + (lv / cp) * rp[k + 1]
                b = lv * lv * rp[k + 1] / (cp * rv * Tp[k + 1] * Tp[k + 1])
                Tp[k] = Tp[k + 1] + a / (1.0 + b) * (lpf[k, c] - lpf[k + 1, c]) / 2.0
                if Tp[k] < tmin and nocape:
                    break
                es = _es_agua_escalar(Tp[k])
                rp[k] = rd * es / rv / ((pf[k, c] + pf[k + 1, c]) / 2.0 - es)
                a = kappa * Tp[k] + (lv / cp) * rp[k]
                b = lv * lv * rp[k] / (cp * rv * Tp[k] * Tp[k])
                Tp[k] = Tp[k + 1] + a / (1.0 + b) * (lpf[k, c] - lpf[k + 1, c])
                if Tp[k] < tmin and nocape:
                    break
                es = _es_agua_escalar(Tp[k]); rp[k] = rd * es / rv / (pf[k, c] - es)
                tvp = Tp[k] * (1.0 + rp[k] / (1.0 + rp[k]) * (rv / rd - 1.0))
                if tvp < Tv[k]:
                    if not nocape:
                        kLZB = k + 1
                        break
                else:
                    cape += rd * (tvp - Tv[k]) * dlnph[k, c]
                    if nocape:
                        nocape = False; kLFC = k
            if (not nocape) and kLZB < 0:
                kLZB = ktope                            # la flotacion llega al tope de la conveccion
        if nocape or cape <= 0.0 or kLZB < 0:
            continue
        # ---- perfiles de referencia (set_reference_profiles) ----
        for k in range(n):
            Tref[k] = T[k, c]; qref[k] = q[k, c]; dTp[k] = 0.0; dqp[k] = 0.0
        for k in range(kLZB, n):
            Tref[k] = Tp[k]
            eref = rhbm * pf[k, c] * rp[k] / (rp[k] + rd / rv)
            rr = rd * eref / rv / (pf[k, c] - eref)
            qref[k] = rr / (1.0 + rr)
        # ---- Pq, Pt ----
        Pq = 0.0; Pt = 0.0
        for k in range(kLZB, n):
            dqp[k] = -(q[k, c] - qref[k]) * dt / tau
            Pq += dqp[k] * (ph[k, c] - ph[k + 1, c])
            dTp[k] = -(T[k, c] - Tref[k]) * dt / tau
            Pt += (cp / lv) * dTp[k] * (ph[k + 1, c] - ph[k, c])
        Pq /= g; Pt /= g
        if Pq > 0.0 and Pt > 0.0:
            # ---- conveccion profunda ----
            if Pq > Pt:                                # se reduce la tasa de secado (do_change_time_scale_deepconv)
                fac = Pt / Pq
                for k in range(kLZB, n):
                    dqp[k] *= fac
                Pq = Pt
            else:                                      # se corrige la temperatura de referencia (do_change_Tref_deepconv)
                dk = 0.0
                for k in range(kLZB, n):
                    dk -= (dTp[k] + (lv / cp) * dqp[k]) * (ph[k + 1, c] - ph[k, c])
                dk /= (ph[n, c] - ph[kLZB, c])
                for k in range(kLZB, n):
                    dTp[k] += dk
        elif Pt > 0.0:
            # ---- conveccion somera (do_shallow_convection): no llueve ----
            k = kLZB
            while Pq < 0.0 and k <= n - 1:
                Pq -= dqp[k] * (ph[k, c] - ph[k + 1, c]) / g
                k += 1
            ktop = k - 1
            for kk in range(kLZB, ktop):               # por encima del nivel de precipitacion cero: nada
                dTp[kk] = 0.0; dqp[kk] = 0.0
            den = dqp[ktop] * (ph[ktop + 1, c] - ph[ktop, c])
            if Pq > 0.0 and den != 0.0:
                cc = Pq * g / den
                dqp[ktop] *= cc; dTp[ktop] *= cc
                dk = 0.0
                for kk in range(ktop, n):
                    dk += dTp[kk] * (ph[kk, c] - ph[kk + 1, c])
                dk /= (ph[n, c] - ph[ktop, c])
                if ktop != n - 1:
                    for kk in range(ktop, n):
                        dTp[kk] += dk
            else:
                lo = n - 1 if ktop == kLZB else kLZB
                hi = n - 1 if ktop == kLZB else ktop
                for kk in range(lo, hi + 1):
                    dTp[kk] = 0.0; dqp[kk] = 0.0
            Pq = 0.0
        else:
            continue
        for k in range(n):
            dT[k, c] = dTp[k]; dq[k, c] = dqp[k]
        prec[c] = Pq
    return dT, dq, prec


def _dlnp(p_abajo, p_arriba, p_centro):
    """ln(p_abajo/p_arriba) de una capa; en la capa mas alta (p_arriba = 0) se usa el doble del tramo
    entre su centro y su base (Isca obtendria infinito en ese caso limite)."""
    if p_arriba > 0.0:
        return math.log(p_abajo / p_arriba)
    return 2.0 * math.log(p_abajo / p_centro)


def _lv_kirchhoff(T):
    return L_V - (C_PL - C_PV) * (T - T0_SAT)


# Tabla de e_s sobre agua cada 0,005 K entre 100 y 400 K (como la tabla escomp de Isca), con
# interpolacion lineal: error relativo < 3e-7 frente a la formula exacta, bastante mas rapido.
T_TABLA_MIN, T_TABLA_PASO, N_TABLA = 100.0, 0.005, 60001
_TABLA_ES = None


def _construir_tabla_es():
    T = T_TABLA_MIN + T_TABLA_PASO * np.arange(N_TABLA)
    return es_agua(T)


_TABLA_ES = _construir_tabla_es()
_TABLA_ES_HIELO = es_hielo(T_TABLA_MIN + T_TABLA_PASO * np.arange(N_TABLA))


def _es_agua_escalar(T):
    x = (T - T_TABLA_MIN) / T_TABLA_PASO
    if x <= 0.0:
        x = 0.0
    if x >= N_TABLA - 1.000001:
        x = N_TABLA - 1.000001
    i = int(x)
    f = x - i
    return _TABLA_ES[i] * (1.0 - f) + _TABLA_ES[i + 1] * f


if HAY_NUMBA:
    _dlnp = njit(cache=True)(_dlnp)
    _lv_kirchhoff = njit(cache=True)(_lv_kirchhoff)
    _es_agua_escalar = njit(cache=True)(_es_agua_escalar)
    # en paralelo sobre las columnas (cada columna es independiente): usa los nucleos disponibles
    _sbm_columnas = njit(cache=True, parallel=True)(_sbm_columnas)


_CACHE_GEO = {}


def _geometria_conveccion(pf, ph):
    """Logaritmos y potencias de la presion: se calculan una vez y se reutilizan mientras pf y ph no cambien.
    v3.1-pre9: la cache compara el CONTENIDO completo de pf y ph (np.array_equal con una copia guardada).
    Antes la clave era id() de los arrays y tres valores sueltos: con la presion en superficie variable del
    nucleo dinamico (I16) Python puede reutilizar las direcciones de memoria y p_s puede cambiar solo en celdas
    interiores, y la conveccion habria usado en silencio la geometria de un paso anterior."""
    if _CACHE_GEO:
        pf_c, ph_c = _CACHE_GEO["pf"], _CACHE_GEO["ph"]
        vigente = (pf_c.shape == pf.shape and ph_c.shape == ph.shape
                   and np.array_equal(pf_c, pf) and np.array_equal(ph_c, ph))
    else:
        vigente = False
    if not vigente:
        n = pf.shape[0]
        pf2 = np.ascontiguousarray(pf.reshape(n, -1)); ph2 = np.ascontiguousarray(ph.reshape(n + 1, -1))
        lpf = np.log(pf2)
        pk = (pf2 / P0) ** KAPPA
        with np.errstate(divide="ignore"):
            dl = np.where(ph2[:-1] > 0, np.log(ph2[1:] / np.where(ph2[:-1] > 0, ph2[:-1], 1.0)),
                          2.0 * np.log(ph2[1:] / pf2))
        _CACHE_GEO.clear()
        _CACHE_GEO["pf"] = np.array(pf, copy=True); _CACHE_GEO["ph"] = np.array(ph, copy=True)
        _CACHE_GEO["geo"] = (pf2, ph2, np.ascontiguousarray(lpf), np.ascontiguousarray(pk), np.ascontiguousarray(dl))
    return _CACHE_GEO["geo"]


def conveccion_humeda(T, q, pf, ph, dt, gravedad, rh_ref=None, tau=None):
    """Betts-Miller simplificado en todas las columnas. T, q: (n, F, C); pf: (n, F, C);
    ph: (n+1, F, C). Devuelve (T nueva, q nueva, precipitacion kg/m2 en el paso, (F, C))."""
    n = T.shape[0]
    forma = T.shape[1:]
    rh_ref = RH_REF_CONVECCION if rh_ref is None else rh_ref      # se leen al llamar (pruebas de sensibilidad)
    tau = TAU_CONVECCION if tau is None else tau
    geo = _geometria_conveccion(pf, ph)
    dT, dq, prec = _sbm_columnas(np.ascontiguousarray(T.reshape(n, -1)), np.ascontiguousarray(q.reshape(n, -1)),
                                 *geo, float(dt), float(tau), float(rh_ref), float(gravedad), CP_AIRE, R_AIRE, R_V, L_V,
                                 T_MIN_LCL, SIGMA_TOPE_CONVECCION)
    return T + dT.reshape(T.shape), q + dq.reshape(q.shape), prec.reshape(forma)


# ================================================================
# CONDENSACION DE GRAN ESCALA (Isca lscale_cond, do_evap = falso ✅)
# ================================================================
def condensacion_gran_escala(T, q, pf, dp, gravedad, rh=None):
    """Donde q > rh*q_s(T): se condensa dq = (q - rh*q_s)/(1 + (L/cp)*rh*dq_s/dT) (un paso de Newton,
    como Isca) y el calor latente calienta la capa. Sin reevaporacion. Devuelve (T, q, precipitacion
    kg/m2 del paso, (F, C)). Conserva cp*T + L*q exactamente."""
    rh = RH_CONDENSACION if rh is None else rh
    qs, dqs = qs_y_derivada(T, pf)
    exceso = q - rh * qs
    cond = np.where(exceso > 0.0, exceso / (1.0 + (L_V / CP_AIRE) * rh * dqs), 0.0)
    T = T + (L_V / CP_AIRE) * cond
    q = q - cond
    return T, q, (cond * dp).sum(axis=0) / gravedad


# ================================================================
# LLUVIA O NIEVE
# ================================================================
def fraccion_nieve(T_aire_2m):
    """Fraccion de la precipitacion que cae como nieve (rampa del CLM5 ✅)."""
    return np.clip((T_TODO_LLUVIA - T_aire_2m) / (T_TODO_LLUVIA - T_TODO_NIEVE), 0.0, 1.0)


def beta_cubo(W):
    """Disponibilidad de agua del suelo para evaporar (Manabe 1969)."""
    return np.minimum(1.0, W / (UMBRAL_BETA * W_CAMPO))


def fraccion_cubierta_nieve(S):
    return S / (S + S_CUBRE)


# ================================================================
# ALBEDO SEGUN EL ESPECTRO DE LA ESTRELLA (a + b), interruptor I15
# ================================================================
# La nieve y el hielo reflejan mucho mas en el visible (< 0,7 um) que en el
# infrarrojo cercano. Albedos por bandas de CICE/Icepack (codigo fuente ✅,
# icepack_parameters.F90, esquema CCSM3):
#   nieve fria: visible 0,98, infrarrojo cercano 0,70;
#   nieve fundiendose (1 C antes de la fusion): -0,10 y -0,15;
#   hielo desnudo: 0,78 y 0,36.
# El reparto de la luz de la estrella entre las dos bandas sale de la ley de
# Planck con su temperatura (S3N ~5420 K, Sol 5772 K): es la luz que llega
# arriba de la atmosfera (aproximacion: la atmosfera absorbe algo mas de
# infrarrojo cercano, que en la Tierra sube algo el albedo efectivo).
#   - Nieve: se usa directamente el albedo de banda ancha calculado asi.
#   - Hielo marino: M3N usa 0,65 OBSERVADO en la Tierra (SHEBA); se multiplica
#     por el cociente (albedo con S3N)/(albedo con el Sol) del calculo por
#     bandas, para conservar la calibracion terrestre y añadir solo el efecto
#     de la estrella.
# Comprobacion independiente: Shields et al. (2013, Astrobiology 13, 715,
# tabla 2 ✅) con espectros completos: nieve 0,796 (estrella G) -> 0,748
# (estrella K). Su efecto es mayor que el de dos bandas: este calculo es
# CONSERVADOR (subestima la reduccion). Limitacion declarada.
ALBEDO_BANDAS = {"nieve_fria": (0.98, 0.70), "nieve_fundiendo": (0.88, 0.55), "hielo": (0.78, 0.36)}
LONGITUD_CORTE = 0.7e-6
T_SOL = 5772.0


def fraccion_visible(T_estrella, corte=LONGITUD_CORTE):
    """Fraccion de la emision de un cuerpo negro por debajo de 'corte' (integral de Planck por la serie
    exacta de la funcion de Planck acumulada)."""
    c2 = 1.438776877e-2                                # h*c/k (m K)
    x = c2 / (corte * T_estrella)
    suma = 0.0
    for k in range(1, 200):
        suma += math.exp(-k * x) * (x ** 3 / k + 3 * x ** 2 / k ** 2 + 6 * x / k ** 3 + 6 / k ** 4)
    return 15.0 / math.pi ** 4 * suma                  # fraccion por DEBAJO de 'corte' (x grande = onda corta)


def albedo_banda_ancha(clave, T_estrella):
    f = fraccion_visible(T_estrella)
    vis, nir = ALBEDO_BANDAS[clave]
    return f * vis + (1 - f) * nir


def albedos_estrella(T_estrella):
    """Albedos de banda ancha para la estrella dada: nieve fria, nieve fundiendose y el factor que
    multiplica al albedo del hielo marino observado en la Tierra."""
    return {"nieve_fria": albedo_banda_ancha("nieve_fria", T_estrella),
            "nieve_fundiendo": albedo_banda_ancha("nieve_fundiendo", T_estrella),
            "factor_hielo": albedo_banda_ancha("hielo", T_estrella) / albedo_banda_ancha("hielo", T_SOL)}


# ================================================================
# SUELO: PROPIEDADES TERMICAS SEGUN SU AGUA (interruptor I14)
# ================================================================
# Formulas del CLM5 (nota tecnica ✅; Farouki 1981, de Vries 1963, Johansen),
# para un suelo mineral "medio" (50 % arena, 50 % arcilla) con porosidad 0,43:
#   conductividad seca  l_seco = (0,135 rho_d + 64,7)/(2700 - 0,947 rho_d), rho_d = 2700 (1 - porosidad)
#   solidos             l_s = (8,80 * 50 + 2,92 * 50) / 100;  c_s = (2,128*50 + 2,385*50)/100 * 1e6
#   saturada            l_sat = l_s^(1-porosidad) * l_agua^porosidad (l_agua = 0,57)
#   mezcla              l = Ke l_sat + (1 - Ke) l_seco, Ke = log10(Sr) + 1 (>= 0)
#   capacidad           C = c_s (1 - porosidad) + 4,188e6 * porosidad * Sr
# Sr (saturacion) = W / W_CAMPO del cubo (aproximacion: el cubo representa la
# capa activa del suelo). Resultado: seco l = 0,22 W/m/K, C = 1,29e6 J/m3/K
# (inercia ~530); saturado l = 2,1, C = 3,1e6 (inercia ~2560). Antes (v2.4):
# un unico valor de "roca densa", 2500.
POROSIDAD = 0.43
_RHO_D = 2700.0 * (1 - POROSIDAD)
LAMBDA_SECO = (0.135 * _RHO_D + 64.7) / (2700.0 - 0.947 * _RHO_D)
LAMBDA_SOLIDOS = (8.80 * 50 + 2.92 * 50) / 100
LAMBDA_SAT = LAMBDA_SOLIDOS ** (1 - POROSIDAD) * 0.57 ** POROSIDAD
C_SOLIDOS = (2.128 * 50 + 2.385 * 50) / 100 * 1e6
C_SECO = C_SOLIDOS * (1 - POROSIDAD)
C_AGUA_VOL = 4.188e6


def propiedades_suelo(W):
    """(conductividad W/m/K, capacidad J/m3/K) del suelo con agua W (kg/m2) en el cubo."""
    Sr = np.clip(W / W_CAMPO, 0.0, 1.0)
    Ke = np.where(Sr > 0.1, np.log10(np.maximum(Sr, 1e-12)) + 1.0, 0.0)
    lam = Ke * LAMBDA_SAT + (1 - Ke) * LAMBDA_SECO
    C = C_SECO + C_AGUA_VOL * POROSIDAD * Sr
    return lam, C


# Geometria de las capas del suelo con I14: la de la difusividad del suelo SATURADO (la mayor), para
# que el fondo de la columna quede por debajo de la onda estacional tambien en el suelo humedo
# (con suelo seco la onda penetra menos y la columna sobra un poco, sin error).
K_GEOMETRIA_SUELO = LAMBDA_SAT / (C_SECO + C_AGUA_VOL * POROSIDAD)


def columna_suelo(es_tierra, W, capacidades_base, conductancias_base):
    """Capacidades (J/m2/K) y conductancias (W/m2/K) de las columnas: en tierra, del suelo con su agua
    W (kg/m2); en el agua, las de siempre (capacidades_base, conductancias_base)."""
    from fase2_inercia_multicapa import calcular_limites_capas
    n = capacidades_base.shape[-1]
    lam, C = propiedades_suelo(W)
    limites = np.concatenate(([0.0], calcular_limites_capas(K_GEOMETRIA_SUELO, n)))
    espesores = np.diff(limites)
    distancias = np.diff((limites[:-1] + limites[1:]) / 2.0)
    cap = capacidades_base.copy(); con = conductancias_base.copy()
    for i in range(n):
        cap[..., i] = np.where(es_tierra, C * espesores[i], capacidades_base[..., i])
    for i in range(n - 1):
        con[..., i] = np.where(es_tierra, lam / distancias[i], conductancias_base[..., i])
    return cap, con

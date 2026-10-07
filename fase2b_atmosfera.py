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
#     Louis, Tiedtke y Geleyn (1982), resuelto de forma implicita.
#   - Ajuste convectivo CL <-> TR (Manabe y Wetherald 1967), 6.5 K/km.
#   - Difusion horizontal reubicada: parte atmosferica en la TR, parte
#     oceanica en la capa de mezcla del oceano (Trenberth y Caron 2001).
#   - Altitud dentro de la fisica (aire de referencia mas frio con la
#     altura) y temperatura del aire a 2 m como temperatura de referencia.
#
# DIEZ INTERRUPTORES (INTERRUPTORES_FASE2B: I1-I7 de la Fase 2b, I8 de la
# v2.2c, I9 del hielo marino de la Fase 3, I10 de la atmosfera de varias
# capas de la v3.0; con I10 apagado y el resto encendidos, el modelo
# reproduce la v2.4.2 bit a bit). Con todos apagados, el
# modelo reproduce la Fase 2 (fase2_combinado.py, v2.2) -- es la
# prueba V0 del diseño. Restriccion: sin capacidad de la atmosfera
# (I2 apagado) la atmosfera es "instantanea" (modelo de una capa
# clasico) y solo admite I1; I3, I4 e I5 necesitan I2.

import math
import os
import numpy as np
import scipy.sparse as sp
import fase6_equilibrio as EQ
from scipy.sparse.linalg import splu

from parametros import CONSTANTE_SB, ROTACION_PERIODO, P3N_GRAVEDAD, DURACION_HORA, FACTOR_ROTACION_D, S3N_TEMPERATURA
from parametros import P3N_RADIO, OMEGA_SIDERAL
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
    "calor_sensible": True,         # I3: calor sensible suelo <-> CL (Louis, Tiedtke y Geleyn 1982)
    "capa_limite_radiativa": True,  # I4: la CL absorbe/emite infrarrojo (eps_b > 0)
    "ajuste_convectivo": True,      # I5: ajuste convectivo CL <-> TR
    "difusion_reubicada": True,     # I6: difusion en TR + oceano, no en la piel del suelo
    "altitud_en_fisica": True,      # I7: altitud dentro de la fisica, no a posteriori
    "albedo_oceano_solar": True,    # I8 (v2.2c): albedo del agua segun la altura del sol (Fresnel + Cox-Munk)
    "hielo_marino": True,           # I9 (Fase 3): hielo marino termodinamico (Semtner 1976; Wagner y Eisenman 2015)
    "atmosfera_multicapa": False,   # I10 (v3.0): atmosfera de N capas en coordenada sigma (fase30_multicapa.py, DISENO_V3.0.md).
                                    # APAGADO hasta validarla y calibrarla en modo Tierra (calibrar_v30.py)
    # v3.1 (Fase 5a, fase31_agua.py, DISENO_V3.1.md). Todos APAGADOS hasta validarlos; necesitan I10:
    "ciclo_agua": False,            # I11: vapor por capas, evaporacion, transporte, condensacion, lluvia/nieve, cubo, nieve
    "conveccion_humeda": False,     # I12: Betts-Miller simplificado (Frierson 2007) + ajuste seco al adiabatico seco
    "suelo_termico_agua": False,    # I13: propiedades termicas del suelo segun su agua (CLM5; Farouki 1981)
    "albedo_espectral": False,      # I14: albedo de la nieve y del hielo con el espectro de la estrella
    "vapor_radiativo": False,       # I15 (PROTOTIPO): el infrarrojo depende del vapor del modelo (Byrne y O'Gorman 2013)
    # v3.1-pre11 (Fase 6.3, DISENO_FASE6_3.md §6.7). APAGADO hasta validarlo; necesita I10 e I11:
    "nucleo_dinamico": False,       # I16: nucleo dinamico propio (viento real) en vez de la difusion del aire y del vapor
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
GRADIENTE_ADIABATICO = P3N_GRAVEDAD / CP_AIRE                   # K/m, ~9.0 K/km


def diferencia_critica_cl_tr(gravedad, presion_superficie=PRESION_SUPERFICIE, presion_techo_cl=PRESION_TECHO_CL):
    """
    v2.2c -- CORRIGE un error de la v2.2b. La diferencia maxima de
    temperatura entre la CL y la TR antes de que actue la conveccion debe
    medirse con las MISMAS definiciones de temperatura con las que se
    calibraron las emisividades: temperatura media PONDERADA EN MASA de
    cada capa. En v2.2b se uso 6.5 K/km x (distancia entre los centros
    de masa de las capas) = 28.7 K, pero la media en masa de la TR (que
    incluye la troposfera alta y la estratosfera, mucho mas frias) no es
    la temperatura de su centro de masa. Aqui se calcula como en la
    calibracion: perfil de la Atmosfera Estandar (6.5 K/km hasta la
    tropopausa a 216.65 K, isoterma hasta 20 km y +1 K/km por encima),
    en equilibrio hidrostatico con la gravedad del planeta, y diferencia
    entre las medias en masa de 1000-900 hPa y 900-0 hPa. Tierra: 39.0 K.
    """
    dz = 5.0
    z = np.arange(0, 60000, dz)
    T = np.where(z < 11000, 288.15 - 0.0065 * z, np.where(z < 20000, 216.65, 216.65 + 0.001 * (z - 20000)))
    lnp = math.log(presion_superficie) - np.cumsum(gravedad / (R_AIRE * T)) * dz
    p = np.exp(lnp)
    dm = -np.gradient(p)
    cl = p >= presion_techo_cl
    return float((T[cl] * dm[cl]).sum() / dm[cl].sum() - (T[~cl] * dm[~cl]).sum() / dm[~cl].sum())


DT_CRITICA_CL_TR = diferencia_critica_cl_tr(P3N_GRAVEDAD)

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
VIENTO = 5.0                     # m/s, provisional hasta la Fase 6 (con I16, el viento real del nucleo)
TAU_HIPERDIFUSION_DIAS = 0.5     # I16: amortiguamiento de la onda mas corta (criterio de CESM, DISENO_FASE6_2.md §7)
CORRECTOR_ENERGIA_I16 = True     # v3.1-pre15: corrector global de energia de la dinamica (decision 1.8, §6.13)
KARMAN = 0.4
Z_REF = 10.0                     # m, altura de referencia del aire
Z0M = {TIERRA: 0.01, AGUA: 2e-4}       # suelo desnudo (Carlos 02/10); oceano ~Charnock a 5 m/s
Z0H = {TIERRA: 0.001, AGUA: 2e-5}      # z0h = z0m/10
C_HN_OCEANO = 1.1e-3             # Large y Pond 1982; Smith 1988
LOUIS_B, LOUIS_C, LOUIS_D = 5.0, 5.0, 5.0

DEPURAR = False
# v3.1: funcion opcional f(año, estado) que se llama al acabar cada año (seguimiento de simulaciones largas)
AL_ACABAR_ANO = None
# v3.1-pre12: funcion opcional f(año, paso, pasos_del_año, registrado) que se llama tras cada paso de la fisica
# (barra de progreso de las simulaciones largas; no cambia nada de la simulacion)
AL_PASO = None

# Difusion reubicada (I6), v2.2c: coeficientes PROPIOS de la atmosfera
# (sobre la troposfera) y del oceano (sobre su capa de mezcla). El 0.55
# de la literatura (North 1975; Williams y Kasting 1997) multiplica el
# gradiente de temperatura de SUPERFICIE en un modelo de una sola
# temperatura; aplicado a la troposfera (gradiente ~la mitad) transportaba
# ~1/3 menos de lo debido. Se recalibran repitiendo el procedimiento de
# la literatura con la estructura de M3N: en "modo Tierra" (orbita,
# inclinacion, gravedad y continentes terrestres) se ajustan para
# reproducir Trenberth y Caron (2001): transporte atmosferico maximo
# 5.0 PW (~43 grados) y, a 35 grados, 78 % (norte) / 92 % (sur) del
# total por la atmosfera. Segun Williams y Kasting (1997), D depende de
# presion, composicion y rotacion, no del radio: el mismo valor vale
# para P3N (1 bar, aire tipo Tierra, dia de 24 h).
# v2.4.3: D depende de la rotacion, D proporcional a 1/Omega^2 (Williams y
# Kasting 1997). D_ATMOSFERA es el valor de la Tierra (24 h); la simulacion
# lo multiplica por parametros.FACTOR_ROTACION_D. D_OCEANO no se escala (el
# escalado de Williams y Kasting es para la atmosfera); pendiente de revisar.
#
# RESULTADO DE LA CALIBRACION (02/10/2026, mapa terrestre a 5 grados,
# sin nubes ni humedad):
#   D_atm  D_oc | pico atm. N/S (PW) | total a 35 N/S | % oceano 35 N/S
#   1.2   0.16 |   3.29 / 2.85      |  3.88 / 3.92   |  15 / 27
#   2.4   0.12 |   4.01 / 3.51      |  4.46 / 4.32   |  10 / 19
#   2.4   0.24 |   3.83 / 3.20      |  4.59 / 4.60   |  16 / 31
#   3.6   0.14 |   4.26 / 3.68      |  4.76 / 4.58   |  10 / 20
#   (Trenberth y Caron 2001: pico atm. 5.0 PW; oceano 22 % / 8 % a 35)
# El transporte SATURA: por encima de D_atm ~2.4, subir D casi no mueve
# mas calor (solo aplana la temperatura de la troposfera). El modelo no
# puede llegar a 5.0 PW porque le falta lo que en la Tierra aumenta ese
# transporte: el calor latente del vapor de agua y las nubes (Fase 5).
# Se elige el inicio de la saturacion, D_atm = 2.4 (el transporte
# atmosferico llega al ~80 % del observado), y D_oc = 0.12, que da un
# 14.5 % de transporte oceanico medio a 35 grados (observado: 15 %).
D_ATMOSFERA = 2.4
D_OCEANO = 0.12

# v3.0 (I10): coeficientes de difusion, CALIBRADOS en modo Tierra (05/10/2026,
# calibrar_v30.py, 18 simulaciones; DISENO_V3.0.md seccion 13) frente a
# Trenberth y Caron (2001): la media de los dos hemisferios del transporte
# atmosferico maximo = 5,0 PW y la media de la parte del oceano a 35 grados
# = 15 % (observado 22 % N / 8 % S). Sin vapor de agua todavia: incluyen "de
# prestado" el transporte de calor latente; se recalibran en la v3.1.
#   - D_ATMOSFERA_V30: difusion de la energia estatica seca cerca de la
#     superficie (opcion C', seccion 6), en la Tierra a 24 h. Para P3N se
#     multiplica por parametros.FACTOR_ROTACION_D (D proporcional a 1/Omega^2).
#     (Antes, 0,55 provisional de North 1975.)
#   - D_OCEANO_V30: difusion del oceano con la atmosfera de N capas (con la de
#     dos capas sigue D_OCEANO = 0,12). No se escala con la rotacion.
D_ATMOSFERA_V30 = 1.55
D_OCEANO_V30 = 0.27

# v3.1: cambio de la media anual del cubo (kg/m2) que cuenta como la tolerancia de convergencia
UMBRAL_CAMBIO_CUBO = 1.0
PASOS_VAPOR = 4      # v3.1: el transporte del vapor se hace cada 4 pasos (con un paso 4 veces mayor)
CLAVES_AGUA_V31 = ("lluvia", "nieve", "evap", "escorrentia", "descarga", "recorte", "agua_columna", "cubo", "nieve_suelo", "conv_latente")


def coeficientes_neutros(tipo_superficie):
    """C_N (momento) y C_HN (calor) neutros por celda."""
    z0m = np.where(tipo_superficie == TIERRA, Z0M[TIERRA], Z0M[AGUA])
    z0h = np.where(tipo_superficie == TIERRA, Z0H[TIERRA], Z0H[AGUA])
    c_n = (KARMAN / np.log(Z_REF / z0m)) ** 2
    c_hn_tierra = KARMAN ** 2 / (np.log(Z_REF / z0m) * np.log(Z_REF / z0h))
    c_hn = np.where(tipo_superficie == TIERRA, c_hn_tierra, C_HN_OCEANO)
    return c_n, c_hn, z0m, z0h


def factor_estabilidad_louis(ri, c_n, z0m):
    """Funcion de estabilidad para el calor de Louis, Tiedtke y Geleyn (1982),
    version de capa superficial con b=c=d=5 (la del ECMWF). En v2.2b se
    citaba por error como Louis (1979); la formula no cambia."""
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


# ================================================================
# ALBEDO DEL OCEANO SEGUN LA ALTURA DEL SOL (v2.2c, interruptor I8)
# ================================================================
# Calculado desde la fisica, con la misma base que Jin et al. (2004):
#   - Reflexion de Fresnel del agua (indice de refraccion n = 1.34,
#     luz no polarizada): 2.1 % con el sol en la vertical, ~100 % rasante.
#   - Rugosidad del mar: cada ola es un "espejo" inclinado; la
#     distribucion de pendientes es la medida por Cox y Munk (1954):
#     pendiente cuadratica media = 0.003 + 0.00512 * viento (m/s),
#     isotropa. Se integra la reflexion de Fresnel sobre todas las
#     facetas iluminadas (se desprecian sombras y reflexiones multiples
#     entre olas, igual que en la formulacion basica de Cox-Munk).
#   - Luz difusa (cielo): reflexion de Fresnel integrada sobre un cielo
#     isotropo = 0.0675.
#   - Luz que sale de debajo del agua (dispersion en el oceano): +0.006
#     (Jin et al. 2004).
#   - Fraccion difusa de la luz que llega al suelo: el haz directo se
#     atenua con un espesor optico mayor (TAU_DIRECTO) que la luz total
#     (TAU_NUEVO = 0.18, directa + difusa); la diferencia es la luz que
#     llega dispersada. TAU_DIRECTO se fija para que con el sol en la
#     vertical el 10 % de la luz sea difusa, valor tipico de cielo
#     despejado (Iqbal 1983). Resulta ~19 % difusa con el sol a 30
#     grados de altura y ~65 % con el sol a 6 grados.
# Resultado (media ponderada por la luz, viento 5 m/s): ~0.05 en el
# ecuador y ~0.2 cerca de los polos, frente al 0.08 fijo de antes.
N_AGUA = 1.34
FRACCION_DIFUSA_CENIT = 0.10
TAU_DIRECTO = -math.log((1 - FRACCION_DIFUSA_CENIT) * math.exp(-TAU_NUEVO))   # ~0.285 con TAU_NUEVO = 0.18
ALBEDO_SUBSUPERFICIE = 0.006
_tabla_albedo_directo = None


def _fresnel(cos_i, n=N_AGUA):
    cos_i = np.clip(cos_i, 0.0, 1.0)
    sin_t = np.sqrt(1 - cos_i ** 2) / n
    cos_t = np.sqrt(1 - sin_t ** 2)
    rs = (cos_i - n * cos_t) / (cos_i + n * cos_t)
    rp = (n * cos_i - cos_t) / (n * cos_i + cos_t)
    return 0.5 * (rs ** 2 + rp ** 2)


def _albedo_directo_cox_munk(mu0, viento, n=301, lim=6.0):
    s2 = 0.003 + 0.00512 * viento
    sig = math.sqrt(s2 / 2)
    z = np.linspace(-lim * sig, lim * sig, n)
    zx, zy = np.meshgrid(z, z)
    p = np.exp(-(zx ** 2 + zy ** 2) / (2 * sig ** 2))
    nz = 1 / np.sqrt(1 + zx ** 2 + zy ** 2)
    cos_chi = np.clip((-zx * math.sqrt(1 - mu0 ** 2) + zy * 0.0 + mu0) * nz, 0, None)
    peso = p / nz * cos_chi
    return float((peso * _fresnel(cos_chi)).sum() / peso.sum())


def _albedo_difuso():
    mu = np.linspace(0, 1, 20001)[1:]
    return float(2 * np.trapezoid(_fresnel(mu) * mu, mu))


ALBEDO_DIFUSO_AGUA = _albedo_difuso()   # 0.0675


def tabla_albedo_directo():
    global _tabla_albedo_directo
    if _tabla_albedo_directo is None:
        mus = np.linspace(0.005, 1.0, 200)
        _tabla_albedo_directo = (mus, np.array([_albedo_directo_cox_munk(m, VIENTO) for m in mus]))
    return _tabla_albedo_directo


def albedo_oceano(mu, transmitancia):
    """Albedo del agua para la luz que llega con coseno cenital mu y
    transmitancia atmosferica dada (arrays de las celdas de dia)."""
    mus, tabla = tabla_albedo_directo()
    directo = np.interp(mu, mus, tabla)
    directa = transmitancia ** (TAU_DIRECTO / TAU_NUEVO)      # = exp(-TAU_DIRECTO * masa)
    fraccion_difusa = 1 - directa / np.maximum(transmitancia, 1e-300)
    return (1 - fraccion_difusa) * directo + fraccion_difusa * ALBEDO_DIFUSO_AGUA + ALBEDO_SUBSUPERFICIE


# ================================================================
# FASE 3: HIELO MARINO (interruptor I9). Ver DISENO_FASE3.md.
# ================================================================
# Modelo termodinamico "de capa cero" (Semtner 1976) en la formulacion
# por ENTALPIA de Wagner y Eisenman (2015): cada celda de oceano guarda
# una energia E (J/m2) respecto al agua a punto de congelarse:
#   E >= 0 -> agua liquida a T = T_CONGELACION + E / C0 (capa de mezcla)
#   E <  0 -> hielo de espesor h = -E / (RHO_HIELO * L_FUSION), con el
#             agua de debajo a T_CONGELACION
# Congelar exige quitar calor latente y fundir aportarlo: la energia se
# conserva por construccion. El hielo AISLA: el calor sube por
# conduccion, F = K_HIELO * (T_congelacion - T_superficie) / h, y la
# temperatura de la superficie del hielo sale del balance de energia
# con la atmosfera (luz, infrarrojo, calor sensible), con un maximo de
# 0 C (si sobra energia, el hielo se funde por arriba).
T_CONGELACION = 273.15 - 1.8     # agua de mar ~34 psu (Tf = -0.054 C/psu x S; CAM6 usa -1.8 C)
K_HIELO = 2.034                  # W/m/K, hielo puro (Maykut y Untersteiner 1971; CAM6/CICE)
RHO_HIELO = 917.0                # kg/m3
L_FUSION = 3.34e5                # J/kg
RHO_L_HIELO = RHO_HIELO * L_FUSION   # ~3.06e8 J/m3 (WE15: 9.5 W yr m-3 = 3.0e8)
ALBEDO_HIELO = 0.65              # hielo desnudo (SHEBA, Perovich et al. 2002); sin nieve hasta la Fase 5
ESPESOR_ALBEDO_HIELO = 0.5       # m: por debajo, el hielo fino es mas oscuro (CAM6: valores de hielo "mas grueso de 0.5 m");
                                 # transicion LINEAL agua -> hielo entre 0 y 0.5 m (forma supuesta, documentada)
ESPESOR_MINIMO_CONDUCCION = 0.01 # m: evita dividir por cero con hielo recien formado
ANOS_SALTO_HIELO = (6, 12, 18)   # años en que se acelera el hielo grueso hacia su equilibrio (solo con acelerar=True)
FLUJO_OCEANO_PROFUNDO = 4.0      # W/m2 hacia la base del hielo (Wagner y Eisenman 2015). Se RETIRA en igual
                                 # cantidad, repartida por igual, de todo el oceano (v2.4.1; antes solo del libre):
                                 # redistribuye energia, no la crea.


def nieve_permanente_posible(registro_media_C, es_tierra, dias_mes=None):
    """
    DIAGNOSTICO, sin efecto en la fisica (decision de Carlos, 02/10/2026:
    la nieve en tierra llega con la precipitacion de la Fase 5). Marca
    las celdas de tierra donde, SI nevara, la nieve no llegaria a
    fundirse nunca: el "mes" mas calido tiene media < 0 C, que es la
    definicion del clima de casquete glaciar (EF) de Koppen. El año de
    P3N se divide en 12 "meses" iguales (con el dia de 19,84 h, 326/12 =
    27,2 dias; v2.4.2: se calcula con el numero real de dias del año).
    """
    n = registro_media_C.shape[0]
    dias_mes = n / 12 if dias_mes is None else dias_mes
    medias = np.array([registro_media_C[int(round(k * dias_mes)):int(round((k + 1) * dias_mes))].mean(axis=0) for k in range(12)])
    return es_tierra & (medias.max(axis=0) < 0.0)


def preparar_luz(albedo_grid, tau, reparto_nuevo, mascara_agua=None, albedo_solar=False, factor_masa=None,
                 albedo_hielo=None):
    """
    Devuelve luz(toa, decl, ang_h_lon0) -> (absorbida_suelo, absorbida_atmosfera).
    Misma geometria y masa de aire que preparar_irradiancia_absorbida_
    rejilla() (fase1_geografia.py). reparto_nuevo=False: esquema antiguo
    (lo atenuado se pierde; atmosfera = 0). factor_masa (v3.0): p_superficie/p0
    de cada celda -- sobre una montaña la luz atraviesa menos aire.
    """
    lat_rad = np.radians(LATITUDES_GRADOS).reshape(-1, 1)
    sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
    lon_rad = np.radians(LONGITUDES_GRADOS)
    uno_menos_albedo = 1 - albedo_grid
    if factor_masa is None:
        t_subida = math.exp(-tau * FACTOR_DIFUSIVIDAD)
    else:
        tau_celda = tau * factor_masa
        t_subida_celda = np.exp(-tau_celda * FACTOR_DIFUSIVIDAD)

    a_hielo = ALBEDO_HIELO if albedo_hielo is None else albedo_hielo   # v3.1 (I15): con el espectro de la estrella

    def luz(toa, declinacion, ang_h_lon0, peso_hielo=None, nieve=None, factor_masa_paso=None, viento=None):
        # nieve (v3.1): (fraccion cubierta, albedo de la nieve), arrays (F, C), o None
        # v3.1-pre11 (I16): factor_masa_paso = p_s/p0 de ESTE paso (la p_s del nucleo cambia) y viento (F, C)
        # para el albedo del oceano (Cox y Munk con el viento local). Sin ellos, exactamente como antes.
        if factor_masa_paso is not None:
            tc = tau * factor_masa_paso
            tsc = np.exp(-tc * FACTOR_DIFUSIVIDAD)
        elif factor_masa is not None:
            tc, tsc = tau_celda, t_subida_celda
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
            transmitancia = np.exp(-masa * tau) if (factor_masa is None and factor_masa_paso is None) else np.exp(-masa * tc[dia])
            inst = toa * np.cos(cenital[dia])
            llega = inst * transmitancia
            if albedo_solar:
                if viento is None:
                    alb = np.where(mascara_agua[dia], albedo_oceano(np.cos(cenital[dia]), transmitancia), albedo_grid[dia])
                else:
                    import fase6_superficie as _SUP
                    alb = np.where(mascara_agua[dia], _SUP.albedo_oceano_viento(np.cos(cenital[dia]), transmitancia,
                                                                                  viento[dia]), albedo_grid[dia])
            else:
                alb = albedo_grid[dia]
            if nieve is not None:
                fn = nieve[0][dia]
                alb = np.where(fn > 0, (1 - fn) * alb + fn * nieve[1][dia], alb)
            if peso_hielo is not None:
                w = peso_hielo[dia]
                alb = np.where(w > 0, (1 - w) * alb + w * a_hielo, alb)
                suelo[dia] = llega * (1 - alb)
            elif nieve is not None:
                suelo[dia] = llega * (1 - alb)
            elif albedo_solar:
                suelo[dia] = llega * (1 - alb)
            else:
                suelo[dia] = llega * uno_menos_albedo[dia]
            if reparto_nuevo:
                reflejada = llega * alb
                ts = t_subida if (factor_masa is None and factor_masa_paso is None) else tsc[dia]
                atmos[dia] = FRACCION_ATENUADO_ABSORBIDA * (inst * (1 - transmitancia) + reflejada * (1 - ts))
        return suelo, atmos

    return luz


def simular_fase2b(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
    interruptores=None, K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.015,
    eps_cl=EMISIVIDAD_CL, eps_tr=EMISIVIDAD_TR, d_atmosfera=None, d_oceano=None,
    acelerar=True, estado_inicial="libre", n_capas_atm=None, archivo_estado=None, guardar_al_terminar=False,
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
    if ROTACION_PERIODO % paso_tiempo != 0:
        # v2.4.2: los registros diarios cortan el año en dias de pasos enteros
        raise ValueError(f"El dia ({ROTACION_PERIODO} s) no es un multiplo exacto del paso de tiempo ({paso_tiempo} s)")
    atm = I["capacidad_atmosfera"]
    if not atm and (I["calor_sensible"] or I["capa_limite_radiativa"] or I["ajuste_convectivo"]):
        raise ValueError("calor_sensible, capa_limite_radiativa y ajuste_convectivo necesitan capacidad_atmosfera")
    if I["difusion_reubicada"] and not atm:
        raise ValueError("difusion_reubicada necesita capacidad_atmosfera (la difusion atmosferica actua sobre la TR)")

    multi = I.get("atmosfera_multicapa", False)
    if multi and not all(I[k] for k in ("luz_absorbida", "capacidad_atmosfera", "calor_sensible",
                                          "capa_limite_radiativa", "ajuste_convectivo", "difusion_reubicada",
                                          "altitud_en_fisica")):
        raise ValueError("atmosfera_multicapa (I10) necesita I1-I7 encendidos (sustituye a sus versiones de dos capas)")
    agua_on = I.get("ciclo_agua", False)
    if any(I.get(k, False) for k in ("ciclo_agua", "conveccion_humeda", "suelo_termico_agua", "vapor_radiativo")) and not (multi and agua_on):
        raise ValueError("I11-I13 e I15 necesitan atmosfera_multicapa (I10) y ciclo_agua (I11)")
    if I.get("albedo_espectral", False) and not multi:
        raise ValueError("albedo_espectral (I14) necesita atmosfera_multicapa (I10)")
    din_on = I.get("nucleo_dinamico", False)
    if din_on and not (multi and agua_on):
        raise ValueError("nucleo_dinamico (I16) necesita atmosfera_multicapa (I10) y ciclo_agua (I11)")
    hielo_on = I.get("hielo_marino", False)
    if hielo_on and not (atm and I["calor_sensible"] and I["luz_absorbida"]):
        raise ValueError("hielo_marino necesita capacidad_atmosfera, calor_sensible y luz_absorbida")

    sigma = CONSTANTE_SB
    es_tierra = tipo_superficie == TIERRA
    es_agua = ~es_tierra
    albedo_grid = np.where(es_tierra, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    tau = TAU_NUEVO if I["luz_absorbida"] else profundidad_optica
    if I["albedo_oceano_solar"] and not I["luz_absorbida"]:
        raise ValueError("albedo_oceano_solar necesita luz_absorbida (usa el reparto de la luz de Wild 2019)")
    if multi:
        import fase30_multicapa as M30
        import fase31_agua as AG
        # v3.1: con la conveccion humeda (I12), el ajuste seco va al adiabatico SECO g/cp (fisica, a);
        # sin ella, al gradiente provisional de 6,5 K/km de la v3.0
        COL = M30.Columna(altitud_metros, gravedad=P3N_GRAVEDAD, n=n_capas_atm or M30.N_CAPAS_ATM,
                          gradiente=(P3N_GRAVEDAD / CP_AIRE) if I.get("conveccion_humeda", False) else M30.GRADIENTE_CRITICO,
                          presion_capa="sb81" if I.get("nucleo_dinamico", False) else "media")   # I16: decision 1.2
        ALB_ESTRELLA = AG.albedos_estrella(S3N_TEMPERATURA if I.get("albedo_espectral", False) else AG.T_SOL)
        luz = preparar_luz(albedo_grid, tau, I["luz_absorbida"], ~es_tierra, I["albedo_oceano_solar"],
                           factor_masa=COL.ps / M30.P0,
                           albedo_hielo=ALBEDO_HIELO * ALB_ESTRELLA["factor_hielo"] if I.get("albedo_espectral", False) else None)
    else:
        luz = preparar_luz(albedo_grid, tau, I["luz_absorbida"], ~es_tierra, I["albedo_oceano_solar"])
    eps_b = eps_cl if I["capa_limite_radiativa"] else 0.0
    eps_t = eps_tr if I["capa_limite_radiativa"] else emisividad

    capacidades_reales, conductancias = construir_columna_rejilla(tipo_superficie, inercia_por_tipo, K_difusividad, n_capas)

    # ---- difusion horizontal ----
    D_grid = np.full((FILAS, COLUMNAS), D)
    L_total = construir_matriz_difusion(D_grid)
    if I["difusion_reubicada"]:
        # v2.4.3: D_ATMOSFERA esta calibrado en la Tierra (24 h); se escala con la rotacion de P3N
        d_atm = D_ATMOSFERA * FACTOR_ROTACION_D if d_atmosfera is None else d_atmosfera
        d_oc = (D_OCEANO_V30 if multi else D_OCEANO) if d_oceano is None else d_oceano
        if not multi:     # v3.1: con I10 no se usa (ahorra la factorizacion)
            L_tr = construir_matriz_difusion(np.full((FILAS, COLUMNAS), d_atm))
            fact_tr = splu((sp.diags(np.full(FILAS * COLUMNAS, CAPACIDAD_TR / paso_tiempo)) - L_tr).tocsc())
        L_oc = construir_matriz_difusion_enmascarada(d_oc, ~es_tierra)
    if multi:
        d_atm30 = (D_ATMOSFERA_V30 * FACTOR_ROTACION_D) if d_atmosfera is None else d_atmosfera
        L_atm30 = construir_matriz_difusion(np.full((FILAS, COLUMNAS), d_atm30))
        cap_tr = COL.cap_transporte.flatten()
        fact_atm30 = splu((sp.diags(cap_tr / paso_tiempo) - L_atm30).tocsc())
        gz_cp = (P3N_GRAVEDAD * np.maximum(altitud_metros, 0.0) / CP_AIRE)

    def preparar_superficie(capacidades):
        """Todo lo que depende de la capacidad de la superficie."""
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

    # Estado del hielo (Fase 3): espesor h (m) y temperatura de su superficie.
    HIELO = {"h": np.zeros((FILAS, COLUMNAS)), "Ts": np.full((FILAS, COLUMNAS), T_CONGELACION)}
    if isinstance(estado_inicial, dict):
        # v3.1-pre12: ARRANQUE CALIENTE (DISENO_FASE6_3.md §6.3 y §6.10): el oceano, el suelo y el hielo salen
        # de otra simulacion en equilibrio (el modelo de 2 capas con el mismo mapa y los mismos parametros); el
        # aire arranca en reposo con el perfil inicial de siempre, ahora sobre esa superficie.
        T_col0 = np.asarray(estado_inicial["T_col"], dtype=float)
        if T_col0.shape != T_col.shape:
            raise ValueError(f"arranque caliente: T_col {T_col0.shape} no encaja con esta simulacion {T_col.shape}")
        T_col = T_col0.copy()
        HIELO = {"h": np.array(estado_inicial["HIELO"]["h"], dtype=float),
                 "Ts": np.array(estado_inicial["HIELO"]["Ts"], dtype=float)}
        Ts = T_col[..., 0].copy()
        T_tr = (Ts ** 4 / 2) ** 0.25
        T_cl = Ts - desfase_theta
        # el aire parte de la temperatura REAL de la superficie (la del hielo donde lo hay)
        T_sup_inicial = np.where(es_agua & (HIELO["h"] > 0), HIELO["Ts"], T_col[..., 0]) if hielo_on else T_col[..., 0]
    if hielo_on and estado_inicial == "helado":
        # Prueba de biestabilidad: TODO el oceano cubierto de 5 m de hielo
        # y el planeta frio (agua a punto de congelarse, aire acorde).
        HIELO["h"] = np.where(es_agua, 5.0, 0.0)
        T_col[es_agua] = T_CONGELACION
        T_tr = np.minimum(T_tr, 230.0); T_cl = np.minimum(T_cl, 250.0)
    if multi:
        # v3.0: el "T_cl" pasa a ser la atmosfera entera (N capas, F, C) y el
        # "T_tr", su media ponderada en masa (solo para el criterio de
        # convergencia; se recalcula en cada paso).
        T_cl = COL.perfil_inicial(T_sup_inicial if isinstance(estado_inicial, dict) else T_col[..., 0])
        if hielo_on and estado_inicial == "helado":
            T_cl = np.minimum(T_cl, 250.0)
        T_tr = COL.media_masa(T_cl)
        CAPAS_CL = COL.sigma_media > 0.9          # para los diagnosticos "capa limite" / "troposfera"

    # ---- v3.1: estado del agua (I11) ----
    AGUA_EST = None
    if agua_on:
        # vapor inicial: 60 % de la saturacion en la troposfera (sigma > 0,3), casi seco por encima; el
        # suelo, con el cubo a la mitad; sin nieve. El equilibrio se alcanza solo en pocas semanas.
        qs0, _ = AG.qs_y_derivada(T_cl, COL.pm)
        AGUA_EST = {"q": np.where(COL.sigma_media[:, None, None] > 0.3, 0.6 * qs0, np.minimum(0.6 * qs0, 3e-6)),
                "W": np.where(es_tierra, 0.5 * AG.W_CAMPO, 0.0), "S": np.zeros((FILAS, COLUMNAS)),
                "T_sup": T_col[..., 0].copy(), "contador_vapor": 0,
                "W_suma": np.zeros((FILAS, COLUMNAS)), "n_suma": 0}
        CAPACIDADES_BASE = (capacidades_reales.copy(), conductancias.copy())
        if I.get("suelo_termico_agua", False):
            capacidades_reales, conductancias = AG.columna_suelo(es_tierra, AGUA_EST["W"], *CAPACIDADES_BASE)
        m_capa_baja = COL.dp[-1] / P3N_GRAVEDAD                  # kg/m2 de aire de la capa mas baja
        CAPAS_TR_IDX = np.where(COL.capas_transporte)[0]
        N_TR = len(CAPAS_TR_IDX)
        M_tr = (COL.cap_transporte / CP_AIRE).flatten()                # kg/m2 de las capas del transporte
        fact_vapor = splu((sp.diags(M_tr / (PASOS_VAPOR * paso_tiempo)) - L_atm30 / CP_AIRE).tocsc())
        AREA_REL = (np.cos(np.radians(LATITUDES_GRADOS)).reshape(-1, 1) * np.ones((1, COLUMNAS)))
        AREA_OCEANO_REL = (AREA_REL * es_agua).sum()

    # ---- v3.1-pre11: nucleo dinamico (I16), DISENO_FASE6_3.md §6.7 ----
    DIN = None
    if din_on:
        from fase6_nucleo import NucleoSeco
        import fase6_capa_limite as CL
        import fase6_superficie as SUP
        from fase6_acoplamiento import viento_en_centros, tendencia_a_caras
        from fase6_trazadores import transportar_3d
        NUC = NucleoSeco(P3N_RADIO, OMEGA_SIDERAL, P3N_GRAVEDAD, R_AIRE, CP_AIRE, M30.sigma_seminiveles(COL.n),
                         phis=P3N_GRAVEDAD * np.maximum(altitud_metros, 0.0), filas=FILAS)
        assert NUC.Rj.columnas == COLUMNAS and np.allclose(np.degrees(NUC.Rj.phi_c), LATITUDES_GRADOS)
        NUC.preparar_semiimplicito(paso_tiempo / 2)                 # 2 subpasos de la dinamica por paso de la fisica
        NUC.preparar_hiperdifusion(TAU_HIPERDIFUSION_DIAS, calor_rozamiento=True, correccion_presion=True)
        AREA_NUC = NUC.Rj.area
        DSIG3 = NUC.dsig[:, None, None]
        DIN = {"estado0": (COL.ps.copy(), np.array(T_cl, copy=True), np.zeros((COL.n, FILAS, COLUMNAS)),
                           np.zeros((COL.n, FILAS - 1, COLUMNAS))),
               "ant": None, "act": None, "F": None, "m_tr": None, "paridad": True,
               "v_a": np.full((FILAS, COLUMNAS), VIENTO)}

    def dinamica_y_transporte(q):
        """I16: los 2 subpasos del nucleo (con la tendencia de la fisica del paso anterior como forzamiento
        constante, decision 1.1), el transporte del vapor con los flujos de masa de cada subpaso (§5) y la
        conciliacion con la masa del nucleo (el agua se conserva exactamente). Deja la columna con la p_s del
        nivel n (donde vive el vapor: agua exacta) y devuelve la T del nivel n-1 (la fisica se evalua ahi por
        estabilidad, como CAM) y el viento de ese nivel en los centros."""
        F = DIN["F"]
        forz = None if F is None else (lambda ps_, T_, u_, v_: F)
        # balance de energia de la dinamica (diagnostico, §6.9): energia del nivel n antes y despues de los
        # 2 subpasos, menos la del forzamiento de T aplicado (calculada como en energia_total)
        E0 = None
        if DIN["act"] is not None:
            E0 = energia_dinamica(DIN["act"])
            if F is not None:
                E0 = E0 + float(((CP_AIRE * F[1] * paso_tiempo * NUC.dsig[:, None, None] * DIN["act"][0][None]
                                  ).sum(axis=0) / P3N_GRAVEDAD * PESO).sum())
        for _ in range(2):
            if DIN["act"] is None:
                base = DIN["estado0"]
                Fu, Fv, Wv = NUC.flujos_masa(base[0], base[2], base[3])
                DIN["m_tr"] = DSIG3 * base[0][None] * AREA_NUC[None]
                DIN["ant"], DIN["act"] = NUC.arrancar(base, forz)
            else:
                ps_n, _, u_n, v_n = DIN["act"]
                Fu, Fv, Wv = NUC.flujos_masa(ps_n, u_n, v_n)
                DIN["ant"], DIN["act"] = NUC.avanzar(DIN["ant"], DIN["act"], forz)
            dts = NUC.si_dt
            q, DIN["m_tr"] = transportar_3d(q, DIN["m_tr"], Fu * dts, Fv * dts, Wv[1:-1] * AREA_NUC[None] * dts,
                                            DIN["paridad"])
            DIN["paridad"] = not DIN["paridad"]
        DIN["residuo"] = 0.0 if E0 is None else (energia_dinamica(DIN["act"]) - E0) / paso_tiempo
        # v3.1-pre15: CORRECTOR GLOBAL DE ENERGIA (decision 1.8; DISENO_FASE6_3.md §6.13). Lo que la dinamica y el
        # acoplamiento no conservan en este paso (su residuo, incluida la energia cinetica que quita el rozamiento,
        # mas el calor de rozamiento que devolvio la capa limite) se devuelve como un incremento UNIFORME de T en
        # toda la atmosfera, como el "energy fixer" de CAM (Lauritzen y Williamson 2019, JAMES). Asi el planeta
        # no pierde energia por la numerica y el balance en el tope (N) puede tender a 0 en el equilibrio.
        DIN["fijado"] = 0.0
        if E0 is not None and CORRECTOR_ENERGIA_I16:
            error = DIN["residuo"] + DIN.get("roz", 0.0)                         # W (suma ponderada como acum['abs'])
            masa = float((DIN["act"][0] / P3N_GRAVEDAD * PESO).sum())             # kg (suma ponderada)
            dT_fix = -error * paso_tiempo / (CP_AIRE * masa)
            DIN["ant"] = (DIN["ant"][0], DIN["ant"][1] + dT_fix, DIN["ant"][2], DIN["ant"][3])
            DIN["act"] = (DIN["act"][0], DIN["act"][1] + dT_fix, DIN["act"][2], DIN["act"][3])
            DIN["fijado"] = -error
        ps_act = DIN["act"][0]
        m_din = DSIG3 * ps_act[None] * AREA_NUC[None]
        q = q * DIN["m_tr"] / m_din
        DIN["m_tr"] = m_din
        COL.actualizar_ps(ps_act.copy())
        T_fis = np.array(DIN["ant"][1], copy=True)
        u_c, v_c = viento_en_centros(DIN["ant"][2], DIN["ant"][3])
        v_a = np.sqrt(u_c[-1] ** 2 + v_c[-1] ** 2)
        DIN["v_a"] = v_a
        return T_fis, q, u_c, v_c, v_a

    def energia_dinamica(estado):
        """I16: c_p T + energia cinetica + Phi_s p_s/g de un nivel del nucleo (J/m2, suma ponderada como
        acum['abs']); la misma expresion que energia_total."""
        ps_a, T_a, u_a, v_a_ = estado
        dp_a = NUC.dsig[:, None, None] * ps_a[None]
        e = ((CP_AIRE * T_a + NUC.energia_cinetica(u_a, v_a_)) * dp_a).sum(axis=0) / P3N_GRAVEDAD
        e = e + NUC.phis * ps_a / P3N_GRAVEDAD
        return float((e * PESO).sum())

    def luz_extra():
        """I16: p_s del paso y viento de la superficie para la luz. Sin I16, nada (la luz es la de siempre)."""
        if DIN is None:
            return {}
        return {"factor_masa_paso": COL.ps / M30.P0, "viento": DIN["v_a"]}

    def nieve_actual():
        """(fraccion cubierta, albedo) de la nieve en tierra, para la luz (v3.1)."""
        if AGUA_EST is None:
            return None
        f = np.where(es_tierra, AG.fraccion_cubierta_nieve(AGUA_EST["S"]), 0.0)
        fundiendo = np.clip(AGUA_EST["T_sup"] - (AG.T_FUSION - 1.0), 0.0, 1.0)     # CICE: rampa en el ultimo grado
        alb = ALB_ESTRELLA["nieve_fria"] + (ALB_ESTRELLA["nieve_fundiendo"] - ALB_ESTRELLA["nieve_fria"]) * fundiendo
        return f, alb

    def peso_hielo_actual():
        return np.minimum(1.0, HIELO["h"] / ESPESOR_ALBEDO_HIELO) if hielo_on else None

    def paso(T_col, T_cl, T_tr, s_abs, a_abs, acumular=None):
        C0, resolver_conduccion, fact_sup = SUPERFICIE
        Ts = T_col[..., 0]
        F0 = sigma * Ts ** 4
        if atm:
            Eb = sigma * T_cl ** 4
            Et = sigma * T_tr ** 4
            baja_tr = eps_t * Et
            dlr = (1 - eps_b) * baja_tr + eps_b * Eb
            if hielo_on:
                hay_hielo = es_agua & (HIELO["h"] > 0)
                if hay_hielo.any():
                    # Balance de la superficie del hielo (capa cero): resuelve
                    # luz + infrarrojo - sigma*T^4 - H(T) + K/h (Tf - T) = 0
                    # (Newton), con un maximo de 0 C (fusion por arriba).
                    th = T_cl + desfase_theta
                    kh = K_HIELO / np.maximum(HIELO["h"], ESPESOR_MINIMO_CONDUCCION)
                    T = HIELO["Ts"]
                    tm = 0.5 * (T + th)
                    ri = P3N_GRAVEDAD * Z_REF * (th - T) / (tm * VIENTO ** 2)
                    g_h = PRESION_SUPERFICIE / (R_AIRE * th) * CP_AIRE * c_hn * factor_estabilidad_louis(ri, c_n, z0m) * VIENTO
                    for _ in range(4):
                        f = s_abs + dlr - sigma * T ** 4 - g_h * (T - th) + kh * (T_CONGELACION - T)
                        T = T - f / (-4 * sigma * T ** 3 - g_h - kh)
                    T = np.minimum(T, 273.15)
                    T = np.where(hay_hielo, T, T_CONGELACION)
                    HIELO["Ts"] = T
                    flujo_h_hielo = np.where(hay_hielo, g_h * (T - th), 0.0)
                    F0 = np.where(hay_hielo, sigma * T ** 4, F0)
                    F_superior = s_abs + dlr - F0 - flujo_h_hielo
                    T_cl = T_cl + (paso_tiempo / CAPACIDAD_CL) * flujo_h_hielo
                    # Flujo del oceano profundo a la base del hielo. v2.4.1: la misma
                    # cantidad total se retira de TODO el oceano por igual (hielo
                    # incluido), no solo del libre: asi se conserva la energia
                    # siempre, tambien con el oceano entero helado, y ninguna
                    # celda libre recibe una retirada desproporcionada cuando
                    # quedan pocas (antes: miles de W/m2 en ese caso limite).
                    compensacion = FLUJO_OCEANO_PROFUNDO * (PESO * hay_hielo).sum() / max((PESO * es_agua).sum(), 1e-12)
                else:
                    hay_hielo = None
            sube_cl = (1 - eps_b) * F0 + eps_b * Eb
            olr = (1 - eps_t) * sube_cl + eps_t * Et
            neto_s = s_abs + dlr - F0
            if hielo_on and hay_hielo is not None:
                neto_s = np.where(hay_hielo, F_superior + FLUJO_OCEANO_PROFUNDO, neto_s)
                neto_s = np.where(es_agua, neto_s - compensacion, neto_s)
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
            if hielo_on:
                g = np.where(es_agua & (HIELO["h"] > 0), 0.0, g)
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

        if hielo_on:
            # Entalpia del oceano: E = C0 (T0 - Tf) - rho L h  ->  nuevo (T0, h)
            E = np.where(es_agua, C0 * (T_col[..., 0] - T_CONGELACION) - RHO_L_HIELO * HIELO["h"], 0.0)
            h_nuevo = np.where(es_agua, np.maximum(0.0, -E / RHO_L_HIELO), 0.0)
            nace = (h_nuevo > 0) & (HIELO["h"] <= 0)
            HIELO["Ts"] = np.where(nace, T_CONGELACION, HIELO["Ts"])
            HIELO["h"] = h_nuevo
            T_col[..., 0] = np.where(es_agua, T_CONGELACION + np.maximum(E, 0.0) / C0, T_col[..., 0])
        return T_col, T_cl, T_tr

    def paso_multi(T_col, T_atm, _media, s_abs, a_abs, acumular=None):
        """v3.0 (I10): mismo orden que paso(), con la atmosfera de N capas."""
        C0, resolver_conduccion, fact_sup = SUPERFICIE
        T_atm = T_atm.copy()
        cap1 = COL.cap[-1]
        Ts = T_col[..., 0]
        F0 = sigma * Ts ** 4
        T_aire = COL.aire_superficie(T_atm)
        D_ir, Bh = COL.infrarrojo_bajada(T_atm, T_aire)
        dlr = D_ir[-1]
        hay_hielo = None
        if hielo_on:
            hay_hielo = es_agua & (HIELO["h"] > 0)
            if hay_hielo.any():
                th = T_aire
                kh = K_HIELO / np.maximum(HIELO["h"], ESPESOR_MINIMO_CONDUCCION)
                T = HIELO["Ts"]
                tm = 0.5 * (T + th)
                ri = P3N_GRAVEDAD * Z_REF * (th - T) / (tm * VIENTO ** 2)
                g_h = COL.ps / (R_AIRE * th) * CP_AIRE * c_hn * factor_estabilidad_louis(ri, c_n, z0m) * VIENTO
                for _ in range(4):
                    f = s_abs + dlr - sigma * T ** 4 - g_h * (T - th) + kh * (T_CONGELACION - T)
                    T = T - f / (-4 * sigma * T ** 3 - g_h - kh)
                T = np.minimum(T, 273.15)
                T = np.where(hay_hielo, T, T_CONGELACION)
                HIELO["Ts"] = T
                flujo_h_hielo = np.where(hay_hielo, g_h * (T - th), 0.0)
                F0 = np.where(hay_hielo, sigma * T ** 4, F0)
                F_superior = s_abs + dlr - F0 - flujo_h_hielo
                T_atm[-1] = T_atm[-1] + (paso_tiempo / cap1) * flujo_h_hielo
                compensacion = FLUJO_OCEANO_PROFUNDO * (PESO * hay_hielo).sum() / max((PESO * es_agua).sum(), 1e-12)
            else:
                hay_hielo = None
        U_ir = COL.infrarrojo_subida(Bh, F0)
        olr = U_ir[0]
        neto_s = s_abs + dlr - F0
        if hielo_on and hay_hielo is not None:
            neto_s = np.where(hay_hielo, F_superior + FLUJO_OCEANO_PROFUNDO, neto_s)
            neto_s = np.where(es_agua, neto_s - compensacion, neto_s)
        calor = COL.calentamiento_infrarrojo(U_ir, D_ir) + a_abs[None] * COL.peso_sw
        T_atm = T_atm + paso_tiempo * calor / COL.cap
        if acumular is not None:
            acumular["abs"] += np.sum((s_abs + a_abs) * PESO)
            acumular["olr"] += np.sum(olr * PESO)
            acumular["toa_neto"] += s_abs + a_abs - olr
            acumular["olr_celda"] += olr
            acumular["dlr"] += dlr
            acumular["sw_suelo"] += s_abs
            acumular["sw_atm"] += a_abs
        T_col = T_col.copy()
        T_col[..., 0] = Ts + (paso_tiempo / C0) * neto_s

        # transporte horizontal: oceano (como antes) y atmosfera (opcion C', seccion 6 de DISENO_V3.0.md):
        # difusion implicita de X = s/cp = T_aire + g*z/cp, con la capacidad de las capas que reciben lo
        # que llega; el cambio de X se suma por igual a esas capas (reparto proporcional a la masa).
        T_antes = T_col[..., 0].copy()
        T_col[..., 0] = fact_sup.solve((C0.flatten() / paso_tiempo) * T_col[..., 0].flatten()).reshape(FILAS, COLUMNAS)
        X = (COL.aire_superficie(T_atm) + gz_cp).flatten()
        X_nuevo = fact_atm30.solve((cap_tr / paso_tiempo) * X)
        if acumular is not None:
            # convergencia del transporte horizontal (W/m2) de cada celda, para el diagnostico
            acumular["conv_oceano"] += C0 * (T_col[..., 0] - T_antes) / paso_tiempo
            acumular["conv_atmosfera"] += (cap_tr * (X_nuevo - X)).reshape(FILAS, COLUMNAS) / paso_tiempo
            acumular["n"] += 1
        # misma subida de temperatura en todas las capas que reciben: energia = cap_tr * dX, exacta
        T_atm[COL.capas_transporte] += (X_nuevo - X).reshape(FILAS, COLUMNAS)[None]

        # calor sensible suelo <-> capa mas baja (implicito, cerrado, conserva energia)
        Ts = T_col[..., 0]
        theta_aire = COL.aire_superficie(T_atm)
        dif = Ts - theta_aire
        theta_media = 0.5 * (Ts + theta_aire)
        ri = P3N_GRAVEDAD * Z_REF * (theta_aire - Ts) / (theta_media * VIENTO ** 2)
        c_h = c_hn * factor_estabilidad_louis(ri, c_n, z0m)
        rho = COL.ps / (R_AIRE * theta_aire)
        g = rho * CP_AIRE * c_h * VIENTO
        if hielo_on:
            g = np.where(es_agua & (HIELO["h"] > 0), 0.0, g)
        a_s = C0 / paso_tiempo
        a_b = cap1 / paso_tiempo / COL.factor_superficie     # el aire de superficie cambia factor_superficie veces lo que la capa
        dif_nueva = dif / (1 + g * (1 / a_s + 1 / a_b))
        flujo = g * dif_nueva
        T_col[..., 0] = Ts - flujo / a_s
        T_atm[-1] = T_atm[-1] + flujo / (cap1 / paso_tiempo)

        # ajuste convectivo seco de toda la columna (conserva la energia)
        T_atm = COL.ajuste_convectivo(T_atm)

        T_col = resolver_conduccion(T_col)
        if hielo_on:
            E = np.where(es_agua, C0 * (T_col[..., 0] - T_CONGELACION) - RHO_L_HIELO * HIELO["h"], 0.0)
            h_nuevo = np.where(es_agua, np.maximum(0.0, -E / RHO_L_HIELO), 0.0)
            nace = (h_nuevo > 0) & (HIELO["h"] <= 0)
            HIELO["Ts"] = np.where(nace, T_CONGELACION, HIELO["Ts"])
            HIELO["h"] = h_nuevo
            T_col[..., 0] = np.where(es_agua, T_CONGELACION + np.maximum(E, 0.0) / C0, T_col[..., 0])
        if acumular is not None:
            acumular["T_atm"] += T_atm
        return T_col, T_atm, COL.media_masa(T_atm)

    def paso_agua(T_col, T_atm, _media, s_abs, a_abs, acumular=None):
        """v3.1 (I11-I14): el paso de la v3.0 (paso_multi) con el ciclo del agua. Orden:
        radiacion (y hielo, con sublimacion) -> superficie -> transporte (oceano; aire seco; vapor) ->
        calor sensible y evaporacion -> conveccion humeda -> condensacion de gran escala -> ajuste
        seco (que mezcla tambien el vapor) -> lluvia/nieve -> cubo, nieve, escorrentia, descarga
        glaciar -> conduccion -> hielo. Conserva la energia (cp*T + L*q + entalpias de la superficie,
        del hielo y de la nieve) y el agua."""
        C0, resolver_conduccion, fact_sup = SUPERFICIE
        dt = paso_tiempo
        T_atm = T_atm.copy()
        q = AGUA_EST["q"].copy()
        W = AGUA_EST["W"].copy()
        S = AGUA_EST["S"].copy()
        cap1 = COL.cap[-1]
        m1 = m_capa_baja
        if DIN is not None:
            T_atm, q, u_c, v_c, v_a = dinamica_y_transporte(q)
            cap1 = COL.cap[-1]
            m1 = COL.dp[-1] / P3N_GRAVEDAD
            T_fis0 = T_atm.copy()
            z_a = CL.alturas(T_atm, COL.ph, COL.pm, R_AIRE, P3N_GRAVEDAD)[0][-1]
            flujo_h_hielo = np.zeros((FILAS, COLUMNAS))
        if I.get("vapor_radiativo", False):
            COL.actualizar_tau_vapor(q, AG.A_LW_SECO, AG.B_LW_VAPOR)
        Ts = T_col[..., 0]
        F0 = sigma * Ts ** 4
        T_aire = COL.aire_superficie(T_atm)
        D_ir, Bh = COL.infrarrojo_bajada(T_atm, T_aire)
        dlr = D_ir[-1]
        cero = np.zeros((FILAS, COLUMNAS))
        E_hielo = cero
        hay_hielo = None
        if hielo_on:
            hay_hielo = es_agua & (HIELO["h"] > 0)
            if hay_hielo.any():
                th = T_aire
                kh = K_HIELO / np.maximum(HIELO["h"], ESPESOR_MINIMO_CONDUCCION)
                T = HIELO["Ts"]
                tm = 0.5 * (T + th)
                rho_h = COL.ps / (R_AIRE * th)
                if DIN is None:
                    ri = P3N_GRAVEDAD * Z_REF * (th - T) / (tm * VIENTO ** 2)
                    fe = factor_estabilidad_louis(ri, c_n, z0m)
                    g_h = rho_h * CP_AIRE * c_hn * fe * VIENTO
                    g_qh = rho_h * c_hn * fe * VIENTO                      # kg/m2/s por unidad de q
                else:
                    # I16: coeficientes a z_a con el viento real. ⚠️ Sobre el hielo, la rugosidad del oceano
                    # (Charnock), como la v3.1 usaba los coeficientes del oceano: pendiente de una fuente
                    cf = SUP.coeficientes(z_a, v_a, th, T, np.ones((FILAS, COLUMNAS), bool), Z0M[TIERRA], Z0H[TIERRA],
                                          P3N_GRAVEDAD)
                    g_h = rho_h * CP_AIRE * cf["c_h"] * v_a
                    g_qh = rho_h * cf["c_h"] * v_a
                q_a = q[-1]
                for _ in range(4):
                    qsi, dqsi = AG.qs_y_derivada(T, COL.ps, hielo=True)
                    f = (s_abs + dlr - sigma * T ** 4 - g_h * (T - th) + kh * (T_CONGELACION - T)
                         - AG.L_S * g_qh * (qsi - q_a))
                    T = T - f / (-4 * sigma * T ** 3 - g_h - kh - AG.L_S * g_qh * dqsi)
                T = np.minimum(T, 273.15)
                T = np.where(hay_hielo, T, T_CONGELACION)
                HIELO["Ts"] = T
                qsi, _ = AG.qs_y_derivada(T, COL.ps, hielo=True)
                E_hielo = np.where(hay_hielo, g_qh * (qsi - q_a), 0.0)   # sublimacion (< 0: escarcha)
                flujo_h_hielo = np.where(hay_hielo, g_h * (T - th), 0.0)
                F0 = np.where(hay_hielo, sigma * T ** 4, F0)
                F_superior = s_abs + dlr - F0 - flujo_h_hielo - AG.L_S * E_hielo
                T_atm[-1] = T_atm[-1] + (dt / cap1) * flujo_h_hielo
                q[-1] = q[-1] + E_hielo * dt / m1
                # la masa que se sublima sale del hielo (su entalpia, -L_f por kg, sube sola)
                HIELO["h"] = np.where(hay_hielo, np.maximum(HIELO["h"] - E_hielo * dt / RHO_HIELO, 0.0), HIELO["h"])
                compensacion = FLUJO_OCEANO_PROFUNDO * (PESO * hay_hielo).sum() / max((PESO * es_agua).sum(), 1e-12)
            else:
                hay_hielo = None
        U_ir = COL.infrarrojo_subida(Bh, F0)
        olr = U_ir[0]
        neto_s = s_abs + dlr - F0
        if hielo_on and hay_hielo is not None:
            neto_s = np.where(hay_hielo, F_superior + FLUJO_OCEANO_PROFUNDO, neto_s)
            neto_s = np.where(es_agua, neto_s - compensacion, neto_s)
        calor = COL.calentamiento_infrarrojo(U_ir, D_ir) + a_abs[None] * COL.peso_sw
        T_atm = T_atm + dt * calor / COL.cap
        if acumular is not None:
            acumular["abs"] += np.sum((s_abs + a_abs) * PESO)
            acumular["olr"] += np.sum(olr * PESO)
            acumular["toa_neto"] += s_abs + a_abs - olr
            acumular["olr_celda"] += olr
            acumular["dlr"] += dlr
            acumular["sw_suelo"] += s_abs
            acumular["sw_atm"] += a_abs
        T_col = T_col.copy()
        T_col[..., 0] = Ts + (dt / C0) * neto_s

        # ---- transporte horizontal ----
        T_antes = T_col[..., 0].copy()
        T_col[..., 0] = fact_sup.solve((C0.flatten() / dt) * T_col[..., 0].flatten()).reshape(FILAS, COLUMNAS)
        if DIN is None:
            X = (COL.aire_superficie(T_atm) + gz_cp).flatten()
            X_nuevo = fact_atm30.solve((cap_tr / dt) * X)
            T_atm[COL.capas_transporte] += (X_nuevo - X).reshape(FILAS, COLUMNAS)[None]
        else:
            X = X_nuevo = np.zeros(FILAS * COLUMNAS)          # I16: el aire lo mueve el nucleo, no la difusion
        # vapor: cada capa del transporte (sigma >= 0,25, las mismas que reciben el calor seco) difunde su
        # propia humedad con la MISMA difusividad de remolinos que el calor, kappa = D / (cp * M_tr)
        # (M_tr = masa de esas capas): implicito, monotono (nunca negativo) y conservativo capa a capa
        # (en coordenada sigma la fraccion de masa de cada capa es la misma en todas las celdas).
        # Matriz fija: se factoriza una vez (fact_vapor).
        # Se hace cada PASOS_VAPOR pasos con un paso de PASOS_VAPOR*dt (implicito: estable con cualquier
        # paso; el tiempo de mezcla por remolinos es de dias, mucho mayor que esos ~1 h).
        AGUA_EST["contador_vapor"] += 1
        if DIN is not None:
            dW = np.zeros((FILAS, COLUMNAS))                   # I16: el vapor lo transporta el nucleo
        elif AGUA_EST["contador_vapor"] >= PASOS_VAPOR:
            AGUA_EST["contador_vapor"] = 0
            qk = q[CAPAS_TR_IDX].reshape(N_TR, -1).T                   # (celdas, capas)
            qk_nuevo = fact_vapor.solve(np.asfortranarray((M_tr / (PASOS_VAPOR * dt))[:, None] * qk))
            dq_tr = (qk_nuevo - qk).T.reshape(N_TR, FILAS, COLUMNAS)
            q[CAPAS_TR_IDX] = q[CAPAS_TR_IDX] + dq_tr
            dW = (dq_tr * COL.dp[CAPAS_TR_IDX]).sum(axis=0) / P3N_GRAVEDAD
        else:
            dW = np.zeros((FILAS, COLUMNAS))
        recorte = np.zeros((FILAS, COLUMNAS))
        if acumular is not None:
            acumular["conv_oceano"] += C0 * (T_col[..., 0] - T_antes) / dt
            acumular["conv_atmosfera"] += (cap_tr * (X_nuevo - X)).reshape(FILAS, COLUMNAS) / dt + AG.L_V * dW / dt
            acumular["conv_latente"] += AG.L_V * dW / dt
            acumular["n"] += 1

        # ---- calor sensible (implicito, como la v3.0) ----
        Ts = T_col[..., 0]
        theta_aire = COL.aire_superficie(T_atm)
        dif = Ts - theta_aire
        theta_media = 0.5 * (Ts + theta_aire)
        if DIN is None:
            ri = P3N_GRAVEDAD * Z_REF * (theta_aire - Ts) / (theta_media * VIENTO ** 2)
            c_h = c_hn * factor_estabilidad_louis(ri, c_n, z0m)
            viento_s = VIENTO
        else:
            COEF = SUP.coeficientes(z_a, v_a, theta_aire, Ts, es_agua, Z0M[TIERRA], Z0H[TIERRA], P3N_GRAVEDAD)
            c_h = COEF["c_h"]
            viento_s = v_a
        rho = COL.ps / (R_AIRE * theta_aire)
        g = rho * CP_AIRE * c_h * viento_s
        sin_hielo = ~(es_agua & (HIELO["h"] > 0)) if hielo_on else np.ones((FILAS, COLUMNAS), bool)
        g = np.where(sin_hielo, g, 0.0)
        a_s = C0 / dt
        a_b = cap1 / dt / COL.factor_superficie
        dif_nueva = dif / (1 + g * (1 / a_s + 1 / a_b))
        flujo = g * dif_nueva
        T_col[..., 0] = Ts - flujo / a_s
        T_atm[-1] = T_atm[-1] + flujo / (cap1 / dt)

        # ---- evaporacion (oceano libre, suelo) y sublimacion (nieve en tierra) ----
        # E = rho*C_E*U*beta*(q_s(T_sup) - q_aire), C_E = C_H (Frierson 2007; Isca). Implicita en el
        # vapor de la capa baja (DISENO_FASE5A 3.3), explicita en la temperatura de la superficie.
        Ts = T_col[..., 0]
        g_q = np.where(sin_hielo, rho * c_h * viento_s, 0.0)
        qs_l, _ = AG.qs_y_derivada(Ts, COL.ps)
        qs_i, _ = AG.qs_y_derivada(Ts, COL.ps, hielo=True)
        f_n = np.where(es_tierra, AG.fraccion_cubierta_nieve(S), 0.0)
        beta = np.where(es_tierra, AG.beta_cubo(W), 1.0)
        a_l = g_q * (1 - f_n) * beta                              # suelo / oceano (vapor sobre agua)
        a_i = g_q * f_n                                           # nieve (vapor sobre hielo)
        q_a = q[-1]
        qa_n = (q_a + dt / m1 * (a_l * qs_l + a_i * qs_i)) / (1 + dt / m1 * (a_l + a_i))
        E_l = a_l * (qs_l - qa_n)
        E_i = a_i * (qs_i - qa_n)
        # limites: no se evapora mas agua del cubo ni mas nieve de las que hay
        E_l = np.where(es_tierra & (E_l > 0), np.minimum(E_l, W / dt), E_l)
        E_i = np.where(E_i > 0, np.minimum(E_i, S / dt), E_i)
        q[-1] = q_a + (E_l + E_i) * dt / m1
        T_col[..., 0] = T_col[..., 0] - dt / C0 * (AG.L_V * E_l + AG.L_S * E_i)
        W = np.where(es_tierra, W - E_l * dt, 0.0)
        S = S - E_i * dt

        # ---- I16: capa limite (mezcla vertical y rozamiento con la superficie), antes de la conveccion
        # (orden de lo lento a lo rapido, Beljaars et al. 2004, citado en Gross et al. 2018) ----
        if DIN is not None:
            T_sup_bl = np.where(hay_hielo, HIELO["Ts"], T_col[..., 0]) if hay_hielo is not None else T_col[..., 0]
            f2 = lambda x: x.reshape(x.shape[0], -1)
            u1, v1, T1, q1, dBL = CL.paso_capa_limite(
                f2(u_c), f2(v_c), f2(T_atm), f2(COL.ph), f2(COL.pm), T_sup_bl.ravel(), COEF["z0m"].ravel(),
                COEF["z0h"].ravel(), (flujo + flujo_h_hielo).ravel(), dt, P3N_GRAVEDAD, R_AIRE, CP_AIRE, q=f2(q))
            T_atm = T1.reshape(T_atm.shape)
            q = q1.reshape(q.shape)
            F_u, F_v = tendencia_a_caras((u1.reshape(u_c.shape) - u_c) / dt, (v1.reshape(v_c.shape) - v_c) / dt)
            DIN["roz"] = float((dBL["calor_rozamiento"].reshape(FILAS, COLUMNAS) * PESO).sum())
            if acumular is not None:
                acumular["residuo_dinamica"] += DIN["residuo"]
                acumular["calor_rozamiento"] += DIN["roz"]
                acumular["energia_corregida"] += DIN["fijado"]

        # ---- conveccion humeda, condensacion de gran escala y ajuste seco ----
        P = np.zeros((FILAS, COLUMNAS))
        if I.get("conveccion_humeda", False):
            T_atm, q, Pc = AG.conveccion_humeda(T_atm, q, COL.pm, COL.ph, dt, P3N_GRAVEDAD)
            P = P + Pc
        T_atm, q, Pl = AG.condensacion_gran_escala(T_atm, q, COL.pm, COL.dp, P3N_GRAVEDAD)
        P = P + Pl
        T_atm, q = COL.ajuste_convectivo_vapor(T_atm, q)

        # ---- lluvia o nieve (rampa del CLM5 con el aire a 2 m) ----
        T_sup = np.where(hay_hielo, HIELO["Ts"], T_col[..., 0]) if hay_hielo is not None else T_col[..., 0]
        T2m = T_sup + fraccion_2m * (COL.aire_superficie(T_atm) - T_sup)
        Ps = P * AG.fraccion_nieve(T2m)
        Pr = P - Ps
        T_atm[-1] = T_atm[-1] + AG.L_F * Ps / cap1               # congelarse suelta L_f en el aire bajo
        hielo_aqui = (es_agua & (HIELO["h"] > 0)) if hielo_on else np.zeros((FILAS, COLUMNAS), bool)
        oceano_libre = es_agua & ~hielo_aqui
        # nieve sobre el oceano libre: se funde con calor del oceano; sobre el hielo: se suma a su masa
        T_col[..., 0] = np.where(oceano_libre, T_col[..., 0] - AG.L_F * Ps / C0, T_col[..., 0])
        if hielo_on:
            HIELO["h"] = np.where(hielo_aqui, HIELO["h"] + Ps / RHO_HIELO, HIELO["h"])
        W = np.where(es_tierra, W + Pr, 0.0)
        S = np.where(es_tierra, S + Ps, 0.0)
        # fusion de la nieve en tierra: la energia por encima de 0 C del suelo superficial la funde
        fusion = np.where(es_tierra & (S > 0) & (T_col[..., 0] > AG.T_FUSION),
                          np.minimum(S, C0 * (T_col[..., 0] - AG.T_FUSION) / AG.L_F), 0.0)
        T_col[..., 0] = T_col[..., 0] - fusion * AG.L_F / C0
        S = S - fusion
        W = W + fusion
        escorrentia = np.maximum(W - AG.W_CAMPO, 0.0)
        W = W - escorrentia
        descarga = np.where(es_tierra, np.maximum(S - AG.S_MAX, 0.0), 0.0)
        S = S - descarga
        if descarga.any():
            # el hielo que sobra del tope llega al oceano y se funde: cuesta L_f, repartido por todo el oceano
            energia = AG.L_F * (descarga * AREA_REL).sum()
            T_col[..., 0] = np.where(es_agua, T_col[..., 0] - energia / AREA_OCEANO_REL / C0, T_col[..., 0])

        T_col = resolver_conduccion(T_col)
        if hielo_on:
            E = np.where(es_agua, C0 * (T_col[..., 0] - T_CONGELACION) - RHO_L_HIELO * HIELO["h"], 0.0)
            h_nuevo = np.where(es_agua, np.maximum(0.0, -E / RHO_L_HIELO), 0.0)
            nace = (h_nuevo > 0) & (HIELO["h"] <= 0)
            HIELO["Ts"] = np.where(nace, T_CONGELACION, HIELO["Ts"])
            HIELO["h"] = h_nuevo
            T_col[..., 0] = np.where(es_agua, T_CONGELACION + np.maximum(E, 0.0) / C0, T_col[..., 0])

        if DIN is not None:
            # tendencia de la fisica: se aplica como forzamiento en los 2 subpasos del paso siguiente (1.1)
            DIN["F"] = (np.zeros((FILAS, COLUMNAS)), (T_atm - T_fis0) / dt, F_u, F_v)
        AGUA_EST["q"] = q; AGUA_EST["W"] = W; AGUA_EST["S"] = S; AGUA_EST["T_sup"] = T_col[..., 0].copy()
        AGUA_EST["W_suma"] += W; AGUA_EST["n_suma"] += 1
        if acumular is not None:
            acumular["T_atm"] += T_atm
            for k, v in (("lluvia", Pr), ("nieve", Ps), ("evap", (E_l + E_i) * dt + E_hielo * dt),
                         ("escorrentia", escorrentia), ("descarga", descarga), ("recorte", recorte),
                         ("agua_columna", (q * COL.dp).sum(axis=0) / P3N_GRAVEDAD), ("cubo", W), ("nieve_suelo", S)):
                acumular[k] += v
        return T_col, T_atm, COL.media_masa(T_atm)

    if multi:
        paso = paso_multi
    if agua_on:
        paso = paso_agua
        pasos_dia_agua = round(ROTACION_PERIODO / paso_tiempo)

    PESO = np.cos(np.radians(LATITUDES_GRADOS)).reshape(-1, 1) * np.ones((1, COLUMNAS))

    SUPERFICIE = preparar_superficie(capacidades_reales)

    def entalpia(T_col):
        """Energia del oceano respecto al agua a punto de congelarse (J/m2);
        0 en tierra. Solo tiene sentido con el hielo activado."""
        if not hielo_on:
            return np.zeros((FILAS, COLUMNAS))
        return np.where(es_agua, SUPERFICIE[0] * (T_col[..., 0] - T_CONGELACION) - RHO_L_HIELO * HIELO["h"], 0.0)

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
    convergido = False
    W_media_anterior = None
    # v3.1-pre12: con I16, medias globales anuales para el criterio de equilibrio (fase6_equilibrio.py)
    serie_equilibrio = {k: [] for k in EQ.CLAVES} if DIN is not None else None
    equilibrio_valores = {}

    # v3.1: PUNTO DE CONTROL (opcional). Con archivo_estado, al acabar cada año se guarda el estado
    # completo; si el archivo ya existe al empezar (y es de ESTA misma simulacion: mismos
    # interruptores, D, mapa, orbita y paso), la simulacion continua desde el ultimo año guardado.
    # El resultado es el mismo que sin interrupcion (lo comprueba test_v31).
    import pickle, hashlib
    huella_inicial = estado_inicial if not isinstance(estado_inicial, dict) else hashlib.sha256(
        np.ascontiguousarray(estado_inicial["T_col"], dtype=float).tobytes()
        + np.ascontiguousarray(estado_inicial["HIELO"]["h"], dtype=float).tobytes()
        + np.ascontiguousarray(estado_inicial["HIELO"]["Ts"], dtype=float).tobytes()).hexdigest()
    huella_estado = hashlib.sha256(pickle.dumps((
        sorted(I.items()), d_atmosfera, d_oceano, D, n_capas, paso_tiempo, tolerancia_convergencia, acelerar,
        huella_inicial, n_capas_atm, np.asarray(tipo_superficie).tobytes(), np.asarray(altitud_metros).tobytes(),
        len(datos_orbita), float(datos_orbita[0][0]), float(datos_orbita[-1][0])))).hexdigest()
    ano_inicio = 0
    if archivo_estado is not None and os.path.exists(archivo_estado):
        with open(archivo_estado, "rb") as f:
            guardado = pickle.load(f)
        if guardado["huella"] != huella_estado:
            raise ValueError(f"{archivo_estado} es de otra simulacion (otros parametros); borralo o usa otro nombre")
        T_col, T_cl, T_tr = guardado["T_col"], guardado["T_cl"], guardado["T_tr"]
        HIELO = guardado["HIELO"]
        historial_r, cambio_anterior, saltos = guardado["historial_r"], guardado["cambio_anterior"], guardado["saltos"]
        W_media_anterior = guardado["W_media_anterior"]
        serie_equilibrio = guardado.get("serie_equilibrio", serie_equilibrio)
        ano_inicio = guardado["ano"]
        if AGUA_EST is not None:
            AGUA_EST.update(guardado["AGUA_EST"])
            if I.get("suelo_termico_agua", False) and W_media_anterior is not None:
                capacidades_reales, conductancias = AG.columna_suelo(es_tierra, W_media_anterior, *CAPACIDADES_BASE)
                SUPERFICIE = preparar_superficie(capacidades_reales)
        if DIN is not None:
            DIN.update(guardado["DIN"])                      # I16: estado del nucleo y forzamiento pendiente
            if DIN["act"] is not None:
                COL.actualizar_ps(DIN["act"][0].copy())
        if DEPURAR:
            print(f"  -> se continua desde el año {ano_inicio} guardado en {archivo_estado}", flush=True)

    def guardar_estado(ano_hecho, registrado=False):
        if archivo_estado is None:
            return
        estado = {"huella": huella_estado, "ano": ano_hecho, "T_col": T_col, "T_cl": T_cl, "T_tr": T_tr,
                  "convergido": convergido, "registrado": registrado,
                  "HIELO": HIELO, "historial_r": historial_r, "cambio_anterior": cambio_anterior, "saltos": saltos,
                  "W_media_anterior": W_media_anterior, "serie_equilibrio": serie_equilibrio,
                  "AGUA_EST": None if AGUA_EST is None else {k: v for k, v in AGUA_EST.items()},
                  "DIN": None if DIN is None else {k: v for k, v in DIN.items()}}
        temporal = f"{archivo_estado}.{os.getpid()}.tmp"
        with open(temporal, "wb") as f:
            pickle.dump(estado, f)
        os.replace(temporal, archivo_estado)
    # temperatura real de la superficie: la del hielo donde lo hay
    def temp_superficie(T_col):
        if hielo_on:
            return np.where(es_agua & (HIELO["h"] > 0), HIELO["Ts"], T_col[..., 0])
        return T_col[..., 0]

    def nuevo_acum():
        """Acumuladores de un año (los mismos que los del año final con registro)."""
        acum = {"abs": 0.0, "olr": 0.0}
        if multi:
            acum.update({k: np.zeros((FILAS, COLUMNAS)) for k in
                         ("toa_neto", "conv_oceano", "conv_atmosfera", "olr_celda", "dlr", "sw_suelo", "sw_atm")})
            acum.update(T_atm=np.zeros((COL.n, FILAS, COLUMNAS)), n=0)
        if DIN is not None:
            acum.update(residuo_dinamica=0.0, calor_rozamiento=0.0, energia_corregida=0.0)
        if agua_on:
            acum.update({k: np.zeros((FILAS, COLUMNAS)) for k in CLAVES_AGUA_V31})
        return acum
    PESO_OCEANO = float((PESO * es_agua).sum())
    PESO_TIERRA = float((PESO * es_tierra).sum())

    for ano in range(ano_inicio, max_anos):
        Ts0, Tt0 = T_col[..., 0].copy(), T_tr.copy()
        estado0 = (T_col.copy(), T_cl.copy(), T_tr.copy())
        E0 = entalpia(T_col)
        Tsup0 = temp_superficie(T_col).copy()
        h0 = HIELO["h"].copy()
        conduccion_media = np.zeros((FILAS, COLUMNAS))
        if agua_on:
            W0_ano = AGUA_EST["W"].copy()
            AGUA_EST["W_suma"][:] = 0.0; AGUA_EST["n_suma"] = 0
        acum_ano = nuevo_acum() if DIN is not None else None
        suma_eq = {"T2m": 0.0, "hielo": 0.0, "agua_suelo": 0.0}
        for p_ano, (toa, decl, ang) in enumerate(datos_orbita):
            s, a = luz(toa, decl, ang, peso_hielo_actual(), nieve_actual(), **luz_extra())
            if AL_PASO is not None:
                AL_PASO(ano + 1, p_ano, len(datos_orbita), False)
            if acum_ano is None:
                T_col, T_cl, T_tr = paso(T_col, T_cl, T_tr, s, a)
            else:
                T_col, T_cl, T_tr = paso(T_col, T_cl, T_tr, s, a, acum_ano)
                # I16: medias globales del año para el criterio de equilibrio (aire a 2 m como temp_referencia)
                Tsup_e = temp_superficie(T_col)
                suma_eq["T2m"] += float(((Tsup_e + fraccion_2m * (COL.aire_superficie(T_cl) - Tsup_e)) * PESO).sum())
                if hielo_on and PESO_OCEANO > 0:
                    suma_eq["hielo"] += float((PESO * (es_agua & (HIELO["h"] > 0))).sum()) / PESO_OCEANO
                if PESO_TIERRA > 0:
                    suma_eq["agua_suelo"] += float((PESO * np.where(es_tierra, AGUA_EST["W"], 0.0)).sum()) / (
                        PESO_TIERRA * AG.W_CAMPO)
            if hielo_on:
                conduccion_media += np.where(HIELO["h"] > 0, K_HIELO * (T_CONGELACION - HIELO["Ts"])
                                             / np.maximum(HIELO["h"], ESPESOR_MINIMO_CONDUCCION), 0.0)
        conduccion_media /= len(datos_orbita)
        cambio_agua = 0.0
        if agua_on:
            # v3.1: el agua del suelo tambien tiene que estar en equilibrio: cambio de su media anual
            W_media = AGUA_EST["W_suma"] / max(AGUA_EST["n_suma"], 1)
            cambio_agua = np.inf if W_media_anterior is None else float(np.max(np.abs(W_media - W_media_anterior)))
            W_media_anterior = W_media
            if I.get("suelo_termico_agua", False):
                # propiedades termicas del suelo con el agua media del año (I13). La temperatura del suelo NO
                # cambia: el agua que entra o sale lo hace a la temperatura del suelo (su calor sensible no se
                # contabiliza, opcion A de la 5a). Ese cambio de "energia" ocurre solo entre años, nunca dentro
                # del año final, en el que se mide el cierre del balance. (Primera version: se conservaba
                # C*(T - 0 C), lo que en suelos muy frios, como la Antartida, daba saltos de decenas de grados.)
                capacidades_reales, conductancias = AG.columna_suelo(es_tierra, W_media, *CAPACIDADES_BASE)
                SUPERFICIE = preparar_superficie(capacidades_reales)
        salto_hielo_este_ano = False
        if hielo_on and acelerar and (ano + 1) in ANOS_SALTO_HIELO:
            # ACELERACION DEL HIELO GRUESO (> 1 m). Sin verano, su espesor
            # tiende a h_eq = k (Tf - Ts) / F_base (McKay 2000), pero se
            # acerca muy despacio (decadas). Con lo medido este año se estima
            # el calor que llega a la base, F_base = conduccion media -
            # rho*L*dh/dt, y se salta a h_eq (limitado a x0.5-x2 del espesor
            # actual). No es fisica nueva: solo acorta el camino; despues se
            # siguen simulando años y el criterio de convergencia es el mismo.
            h1 = HIELO["h"]
            ritmo = (h1 - h0) / (len(datos_orbita) * paso_tiempo)          # m/s
            f_base = conduccion_media - RHO_L_HIELO * ritmo
            h_medio = 0.5 * (h0 + h1)
            valido = es_agua & (h0 > 1.0) & (h1 > 1.0) & (f_base > 1.0) & (np.abs(h1 - h0) > 0.01)
            h_eq = np.where(valido, conduccion_media * h_medio / np.where(valido, f_base, 1.0), h1)
            HIELO["h"] = np.where(valido, np.clip(h_eq, 0.5 * h1, 2.0 * h1), h1)
            historial_r = []
            cambio_anterior = None
            salto_hielo_este_ano = bool(valido.any())
            if DEPURAR and valido.any():
                print(f"  -> salto del hielo grueso en {valido.sum()} celdas: espesor medio "
                      f"{h1[valido].mean():.2f} -> {HIELO['h'][valido].mean():.2f} m", flush=True)
        if hielo_on:
            # Con hielo, el equilibrio se juzga por la temperatura REAL de la
            # superficie (la del hielo donde lo hay) y, en el agua libre y el
            # hielo FINO (< 1 m), por la entalpia expresada en grados de la
            # capa de mezcla (1 cm de hielo = 0.015 C): asi la extension del
            # hielo tambien tiene que estabilizarse. El ESPESOR del hielo
            # grueso puede seguir creciendo unos cm/año durante siglos (en
            # la noche polar casi permanente de P3N solo lo frena el calor que
            # llega desde abajo), pero apenas cambia ya la temperatura:
            # ver DISENO_FASE3.md, seccion 4.1.
            dsup = np.abs(temp_superficie(T_col) - Tsup0)
            fino = es_agua & ((HIELO["h"] < 1.0) | (h0 < 1.0))
            dsup = np.where(fino, np.maximum(dsup, np.abs(entalpia(T_col) - E0) / SUPERFICIE[0]), dsup)
        else:
            dsup = np.abs(T_col[..., 0] - Ts0)
        cambio = max(np.max(dsup), np.max(np.abs(T_tr - Tt0)) if atm else 0.0)
        if agua_on:
            # 1 kg/m2 de cambio de la media anual del cubo (0,7 % de su capacidad) cuenta como la tolerancia
            cambio = max(cambio, cambio_agua / UMBRAL_CAMBIO_CUBO * tolerancia_convergencia)
        en_equilibrio = False
        if DIN is not None:
            # I16 (§6.3): equilibrio por medias globales anuales en una ventana de 5 años. La aceleracion
            # geometrica de la v3.1 NO se usa con I16: su razon entre cambios anuales celda a celda no tiene
            # sentido con tiempo meteorologico (§6.10).
            n_p = len(datos_orbita)
            serie_equilibrio["N"].append((acum_ano["abs"] - acum_ano["olr"]) / (n_p * float(PESO.sum())))
            serie_equilibrio["T2m"].append(suma_eq["T2m"] / (n_p * float(PESO.sum())))
            serie_equilibrio["hielo"].append(suma_eq["hielo"] / n_p if (hielo_on and PESO_OCEANO > 0) else None)
            serie_equilibrio["agua_suelo"].append(suma_eq["agua_suelo"] / n_p if PESO_TIERRA > 0 else None)
            en_equilibrio, equilibrio_valores = EQ.evaluar(serie_equilibrio)
        if atm and acelerar and cambio >= tolerancia_convergencia and DIN is None:
            if cambio_anterior is not None:
                historial_r.append(cambio / cambio_anterior)
            cambio_anterior = cambio
            if len(historial_r) >= 3 and max(historial_r[-3:]) - min(historial_r[-3:]) < 0.03 and historial_r[-1] < 0.9:
                r = historial_r[-1]
                factor = r / (1 - r)
                E1 = entalpia(T_col)
                T_col = T_col + factor * (T_col - estado0[0])
                if DIN is None:          # I16: el aire y el viento NO se extrapolan (DISENO_FASE6_3.md §6.3)
                    T_cl = T_cl + factor * (T_cl - estado0[1])
                    T_tr = T_tr + factor * (T_tr - estado0[2])
                if agua_on:
                    # el agua del suelo tambien se extrapola (acotada); el vapor no (se reajusta en dias)
                    AGUA_EST["W"] = np.where(es_tierra, np.clip(AGUA_EST["W"] + factor * (AGUA_EST["W"] - W0_ano),
                                                                0.0, AG.W_CAMPO), 0.0)
                if hielo_on:
                    # en el oceano se extrapola la ENTALPIA (agua + hielo juntos)
                    E = E1 + factor * (E1 - E0)
                    C0s = SUPERFICIE[0]
                    h_nuevo = np.where(es_agua, np.maximum(0.0, -E / RHO_L_HIELO), 0.0)
                    HIELO["Ts"] = np.where((h_nuevo > 0) & (HIELO["h"] <= 0), T_CONGELACION, HIELO["Ts"])
                    HIELO["h"] = h_nuevo
                    T_col[..., 0] = np.where(es_agua, T_CONGELACION + np.maximum(E, 0.0) / C0s, T_col[..., 0])
                saltos += 1
                historial_r = []
                cambio_anterior = None
                if DEPURAR:
                    print(f"  -> salto con r={r:.3f} (factor {factor:.2f})", flush=True)
        if DEPURAR:
            ds = dsup; dt_ = np.abs(T_tr - Tt0)
            i = np.unravel_index(np.argmax(ds), ds.shape); j = np.unravel_index(np.argmax(dt_), dt_.shape)
            print(f"  año {ano+1}: suelo {ds.max():.4f} en {i} ({'tierra' if es_tierra[i] else 'agua'}), TR {dt_.max():.4f} en {j}", flush=True)
        if AL_ACABAR_ANO is not None:
            # v3.1: seguimiento opcional de simulaciones largas (no cambia nada de la simulacion)
            AL_ACABAR_ANO(ano + 1, {"T_sup": temp_superficie(T_col).copy(), "T_atm": np.array(T_cl).copy(),
                                    "hielo_h": HIELO["h"].copy(), "cambio": cambio,
                                    "q": None if AGUA_EST is None else AGUA_EST["q"].copy(),
                                    "nieve": None if AGUA_EST is None else AGUA_EST["S"].copy(),
                                    "cubo": None if AGUA_EST is None else AGUA_EST["W"].copy(),
                                    **({} if DIN is None else {"equilibrio": (en_equilibrio, equilibrio_valores),
                                                               "serie_equilibrio": serie_equilibrio})})
        # v2.4.1: el año del salto del hielo no puede ser el ultimo (el estado
        # acaba de cambiar y el cambio medido es de antes del salto).
        if (en_equilibrio if DIN is not None else cambio < tolerancia_convergencia) and not salto_hielo_este_ano:
            anos = ano + 1
            convergido = True
            if guardar_al_terminar:
                guardar_estado(anos)           # el estado en equilibrio, antes del año registrado
            break
        guardar_estado(ano + 1)

    # ---- año final con registro ----
    pasos_dia = round(ROTACION_PERIODO / paso_tiempo)
    acum = {"abs": 0.0, "olr": 0.0}
    if multi:
        acum.update({k: np.zeros((FILAS, COLUMNAS)) for k in
                     ("toa_neto", "conv_oceano", "conv_atmosfera", "olr_celda", "dlr", "sw_suelo", "sw_atm")})
        acum.update(T_atm=np.zeros((COL.n, FILAS, COLUMNAS)), n=0)
    if DIN is not None:
        acum.update(residuo_dinamica=0.0, calor_rozamiento=0.0, energia_corregida=0.0)
    if agua_on:
        acum.update({k: np.zeros((FILAS, COLUMNAS)) for k in CLAVES_AGUA_V31})
        reg_agua = {k: [] for k in CLAVES_AGUA_V31}
        dia_previo = {k: acum[k].copy() for k in CLAVES_AGUA_V31}

    def energia_total(T_col, T_atm):
        """v3.1: energia total del planeta (J/m2, ponderada por area como acum['abs']) respecto a una
        referencia fija: suelo y oceano (respecto a 0 C), hielo marino (-rho*L*h), aire (cp*T), vapor
        (L_v*q) y nieve en tierra (-L_f por kg). Solo sirve su CAMBIO: cierre del balance de energia."""
        caps = capacidades_reales
        e = np.where(es_tierra, (caps * (T_col - 273.15)).sum(axis=-1), caps[..., 0] * (T_col[..., 0] - 273.15))
        if hielo_on:
            e = e - np.where(es_agua, RHO_L_HIELO * HIELO["h"], 0.0)
        if DIN is not None and DIN["act"] is not None:
            # I16: c_p T + energia cinetica del nivel n del nucleo, mas Phi_s p_s/g; mas la energia de la
            # tendencia de la fisica pendiente de aplicar (se aplica en los subpasos del paso siguiente)
            ps_a, T_a, u_a, v_a_ = DIN["act"]
            dp_a = NUC.dsig[:, None, None] * ps_a[None]
            e = e + ((CP_AIRE * T_a + NUC.energia_cinetica(u_a, v_a_)) * dp_a).sum(axis=0) / P3N_GRAVEDAD
            e = e + NUC.phis * ps_a / P3N_GRAVEDAD
            if DIN["F"] is not None:
                e = e + (CP_AIRE * DIN["F"][1] * paso_tiempo * dp_a).sum(axis=0) / P3N_GRAVEDAD
        else:
            e = e + (COL.cap * T_atm).sum(axis=0)
        if agua_on:
            e = e + AG.L_V * (AGUA_EST["q"] * COL.dp).sum(axis=0) / P3N_GRAVEDAD - AG.L_F * AGUA_EST["S"]
        return float((e * PESO).sum())
    E_inicio_final = energia_total(T_col, T_cl) if multi else None
    if agua_on:
        agua_inicio = {"columna": float(((AGUA_EST["q"] * COL.dp).sum(axis=0) / P3N_GRAVEDAD * PESO).sum()),
                       "suelo": float(((AGUA_EST["W"] + AGUA_EST["S"]) * PESO).sum())}
    reg = {k: [] for k in ("min", "media", "max", "s_min", "s_media", "s_max", "cl", "tr", "h")}
    # Fase 4 (v2.4): registro HORARIO del año final -- valor instantaneo a
    # cada hora en punto (hora del meridiano 0) del aire a 2 m y de la
    # superficie, en todas las celdas, para los dias completos del año.
    # Solo si el paso de tiempo divide exactamente una HORA DE P3N (v2.4.3:
    # 1/24 del dia, parametros.DURACION_HORA; antes, 3600 s).
    pasos_hora = round(DURACION_HORA / paso_tiempo) if DURACION_HORA % paso_tiempo == 0 else None
    pasos_registro_horario = (len(datos_orbita) // pasos_dia) * pasos_dia
    horario_aire, horario_sup = [], []

    def temp_referencia(T_col, T_cl):
        Ts = temp_superficie(T_col)
        if multi:
            return Ts + fraccion_2m * (COL.aire_superficie(T_cl) - Ts)   # aire a 2 m (v3.0)
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
                reg["h"].append(h_suma / cnt)
                if agua_on:
                    for k in CLAVES_AGUA_V31:
                        reg_agua[k].append((acum[k] - dia_previo[k]).astype(np.float32))
                        dia_previo[k] = acum[k].copy()
            mn = np.full((FILAS, COLUMNAS), np.inf); mx = -mn; suma = np.zeros((FILAS, COLUMNAS))
            smn = mn.copy(); smx = mx.copy(); ssuma = suma.copy(); cl_suma = suma.copy(); tr_suma = suma.copy(); cnt = 0
            h_suma = suma.copy()
        if pasos_hora and p < pasos_registro_horario and p % pasos_hora == 0:
            horario_aire.append(a_celsius(temp_referencia(T_col, T_cl)).astype(np.float32))
            horario_sup.append(a_celsius(temp_superficie(T_col)).astype(np.float32))
        s, a = luz(toa, decl, ang, peso_hielo_actual(), nieve_actual(), **luz_extra())
        if AL_PASO is not None:
            AL_PASO(anos + 1, p, len(datos_orbita), True)
        T_col, T_cl, T_tr = paso(T_col, T_cl, T_tr, s, a, acum)
        tref = temp_referencia(T_col, T_cl)
        mn = np.minimum(mn, tref); mx = np.maximum(mx, tref); suma = suma + tref
        Ts = temp_superficie(T_col)
        smn = np.minimum(smn, Ts); smx = np.maximum(smx, Ts); ssuma = ssuma + Ts
        if multi:
            cl_suma = cl_suma + COL.media_masa(T_cl, CAPAS_CL); tr_suma = tr_suma + COL.media_masa(T_cl, ~CAPAS_CL)
        else:
            cl_suma = cl_suma + T_cl; tr_suma = tr_suma + T_tr
        h_suma = h_suma + HIELO["h"]; cnt += 1
    if cnt == pasos_dia:
        for clave, valor in (("min", mn), ("max", mx), ("s_min", smn), ("s_max", smx)):
            reg[clave].append(a_celsius(valor))
        reg["media"].append(a_celsius(suma / cnt)); reg["s_media"].append(a_celsius(ssuma / cnt))
        reg["cl"].append(cl_suma / cnt - 273.15); reg["tr"].append(tr_suma / cnt - 273.15)
        reg["h"].append(h_suma / cnt)
        if agua_on:
            for k in CLAVES_AGUA_V31:
                reg_agua[k].append((acum[k] - dia_previo[k]).astype(np.float32))

    energia_cierre = None
    if multi:
        # v3.1: cierre del balance de energia del año final: cambio de la energia total frente a lo que
        # entra menos lo que sale por arriba (debe ser ~0; error de redondeo e integracion)
        energia_cierre = (energia_total(T_col, T_cl) - E_inicio_final - (acum["abs"] - acum["olr"]) * paso_tiempo) / (
            acum["abs"] * paso_tiempo)
    energia_din = {}
    if multi and DIN is not None:
        # I16 (§6.9): lo que la dinamica no conserva (incluye la energia cinetica que quita el rozamiento) y el
        # calor de rozamiento que devuelve la capa limite. Su suma es el error de energia de la dinamica y del
        # acoplamiento; quitandolo, el resto del balance debe cerrar al nivel del redondeo (prueba).
        n_area = len(datos_orbita) * float(PESO.sum())
        energia_din = {
            "residuo_dinamica_W_m2": acum["residuo_dinamica"] / n_area,
            "calor_rozamiento_W_m2": acum["calor_rozamiento"] / n_area,
            "error_dinamica_W_m2": (acum["residuo_dinamica"] + acum["calor_rozamiento"]) / n_area,
            "correccion_energia_W_m2": acum["energia_corregida"] / n_area,
            "cierre_sin_dinamica": energia_cierre - (acum["residuo_dinamica"] + acum["calor_rozamiento"]
                                                     + acum["energia_corregida"]) / acum["abs"],
        }
    agua = None
    if agua_on:
        # agua: totales del año final (kg/m2) y registros diarios. Cierres (kg/m2 de media, ponderados por area):
        #   atmosfera: cambio del vapor = evaporacion - precipitacion (+ recortes, que deben ser 0)
        #   tierra:    cambio del cubo y la nieve = precipitacion - evaporacion - escorrentia - descarga
        n_dias_reg = len(reg_agua["lluvia"])
        suma_area = PESO.sum()
        tot = {k: acum[k] for k in ("lluvia", "nieve", "evap", "escorrentia", "descarga", "recorte", "conv_latente")}
        P_tot = tot["lluvia"] + tot["nieve"]
        col_fin = float(((AGUA_EST["q"] * COL.dp).sum(axis=0) / P3N_GRAVEDAD * PESO).sum())
        suelo_fin = float(((AGUA_EST["W"] + AGUA_EST["S"]) * PESO).sum())
        cierre_atm = (col_fin - agua_inicio["columna"] - float(((tot["evap"] - P_tot + tot["recorte"]) * PESO).sum())) / suma_area
        P_tierra = np.where(es_tierra, P_tot, 0.0)
        # evaporacion en tierra = la de su celda (en tierra no hay hielo marino)
        E_tierra = np.where(es_tierra, tot["evap"], 0.0)
        cierre_tierra = (suelo_fin - agua_inicio["suelo"] - float(((P_tierra - E_tierra - tot["escorrentia"]
                                                                   - tot["descarga"]) * PESO).sum())) / suma_area
        ano_s = len(datos_orbita) * paso_tiempo
        dia_s = pasos_dia * paso_tiempo
        agua = {
            # medias anuales (mm/dia = kg/m2/dia) por celda
            "precipitacion": P_tot / ano_s * 86400, "nieve": tot["nieve"] / ano_s * 86400,
            "evaporacion": tot["evap"] / ano_s * 86400, "escorrentia": tot["escorrentia"] / ano_s * 86400,
            "descarga_glaciar": tot["descarga"] / ano_s * 86400,
            "agua_precipitable": acum["agua_columna"] / max(acum["n"], 1),     # kg/m2
            "cubo_medio": acum["cubo"] / max(acum["n"], 1), "nieve_media": acum["nieve_suelo"] / max(acum["n"], 1),
            "transporte_latente_conv": tot["conv_latente"] / max(acum["n"], 1),  # W/m2 (convergencia)
            # registros diarios: flujos en mm/dia de P3N... convertidos a mm por dia terrestre (86400 s)
            "diario_precipitacion": (np.array(reg_agua["lluvia"]) + np.array(reg_agua["nieve"])) / dia_s * 86400,
            "diario_nieve": np.array(reg_agua["nieve"]) / dia_s * 86400,
            "diario_evaporacion": np.array(reg_agua["evap"]) / dia_s * 86400,
            "diario_cubo": np.array(reg_agua["cubo"]) / pasos_dia, "diario_nieve_suelo": np.array(reg_agua["nieve_suelo"]) / pasos_dia,
            "diario_agua_precipitable": np.array(reg_agua["agua_columna"]) / pasos_dia,
            "cierre_agua_atmosfera_kg_m2": cierre_atm, "cierre_agua_tierra_kg_m2": cierre_tierra,
            "recortes_kg_m2": float((tot["recorte"] * PESO).sum() / suma_area),
            "dias_registrados": n_dias_reg,
        }

    flujos = None
    if multi:
        # v3.0: medias anuales (W/m2) para el transporte implicito (diagnostico, ver
        # fase30_multicapa.transporte_meridional)
        flujos = {k: acum[k] / max(acum["n"], 1) for k in
                  ("toa_neto", "conv_oceano", "conv_atmosfera", "olr_celda", "dlr", "sw_suelo", "sw_atm", "T_atm")}
        # geometria vertical, para los diagnosticos del perfil (presion en el centro de cada capa, Pa)
        flujos["p_capas"] = COL.pm.copy()
        flujos["p_superficie"] = COL.ps.copy()
    if guardar_al_terminar and archivo_estado is not None:
        # v3.1-pre12: el estado al acabar el año registrado tambien va al punto de control (como un año mas):
        # llamando otra vez con max_anos = ese año, se simula directamente el SIGUIENTE año registrado. Asi
        # fase6_clima.py encadena los años de la climatologia (§6.10). El punto de control anterior (el del
        # comienzo de este año) se conserva en "<archivo>.previo" hasta que fase6_clima.py guarde el año.
        if os.path.exists(archivo_estado):
            import shutil
            shutil.copyfile(archivo_estado, archivo_estado + ".previo")
        guardar_estado(anos + 1, registrado=True)
    return {
        "anos": anos, "saltos": saltos, "flujos": flujos, "convergido": convergido,
        # v3.1: coeficientes de difusion que se han usado DE VERDAD (W/m2/K, ya escalados con la rotacion)
        "D_usados": ({"atmosfera": float(d_atm30), "oceano": float(d_oc), "esquema": "v3.0 (N capas)"} if multi else
                     {"atmosfera": float(d_atm), "oceano": float(d_oc), "esquema": "dos capas"} if I["difusion_reubicada"] else
                     {"superficie": float(D), "esquema": "Fase 2"}),
        "T_final": a_celsius(temp_referencia(T_col, T_cl)),
        "reg_min": np.array(reg["min"]), "reg_media": np.array(reg["media"]), "reg_max": np.array(reg["max"]),
        "suelo_min": np.array(reg["s_min"]), "suelo_media": np.array(reg["s_media"]), "suelo_max": np.array(reg["s_max"]),
        "cl_media": np.array(reg["cl"]), "tr_media": np.array(reg["tr"]),
        "hielo_espesor": np.array(reg["h"]),
        "nieve_permanente_posible": nieve_permanente_posible(np.array(reg["media"]), es_tierra),
        # (dia*24 + hora, FILAS, COLUMNAS), en C, instantaneo a la hora en punto del meridiano 0
        "horario_aire2m": np.array(horario_aire) if pasos_hora else None,
        "horario_superficie": np.array(horario_sup) if pasos_hora else None,
        "energia": {"absorbido": acum["abs"], "olr": acum["olr"],
                    "diferencia_relativa": abs(acum["abs"] - acum["olr"]) / acum["abs"],
                    "cierre_relativo": energia_cierre, **energia_din},
        "agua": agua,
        # v3.1-pre12: estado final del oceano, el suelo y el hielo (para arrancar en caliente otra simulacion)
        "estado_final": {"T_col": T_col.copy(), "HIELO": {k: np.array(v, copy=True) for k, v in HIELO.items()}},
        "equilibrio": None if DIN is None else {"serie": serie_equilibrio, "valores": equilibrio_valores},
    }

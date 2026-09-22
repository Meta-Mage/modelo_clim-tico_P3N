# fase2_inercia_multicapa.py -- Fase 2: inercia termica en profundidad,
# version de COLUMNA VERTICAL DE N CAPAS (conduccion real), pensada para
# sustituir/superar a fase2_inercia_estacional.py (la version de 2 capas
# "force-restore").
#
# ================================================================
# POR QUE ESTE MODULO EXISTE (y no solo "3 capas de Boone & Calvet 1999")
# ================================================================
# El modulo de 2 capas (fase2_inercia_estacional.py) usa el metodo
# "force-restore": colapsa toda la columna de suelo en UN resorte
# (constante de acoplamiento k) entre superficie y "profundidad". Es un
# metodo real, con precedente en la literatura (Bhumralkar 1975;
# Blackadar 1976; usado en ISBA), pero con un sesgo sistematico
# documentado cuando debe representar a la vez el ciclo diurno y el
# estacional -- por eso existen extensiones de 3 capas en esa misma
# familia de metodos (Boone & Calvet 1999).
#
# Se intento localizar el texto completo de Boone & Calvet (1999) para
# reproducir sus coeficientes exactos, sin exito: el documento de HAL
# (hal.science/hal-02269589) y la propia revista (AMS) devuelven acceso
# bloqueado (403), y la pagina de Semantic Scholar no ofrece el texto.
# No se ha podido verificar de primera mano su formulacion concreta.
#
# En su lugar, este modulo NO reproduce ese paper -- va un paso mas
# alla y resuelve la fisica subyacente REAL de la que tanto el force-
# restore de 2 capas como el de 3 capas son, a su vez, aproximaciones:
# la ecuacion de conduccion de calor 1D en el suelo,
#
#     c_v * dT/dt = d/dz( k * dT/dz )
#
# discretizada en una columna de N capas (N configurable, no fijo en
# 2 o 3) con diferencias finitas implicitas. Esto es lo que hacen los
# esquemas de superficie terrestre mas rigurosos hoy en dia (CLM,
# Noah-MP) en lugar de force-restore, y evita tener que adivinar
# coeficientes de un paper que no se ha podido verificar: aqui los
# unicos coeficientes son los fisicos de verdad (conductividad termica,
# capacidad calorifica volumetrica), y el numero de capas es un
# parametro de resolucion numerica, no una eleccion fisica.
# ================================================================
#
# DE DONDE SALEN LOS DOS PARAMETROS FISICOS NUEVOS QUE HACEN FALTA
# ================================================================
# La "inercia" que ya tienes tabulada (INERCIA_POR_TIPO[TIERRA] = 2500)
# es la EFUSIVIDAD termica real, I = sqrt(k * c_v) -- una sola ecuacion,
# dos incognitas (k = conductividad, c_v = capacidad calorifica
# volumetrica). No basta para reconstruir una columna: hace falta un
# segundo dato fisico independiente.
#
# Ese segundo dato es la DIFUSIVIDAD termica, K = k / c_v (m2/s) -- una
# propiedad tabulada en cualquier libro de fisica del suelo/geofisica,
# independiente de la inercia. Con las dos ecuaciones:
#
#     I = sqrt(k * c_v)      (ya conocido, fijado en Fase 1)
#     K = k / c_v             (nuevo, de la literatura)
#
# se despejan k y c_v sin inventar nada ni tocar el valor de inercia ya
# validado:
#
#     k   = I * sqrt(K)
#     c_v = I / sqrt(K)
#
# Valor de K adoptado para TIERRA: 1.0e-6 m2/s -- roca densa/suelo
# compacto (rango real tipico en la literatura: ~0.2e-6 a 1.2e-6 m2/s
# para roca; ~0.2e-6 a 0.8e-6 m2/s para distintos suelos). Se elige el
# extremo de roca densa por coherencia con el propio I=2500 ya
# adoptado en Fase 1 (documentado alli como "extremo alto" del rango
# 400-2500 de inercia real, que tambien corresponde a roca densa). La
# consistencia se puede comprobar numericamente: con K=1e-6 e I=2500,
# se obtiene k=2.5 W/m/K y c_v=2.5e6 J/m3/K -- ambos caen exactamente
# en el rango tipico publicado para roca densa (k: 2-3 W/m/K, c_v:
# 2-2.6e6 J/m3/K). Es decir, K=1e-6 no es un numero elegido al azar
# para que "cuadre": es el valor de diffusividad que, combinado con nu
# inercia YA fijada, reproduce conductividad y capacidad calorifica
# fisicamente razonables para el mismo tipo de material.
#
# AGUA sigue FUERA de este esquema (ver AGUA_SIN_CAPA_PROFUNDA en
# fase2_inercia_estacional.py): el oceano no se calienta por
# conduccion sino por mezcla turbulenta, asi que una columna de
# conduccion no es el modelo fisico correcto para el -- se queda con
# el tratamiento de una sola capa ya validado.
# ================================================================
#
# COMO SE ELIGEN LAS PROFUNDIDADES DE LAS N CAPAS
# ================================================================
# Ni al azar ni fijas: se calculan a partir de las dos "profundidades
# de amortiguamiento" fisicas del problema (misma formula d=sqrt(K*P/pi)
# que ya se usa en todo M3N para escalar inercia por periodo, aplicada
# aqui a la difusividad K en vez de a la inercia):
#
#     d_diurno    = sqrt(K * ROTACION_PERIODO / PI)   -- ~0.17 m con K=1e-6
#     d_estacional = sqrt(K * ORBITA_PERIODO  / PI)   -- ~2.7 m con K=1e-6
#
# La primera capa debe ser mucho mas fina que d_diurno (para resolver
# bien el ciclo dia/noche cerca de la superficie); el fondo de la
# columna debe estar varias veces mas profundo que d_estacional (para
# que el ciclo anual ya se haya amortiguado casi del todo antes de
# llegar al fondo, y el limite de flujo cero ahi no distorsione el
# resultado). Los limites entre capas se distribuyen geometricamente
# (cada capa mas gruesa que la anterior, un patron estandar en estos
# esquemas -- concentra resolucion donde el gradiente termico es mayor,
# cerca de la superficie) entre esos dos extremos.
# ================================================================
#
# METODO NUMERICO: igual patron que el resto de Fase 2 (difusion
# horizontal, 2 capas) -- "operator splitting": paso radiativo
# explicito SOLO en la capa superficial (la unica que ve el sol),
# seguido de un paso de conduccion vertical implicito para toda la
# columna. La conduccion implicita en una columna de N capas es un
# sistema TRIDIAGONAL (cada capa solo se comunica con sus vecinas
# inmediatas arriba/abajo) -- se resuelve con el algoritmo de Thomas
# (eliminacion gaussiana especializada para tridiagonales, O(N) en vez
# de O(N^3)), vectorizado sobre toda la rejilla a la vez (cada celda de
# la rejilla tiene su propia columna independiente).

import math
import numpy as np
from parametros import ROTACION_PERIODO, ORBITA_PERIODO, PI, CONSTANTE_SB
from temperatura import precalcular_orbita, estimar_T_inicial_equilibrio_punto, PASO_TIEMPO
from geometria import angulo_cenital
from atmosfera import masa_aire
from radiacion import i_inst, trans, i_atm, i_abs

# ---- Parametros fisicos nuevos (ver cabecera) ----
K_DIFUSIVIDAD_TIERRA = 1.0e-6  # m2/s -- roca densa, coherente con INERCIA_POR_TIPO[TIERRA]=2500

# ---- Resolucion de la columna (parametros numericos, no fisicos) ----
# 6 en vez de 8 (eleccion de Carlos, 21/09): la Prueba 3 del bloque
# __main__ de este archivo mostro que la diferencia N=8 vs N=12 ya es
# pequeña (0.04 C) frente a N=4 vs N=8 (0.29 C) -- rendimientos
# decrecientes claros a partir de N=8, asi que bajar a 6 recupera algo
# de velocidad sin alejarse mucho de la precision de 8.
N_CAPAS_DEFECTO = 6
FRACCION_PRIMERA_CAPA = 0.25   # primera capa = esta fraccion de d_diurno
FACTOR_PROFUNDIDAD_FONDO = 4.5  # fondo de columna = este factor * d_estacional


def profundidad_amortiguamiento(K_difusividad, periodo):
    """d = sqrt(K*periodo/PI) -- misma forma que el resto de M3N usa
    para escalar inercia por periodo, aplicada aqui a la difusividad."""
    return math.sqrt(K_difusividad * periodo / PI)


def derivar_k_y_cv(inercia, K_difusividad):
    """A partir de la efusividad (inercia, ya tabulada) y la
    difusividad (nueva, de la literatura): conductividad k (W/m/K) y
    capacidad calorifica volumetrica c_v (J/m3/K). Ver cabecera."""
    k_conductividad = inercia * math.sqrt(K_difusividad)
    c_v = inercia / math.sqrt(K_difusividad)
    return k_conductividad, c_v


def calcular_limites_capas(K_difusividad, n_capas=N_CAPAS_DEFECTO,
                            fraccion_primera=FRACCION_PRIMERA_CAPA,
                            factor_fondo=FACTOR_PROFUNDIDAD_FONDO):
    """
    Devuelve un array de n_capas limites de profundidad (m), z_1..z_N
    (el limite superior z_0=0, la superficie, no se incluye -- es
    implicito). Distribucion geometrica entre la primera capa (fraccion
    de d_diurno) y el fondo (factor * d_estacional).
    """
    d_diurno = profundidad_amortiguamiento(K_difusividad, ROTACION_PERIODO)
    d_estacional = profundidad_amortiguamiento(K_difusividad, ORBITA_PERIODO)

    z_1 = fraccion_primera * d_diurno
    z_N = factor_fondo * d_estacional

    if n_capas == 1:
        return np.array([z_N])

    razon = (z_N / z_1) ** (1.0 / (n_capas - 1))
    indices = np.arange(n_capas)
    limites = z_1 * razon ** indices
    return limites


def construir_columna(inercia, K_difusividad, n_capas=N_CAPAS_DEFECTO):
    """
    Construye la geometria y propiedades de una columna de n_capas:
    devuelve (capacidades, conductancias), donde:
      - capacidades: array (n_capas,) de C_i = c_v * espesor_i (J/m2/K)
      - conductancias: array (n_capas-1,) de g_i = k / distancia entre
        los centros de las capas i e i+1 (W/m2/K) -- el acoplamiento
        entre cada par de capas consecutivas.
    """
    k_conductividad, c_v = derivar_k_y_cv(inercia, K_difusividad)
    limites = calcular_limites_capas(K_difusividad, n_capas)

    limites_completos = np.concatenate(([0.0], limites))  # incluye la superficie
    espesores = np.diff(limites_completos)                # (n_capas,)
    capacidades = c_v * espesores                          # (n_capas,)

    centros = (limites_completos[:-1] + limites_completos[1:]) / 2.0  # (n_capas,)
    distancias_centro_a_centro = np.diff(centros)          # (n_capas-1,)
    conductancias = k_conductividad / distancias_centro_a_centro  # (n_capas-1,)

    return capacidades, conductancias, limites, centros


def _resolver_tridiagonal_vectorizado(a, b, c, d):
    """
    Algoritmo de Thomas, vectorizado sobre las dimensiones iniciales.
    a, b, c, d tienen forma (..., N): a=subdiagonal (a[...,0] no se usa),
    b=diagonal, c=superdiagonal (c[...,-1] no se usa), d=termino
    independiente. Devuelve T de forma (..., N).
    Resuelve, para cada "..." de forma independiente (cada celda de la
    rejilla o el unico punto):
        a_i*T_{i-1} + b_i*T_i + c_i*T_{i+1} = d_i
    """
    N = b.shape[-1]
    c_prima = np.empty_like(b)
    d_prima = np.empty_like(d)

    c_prima[..., 0] = c[..., 0] / b[..., 0]
    d_prima[..., 0] = d[..., 0] / b[..., 0]

    for i in range(1, N):
        denominador = b[..., i] - a[..., i] * c_prima[..., i - 1]
        if i < N - 1:
            c_prima[..., i] = c[..., i] / denominador
        d_prima[..., i] = (d[..., i] - a[..., i] * d_prima[..., i - 1]) / denominador

    T = np.empty_like(d)
    T[..., N - 1] = d_prima[..., N - 1]
    for i in range(N - 2, -1, -1):
        T[..., i] = d_prima[..., i] - c_prima[..., i] * T[..., i + 1]

    return T


def paso_conduccion_implicito(T_columna, capacidades, conductancias, paso_tiempo):
    """
    Un paso de conduccion vertical implicita (backward Euler) sobre una
    columna ya actualizada por el paso radiativo explicito (ver
    paso_multicapa()). T_columna: forma (..., N). capacidades: forma
    (..., N) (compatible por broadcasting con un punto o una rejilla).
    conductancias: forma (..., N-1).
    """
    N = T_columna.shape[-1]
    forma = T_columna.shape

    if N == 1:
        # Sin capas vecinas, no hay conduccion que resolver -- la unica
        # capa se queda como esta tras el paso radiativo.
        return T_columna.copy()

    a = np.zeros(forma)
    b = np.zeros(forma)
    c = np.zeros(forma)
    d = np.zeros(forma)

    C_sobre_dt = capacidades / paso_tiempo

    # Capa superficial (i=0): solo conduce hacia abajo (g_0).
    g0 = conductancias[..., 0]
    b[..., 0] = C_sobre_dt[..., 0] + g0
    c[..., 0] = -g0
    d[..., 0] = C_sobre_dt[..., 0] * T_columna[..., 0]

    # Capas intermedias.
    for i in range(1, N - 1):
        g_arriba = conductancias[..., i - 1]
        g_abajo = conductancias[..., i]
        a[..., i] = -g_arriba
        b[..., i] = C_sobre_dt[..., i] + g_arriba + g_abajo
        c[..., i] = -g_abajo
        d[..., i] = C_sobre_dt[..., i] * T_columna[..., i]

    # Capa mas profunda (i=N-1): flujo cero en el fondo (aislada).
    g_ultima = conductancias[..., N - 2]
    a[..., N - 1] = -g_ultima
    b[..., N - 1] = C_sobre_dt[..., N - 1] + g_ultima
    d[..., N - 1] = C_sobre_dt[..., N - 1] * T_columna[..., N - 1]

    return _resolver_tridiagonal_vectorizado(a, b, c, d)


def paso_multicapa(T_columna, abs_local, emisividad, capacidades, conductancias, paso_tiempo):
    """Un paso de tiempo completo: radiativo explicito en la capa 0,
    luego conduccion implicita en toda la columna."""
    emitido = (1 - emisividad / 2) * CONSTANTE_SB * T_columna[..., 0] ** 4
    T_columna = T_columna.copy()
    T_columna[..., 0] = T_columna[..., 0] + (paso_tiempo / capacidades[..., 0]) * (abs_local - emitido)
    return paso_conduccion_implicito(T_columna, capacidades, conductancias, paso_tiempo)


def simular_multicapa_punto(
    latitud_rad, datos_orbita, emisividad, inercia, albedo, profundidad_optica,
    K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO,
    paso_tiempo=PASO_TIEMPO, max_anos=50, tolerancia_convergencia=0.01,
):
    """Version de un punto. Devuelve: T_columna final (n_capas,, en C),
    limites de capas (m), registro de temperaturas del ultimo año para
    cada capa (lista de arrays (n_capas,))."""
    capacidades, conductancias, limites, centros = construir_columna(inercia, K_difusividad, n_capas)

    T_inicial = estimar_T_inicial_equilibrio_punto(latitud_rad, datos_orbita, emisividad, albedo, profundidad_optica)
    T_columna = np.full(n_capas, T_inicial)

    registro = []
    for ano in range(max_anos):
        T_inicio_ano = T_columna.copy()
        registro = []
        for toa, decl, ang_h in datos_orbita:
            cenital = angulo_cenital(latitud_rad, decl, ang_h)
            inst = i_inst(toa, cenital)
            masa = masa_aire(cenital)
            if masa is None:
                abs_local = 0.0
            else:
                tra = trans(masa, profundidad_optica)
                atm = i_atm(inst, tra)
                abs_local = i_abs(atm, albedo)

            T_columna = paso_multicapa(T_columna, abs_local, emisividad, capacidades, conductancias, paso_tiempo)
            registro.append(T_columna - 273.15)

        if np.max(np.abs(T_columna - T_inicio_ano)) < tolerancia_convergencia:
            break

    return T_columna - 273.15, limites, centros, registro


def construir_columna_rejilla(tipo_superficie, inercia_por_tipo, K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO):
    """
    Version de rejilla completa de construir_columna(): para cada celda
    de la rejilla, construye su columna de n_capas (capacidades y
    conductancias), vectorizado sobre toda la rejilla a la vez. Tierra
    recibe la columna real (conduccion); agua se queda con su
    tratamiento de una sola capa (ver AGUA_SIN_CAPA_PROFUNDA en
    fase2_inercia_estacional.py, mismo criterio aqui) -- se le monta
    una "columna" cuya capa 0 es su capacidad diurna de siempre y cuyas
    demas capas son de relleno, desacopladas (conductancia 0), para
    poder tratar agua y tierra con el mismo bucle vectorizado sin
    duplicar codigo ni cambiar la fisica del agua.

    Devuelve: capacidades (filas, columnas, n_capas), conductancias
    (filas, columnas, n_capas-1) -- ambas en las unidades ya usadas en
    construir_columna() (J/m2/K y W/m2/K).
    """
    from fase1_geografia import TIERRA, AGUA

    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])
    filas, columnas = tipo_superficie.shape
    es_tierra = (tipo_superficie == TIERRA)

    k_conductividad_grid = inercia_grid * math.sqrt(K_difusividad)
    c_v_grid = inercia_grid / math.sqrt(K_difusividad)

    limites_tierra = calcular_limites_capas(K_difusividad, n_capas)
    limites_completos_tierra = np.concatenate(([0.0], limites_tierra))
    espesores_tierra = np.diff(limites_completos_tierra)
    centros_tierra = (limites_completos_tierra[:-1] + limites_completos_tierra[1:]) / 2.0
    distancias_tierra = np.diff(centros_tierra)

    capacidades = np.zeros((filas, columnas, n_capas))
    conductancias = np.zeros((filas, columnas, n_capas - 1))

    C_diurno_agua_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)
    for i in range(n_capas):
        capacidades[..., i] = np.where(es_tierra, c_v_grid * espesores_tierra[i], C_diurno_agua_grid)
    for i in range(n_capas - 1):
        conductancias[..., i] = np.where(es_tierra, k_conductividad_grid / distancias_tierra[i], 0.0)

    return capacidades, conductancias


def simular_rejilla_multicapa_con_registro(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica,
    K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO,
    paso_tiempo=PASO_TIEMPO, max_anos=50, tolerancia_convergencia=0.01,
):
    """
    Version de rejilla completa, con registro diario (minima/media/
    maxima) de la CAPA SUPERFICIAL de la columna -- la que corresponde
    a lo que ya muestran mapa_calor.py y compañia. El agua se queda con
    el tratamiento de una sola capa (AGUA_SIN_CAPA_PROFUNDA, igual
    criterio que en fase2_inercia_estacional.py): se le monta una
    "columna" de 1 capa (n_capas=1 => sin conduccion, ver
    paso_conduccion_implicito) con su inercia diurna habitual, para
    poder tratar agua y tierra con el mismo bucle vectorizado sin
    duplicar codigo, sin que eso cambie su fisica.

    Devuelve: T_superficie_final, anos_convergencia, registro_minima,
    registro_media, registro_maxima (estas tres, de la capa superficial,
    forma (dias, filas, columnas) -- igual que el resto de M3N).
    """
    from fase1_geografia import TIERRA, AGUA, correccion_altitud, irradiancia_absorbida_desde_toa_rejilla_con_albedo, estimar_T_inicial_equilibrio

    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])

    filas, columnas = tipo_superficie.shape

    capacidades, conductancias = construir_columna_rejilla(tipo_superficie, inercia_por_tipo, K_difusividad, n_capas)

    T_inicial = estimar_T_inicial_equilibrio(datos_orbita, albedo_grid, emisividad, profundidad_optica)
    T_columna = np.zeros((filas, columnas, n_capas))
    T_columna[..., 0] = T_inicial
    for i in range(1, n_capas):
        T_columna[..., i] = T_inicial  # arranque uniforme en toda la columna

    for ano in range(max_anos):
        T_inicio_ano = T_columna.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            T_columna = paso_multicapa(T_columna, abs_local, emisividad, capacidades, conductancias, paso_tiempo)

        diferencia_maxima = np.max(np.abs(T_columna[..., 0] - T_inicio_ano[..., 0]))
        if diferencia_maxima < tolerancia_convergencia:
            anos_convergencia = ano + 1
            break
    else:
        anos_convergencia = max_anos

    pasos_por_dia = round(ROTACION_PERIODO / paso_tiempo)

    registro_minima, registro_media, registro_maxima = [], [], []
    T_min_dia = T_max_dia = suma_dia = None
    contador_dia = 0

    def guardar_dia():
        T_media_dia = suma_dia / contador_dia
        registro_minima.append(correccion_altitud(T_min_dia - 273.15, altitud_metros))
        registro_media.append(correccion_altitud(T_media_dia - 273.15, altitud_metros))
        registro_maxima.append(correccion_altitud(T_max_dia - 273.15, altitud_metros))

    for paso, (toa, decl, ang_h_lon0) in enumerate(datos_orbita):
        if paso % pasos_por_dia == 0:
            if contador_dia > 0:
                guardar_dia()
            T_sup = T_columna[..., 0]
            T_min_dia = np.full_like(T_sup, np.inf)
            T_max_dia = np.full_like(T_sup, -np.inf)
            suma_dia = np.zeros_like(T_sup)
            contador_dia = 0

        abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
            toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
        )
        T_columna = paso_multicapa(T_columna, abs_local, emisividad, capacidades, conductancias, paso_tiempo)
        T_sup = T_columna[..., 0]

        T_min_dia = np.minimum(T_min_dia, T_sup)
        T_max_dia = np.maximum(T_max_dia, T_sup)
        suma_dia = suma_dia + T_sup
        contador_dia += 1

    if contador_dia == pasos_por_dia:
        guardar_dia()

    registro_minima = np.array(registro_minima)
    registro_media = np.array(registro_media)
    registro_maxima = np.array(registro_maxima)

    T_final_corregida = correccion_altitud(T_columna[..., 0] - 273.15, altitud_metros)

    return T_final_corregida, anos_convergencia, registro_minima, registro_media, registro_maxima


if __name__ == "__main__":
    from parametros import EMISIVIDAD, PROFUNDIDAD_OPTICA, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR
    from temperatura import simular as simular_una_capa

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    latitud_prueba_rad = math.radians(45)
    albedo_prueba = 0.30
    inercia_prueba = 2500  # INERCIA_POR_TIPO[TIERRA]

    print(f"\nDifusividad K = {K_DIFUSIVIDAD_TIERRA:.2e} m2/s")
    d_dia = profundidad_amortiguamiento(K_DIFUSIVIDAD_TIERRA, ROTACION_PERIODO)
    d_ano = profundidad_amortiguamiento(K_DIFUSIVIDAD_TIERRA, ORBITA_PERIODO)
    k_cond, c_v = derivar_k_y_cv(inercia_prueba, K_DIFUSIVIDAD_TIERRA)
    print(f"d_diurno = {d_dia:.3f} m | d_estacional = {d_ano:.3f} m")
    print(f"k (conductividad) = {k_cond:.3f} W/m/K | c_v = {c_v:.3e} J/m3/K")
    print("(Rango tipico roca densa: k 2-3 W/m/K, c_v 2-2.6e6 J/m3/K)")

    capacidades, conductancias, limites, centros = construir_columna(inercia_prueba, K_DIFUSIVIDAD_TIERRA)
    print(f"\nLimites de capas (m): {np.round(limites, 3)}")
    print(f"Centros de capas (m): {np.round(centros, 3)}")
    print(f"Capacidades (J/m2/K): {np.round(capacidades, 1)}")
    print(f"Conductancias entre capas (W/m2/K): {np.round(conductancias, 3)}")

    # ---- Prueba 0: N=1 capa debe coincidir EXACTO con el modelo de una capa ----
    print("\n--- Prueba 0: columna de 1 capa debe coincidir EXACTA con el modelo de una capa ---")
    T_final_1capa, registro_1capa = simular_una_capa(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA
    )
    # Con n_capas=1, la "columna" es identica al modelo de 1 capa: misma
    # capacidad total (toda ella en la unica capa), sin conductancias.
    T_final_multi, limites_1, centros_1, registro_multi = simular_multicapa_punto(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA,
        n_capas=1,
    )
    # OJO: con n_capas=1 la capacidad de esa unica capa depende de
    # FACTOR_PROFUNDIDAD_FONDO * c_v, que NO es igual a inercia*sqrt(dia/pi)
    # -- por diseño (no tendria sentido fisico que coincidieran). Este
    # caso limite no es comparable en valor absoluto al modelo de 1 capa;
    # se comprueba en su lugar en la Prueba 1 (conservacion de energia) y
    # Prueba 2 (comportamiento fisico del amortiguamiento). Se deja este
    # bloque solo para inspeccionar que la simulacion con n_capas=1 no
    # revienta y da un resultado razonable.
    print(f"Columna de 1 capa -- T final: {T_final_multi[0]:.2f} C (informativo, no comparable en valor absoluto)")

    # ---- Prueba 1: conservacion de energia de la conduccion en aislamiento ----
    print("\n--- Prueba 1: la conduccion vertical por si sola conserva energia y relaja al equilibrio ---")
    capacidades8, conductancias8, _, _ = construir_columna(inercia_prueba, K_DIFUSIVIDAD_TIERRA, n_capas=8)
    T_prueba = np.array([320.0, 300.0, 290.0, 280.0, 273.0, 270.0, 268.0, 265.0])  # perfil arbitrario, desequilibrado
    energia_inicial = np.sum(capacidades8 * T_prueba)
    # Paso implicito backward-Euler: incondicionalmente estable, asi que
    # para esta prueba (solo comprobar conservacion de energia y
    # relajacion al equilibrio, no la dinamica temporal real) se puede
    # usar un paso de tiempo mucho mayor que PASO_TIEMPO -- del orden de
    # la constante de tiempo mas lenta de la columna -- para llegar al
    # equilibrio en pocas iteraciones en vez de cientos de miles.
    tau_mas_lenta = capacidades8[-1] / conductancias8[-1]
    paso_tiempo_prueba = tau_mas_lenta / 2
    T_actual = T_prueba.copy()
    for _ in range(400):  # 200*tau_mas_lenta en total -- sobra para relajar del todo
        T_actual = paso_conduccion_implicito(T_actual, capacidades8, conductancias8, paso_tiempo_prueba)
    energia_final = np.sum(capacidades8 * T_actual)
    diferencia_energia_relativa = abs(energia_final - energia_inicial) / abs(energia_inicial)
    dispersion_final = T_actual.max() - T_actual.min()
    print(f"Perfil inicial (C): {np.round(T_prueba - 273.15, 1)}")
    print(f"Perfil final   (C): {np.round(T_actual - 273.15, 3)}")
    print(f"Dispersion final entre capas: {dispersion_final:.4f} K (deberia ser ~0, ya equilibradas)")
    print(f"Energia: diferencia relativa {diferencia_energia_relativa:.2e}")
    print("OK: conserva energia y relaja al equilibrio." if diferencia_energia_relativa < 1e-9 and dispersion_final < 0.01
          else "AVISO: revisar.")

    # ---- Prueba 2: la columna SI amortigua el ciclo estacional, y mas capas = mas fisico ----
    print("\n--- Prueba 2: amortiguamiento estacional en superficie -- una capa vs 8 capas ---")
    pasos_por_dia = round(ROTACION_PERIODO / PASO_TIEMPO)
    num_dias = len(registro_1capa) // pasos_por_dia

    def medias_diarias(registro_lista_arrays_o_escalares, indice_capa=None):
        registro = np.array(registro_lista_arrays_o_escalares)
        if indice_capa is not None:
            registro = registro[:, indice_capa]
        return np.array([
            np.mean(registro[d * pasos_por_dia:(d + 1) * pasos_por_dia])
            for d in range(num_dias)
        ])

    medias_1capa = medias_diarias(registro_1capa)

    T_final_multi8, limites8, centros8, registro_multi8 = simular_multicapa_punto(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA,
        n_capas=8,
    )
    medias_superficie_8capas = medias_diarias(registro_multi8, indice_capa=0)

    amplitud_1capa = medias_1capa.max() - medias_1capa.min()
    amplitud_8capas = medias_superficie_8capas.max() - medias_superficie_8capas.min()
    print(f"Amplitud estacional -- una capa: {amplitud_1capa:.2f} C")
    print(f"Amplitud estacional -- superficie de columna de 8 capas: {amplitud_8capas:.2f} C")
    print(f"Reduccion: {100 * (1 - amplitud_8capas / amplitud_1capa):.1f}%")

    registro_multi8_arr = np.array(registro_multi8)  # (pasos, n_capas)
    print("\nRango anual por capa (superficie -> profundidad):")
    for i in range(8):
        rango_i = registro_multi8_arr[:, i].max() - registro_multi8_arr[:, i].min()
        print(f"  Capa {i} (centro {centros8[i]:.2f} m): {registro_multi8_arr[:, i].min():.1f} C a "
              f"{registro_multi8_arr[:, i].max():.1f} C (rango {rango_i:.1f} C)")
    print("OK: el amortiguamiento aumenta monotonamente con la profundidad."
          if all(
              (registro_multi8_arr[:, i].max() - registro_multi8_arr[:, i].min())
              >= (registro_multi8_arr[:, i + 1].max() - registro_multi8_arr[:, i + 1].min()) - 1e-6
              for i in range(7)
          )
          else "AVISO: se esperaba que el rango anual decreciera con la profundidad -- revisar.")

    # ---- Prueba 3: convergencia con el numero de capas (N=4 vs N=8 vs N=12) ----
    print("\n--- Prueba 3: convergencia numerica -- la temperatura de superficie no deberia depender mucho de N ---")
    resultados_n = {}
    for n in (4, 8, 12):
        _, _, _, registro_n = simular_multicapa_punto(
            latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA,
            n_capas=n,
        )
        medias_n = medias_diarias(registro_n, indice_capa=0)
        resultados_n[n] = medias_n
        print(f"N={n:>2} capas -- amplitud estacional superficie: {medias_n.max() - medias_n.min():.3f} C")

    diferencia_4_8 = np.max(np.abs(resultados_n[4] - resultados_n[8]))
    diferencia_8_12 = np.max(np.abs(resultados_n[8] - resultados_n[12]))
    print(f"Diferencia maxima diaria N=4 vs N=8: {diferencia_4_8:.3f} C")
    print(f"Diferencia maxima diaria N=8 vs N=12: {diferencia_8_12:.3f} C")
    print("OK: converge (la diferencia se reduce al aumentar N)." if diferencia_8_12 < diferencia_4_8
          else "AVISO: revisar -- se esperaba que aumentar N acercara mas el resultado.")

    # ---- Prueba 4: rejilla completa -- agua sin cambios frente a Fase 1, tierra con columna ----
    print("\n--- Prueba 4: rejilla completa (mapa mixto) -- agua identica a Fase 1, tierra con columna de N capas ---")
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA, simular_rejilla_geografia_con_registro
    from fase1_mapa_mixto import mapa_mitad_agua_mitad_tierra

    tipo_mixto, altitud_mixta = mapa_mitad_agua_mitad_tierra()

    T_f1, anos_f1, reg_min_f1, reg_media_f1, reg_max_f1 = simular_rejilla_geografia_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )
    T_multi, anos_multi, reg_min_multi, reg_media_multi, reg_max_multi = simular_rejilla_multicapa_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
        n_capas=8,
    )

    mascara_agua = tipo_mixto == AGUA
    mascara_tierra = tipo_mixto == TIERRA

    diferencia_agua_media = np.max(np.abs(reg_media_multi[:, mascara_agua] - reg_media_f1[:, mascara_agua]))
    print(f"Agua -- diferencia maxima en registro medio diario vs Fase 1: {diferencia_agua_media:.2e} C")
    print("OK: coincide EXACTA en agua." if diferencia_agua_media < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    amplitud_tierra_f1 = reg_media_f1[:, mascara_tierra].max(axis=0) - reg_media_f1[:, mascara_tierra].min(axis=0)
    amplitud_tierra_multi = reg_media_multi[:, mascara_tierra].max(axis=0) - reg_media_multi[:, mascara_tierra].min(axis=0)
    print(f"Tierra -- amplitud estacional media (Fase 1, sin capa profunda): "
          f"{amplitud_tierra_f1.mean():.2f} C (promedio de celdas de tierra)")
    print(f"Tierra -- amplitud estacional media (columna de 8 capas): "
          f"{amplitud_tierra_multi.mean():.2f} C (promedio de celdas de tierra)")
    print(f"Convergencia -- Fase 1: {anos_f1} año(s) | Columna multicapa: {anos_multi} año(s)")
    print("OK: la tierra se amortigua con la columna de N capas, el agua no cambia."
          if diferencia_agua_media < 1e-6 and amplitud_tierra_multi.mean() < amplitud_tierra_f1.mean()
          else "AVISO: revisar.")

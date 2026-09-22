# fase2_difusion.py -- Fase 2: transporte de calor horizontal (difusion).
#
# Añade difusion lateral de calor entre celdas vecinas a la rejilla de
# Fase 1 (fase1_geografia.py).
#
# ================================================================
# CORRECCION IMPORTANTE (encontrada DESPUES de construir la primera
# version de este archivo, al verificar que la difusion realmente
# suavizaba algo): la primera version incluia un factor 1/R^2 (R =
# radio del planeta) en el operador, como si D fuera una difusividad
# fisica "de verdad" (conductividad) aplicada sobre distancias reales
# en metros. Con ese factor, el efecto de la difusion resultaba
# absurdamente pequeño (cambios de ~1e-7 grados en vez de varios
# grados) porque R^2 en P3N es del orden de 2*10^13 m^2 -- divide
# el termino hasta hacerlo irrelevante. Lo detecte precisamente
# porque comprobe el efecto real en la frontera agua/tierra del mapa
# de prueba y no habia cambio ninguno, cosa que no cuadraba.
#
# La causa: el valor de referencia de la literatura, D ~ 0.55 W/m2/K
# (Budyko-Sellers / North 1975), NO esta definido para un operador en
# distancias fisicas (metros). Esta definido para el operador escrito
# directamente en LATITUD EN RADIANES (una coordenada angular,
# adimensional), sin ningun R^2 de por medio:
#
#     C dT/dt = ... + (D/cos phi) * d/dphi( cos phi * dT/dphi )
#
# (esta es la forma 1D zonal estandar; la generalizacion a 2D con
# longitud añade el termino equivalente en lambda). El R^2 real del
# planeta esta IMPLICITO en el valor calibrado de D -- no se divide
# aparte. Por eso este archivo ahora usa el operador en radianes, sin
# radio_planeta, que es la forma correcta de usar D ~ 0.55 W/m2/K:
#
#     D * [ (1/cos^2 phi) d^2T/dlambda^2  +  (1/cos phi) d/dphi( cos phi dT/dphi ) ]
#
# phi = latitud (rad), lambda = longitud (rad). D es el "coeficiente
# de difusion" (W/m^2/K) -- en el clima real NO es conduccion literal
# por el terreno: es una forma simplificada de representar el
# transporte de calor por la atmosfera (vientos) y las corrientes
# oceanicas, que es lo que de verdad iguala temperaturas entre
# lugares distintos. El valor de referencia habitual en la literatura
# es D ~ 0.55 W/m2/K (un unico valor global, no uno por tierra/agua
# -- ver la nota en DISENO_FASE2.md sobre por que).
#
# CORRECCION 2 (encontrada al comparar el mapa de calor resultante con
# el de antes de la difusion -- Carlos noto que el contraste
# continente/oceano habia desaparecido casi del todo, mas de lo
# esperable): el valor D=0.55 de la literatura esta calibrado con
# datos de la Tierra real (radio ~6.371.000 m). Al usar el operador en
# radianes SIN R^2 (ver arriba), el radio del planeta deja de aparecer
# explicitamente en la formula -- pero el valor D=0.55 sigue siendo,
# implicitamente, "el D que le corresponde a un planeta del tamaño de
# la Tierra". P3N es mas pequeño (radio ~4.677.000 m, un 73% del
# terrestre): para un mismo D en radianes, un grado de longitud cubre
# menos distancia fisica real en un planeta pequeño, asi que el calor
# tiene menos terreno que recorrer para igualar zonas -- el mismo
# numero produce una difusion relativamente mas fuerte en P3N que en
# la Tierra. Eso es justo lo que se notaba en el mapa.
#
# La correccion: reescalar por el cociente de radios al cuadrado,
# D_P3N = D_tierra * (R_P3N/R_tierra)^2 ~= 0.55 * 0.54 ~= 0.30 W/m2/K
# (ver el calculo explicito mas abajo, D_DIFUSION_REFERENCIA). Esto
# asume que el mecanismo fisico real que D representa (transporte por
# vientos y corrientes) tiene una eficiencia intrinseca comparable en
# ambos planetas, y que la diferencia esta solo en el tamaño -- una
# aproximacion razonable de partida (es mejor que reutilizar el numero
# de la Tierra sin mas), pero no una certeza: la dinamica atmosferica
# real de P3N podria ser distinta a la terrestre por otras razones
# (rotacion, insolacion, composicion atmosferica...) que este
# reescalado no captura. Lo señalo como limitacion, igual que antes.
# ================================================================
#
# SOBRE LOS POLOS: el termino meridional se evalua en la FRONTERA entre
# filas (no en el centro), y la frontera norte de la fila 0 es
# exactamente el polo (latitud 90 grados). cos(90 grados) = 0, asi que
# el flujo que cruza esa frontera es automaticamente CERO en esta forma
# -- es la condicion de contorno fisicamente correcta (un punto no
# tiene circunferencia, no puede cruzar nada "a traves" de el en este
# sentido). Por eso este archivo NO reutiliza el envolvido polar de
# C3N: no hace falta aqui, y forzarlo estropearia la conservacion de
# energia exacta que da esta forma. El envolvido de C3N es para otra
# cosa: identificar que celda del mapa corresponde a un pixel que se
# salio del lienzo, que es un problema distinto.
#
# Lo que SI iguala rapidamente las celdas cercanas al polo entre si
# (incluidas las de longitudes opuestas) es el propio termino zonal:
# su factor 1/cos^2(phi) crece mucho cerca de los polos (las celdas
# son geometricamente muy pequeñas en esa direccion alli), asi que la
# difusion este-oeste ya iguala casi al instante todas las celdas de
# una fila cercana al polo -- ese es el mecanismo real por el que un
# "punto caliente" cerca del polo se reparte enseguida por todo el
# entorno del polo, sin necesitar un atajo especial.
#
# SOBRE EL METODO NUMERICO: la primera version de este archivo resolvia
# la difusion de forma EXPLICITA, subdividiendo cada paso_tiempo en
# muchos sub-pasos pequeños para mantener la estabilidad numerica. Esto
# resulto ser computacionalmente inviable: el "problema del polo" (las
# celdas cercanas al polo tienen un espaciado este-oeste real minusculo,
# lo que fuerza un paso de tiempo estable extremadamente pequeño para
# TODA la rejilla bajo un esquema explicito) obligaba a cientos de
# sub-pasos por cada paso_tiempo principal, y una simulacion de un año
# no terminaba ni en varios minutos.
#
# La solucion estandar y correcta para este problema (comun en todos
# los modelos climaticos reales en rejilla lat/lon) es resolver la
# difusion de forma IMPLICITA (metodo de Euler implicito / backward
# Euler), que es incondicionalmente estable: no importa lo grande que
# sea el paso de tiempo, la solucion no diverge numericamente (la
# precision fisica de usar un paso grande es una cuestion aparte de la
# estabilidad, y para encontrar un equilibrio climatico a largo plazo,
# que es lo que hace este modelo, es la eleccion correcta).
#
# Concretamente se usa "operator splitting" (separacion de operadores)
# en cada paso de tiempo:
#   1. Paso radiativo, EXPLICITO (igual que en Fase 0/1):
#        T* = T + (paso_tiempo / C) * (absorbido - emitido)
#   2. Paso difusivo, IMPLICITO:
#        C * (T_nuevo - T*) / paso_tiempo = divergencia_difusion(T_nuevo)
#      Esto es un sistema lineal en T_nuevo (porque divergencia_difusion
#      es lineal en T): (diag(C/paso_tiempo) - L) @ T_nuevo = (C/paso_tiempo) * T*
#      donde L es la matriz dispersa (sparse) que representa el operador
#      de difusion (construida una sola vez, ya que D_grid y C_grid no
#      cambian durante la simulacion). La matriz A = diag(C/paso_tiempo) - L
#      se factoriza UNA SOLA VEZ (scipy.sparse.linalg.splu) al principio
#      de la simulacion, y esa factorizacion se reutiliza para resolver
#      el sistema en cada paso de tiempo, lo que hace cada paso muy
#      rapido (una sustitucion hacia adelante/atras, no una factorizacion
#      nueva cada vez).

import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from parametros import ROTACION_PERIODO, PI, CONSTANTE_SB, P3N_RADIO
from temperatura import precalcular_orbita, t_eq, PASO_TIEMPO as PASO_TIEMPO_POR_DEFECTO
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS, FILAS, COLUMNAS, GRADOS_POR_FILA, GRADOS_POR_COL
from fase1_geografia import (
    TIERRA, AGUA, correccion_altitud, estimar_T_inicial_equilibrio,
    irradiancia_absorbida_desde_toa_rejilla_con_albedo,
)

# ================================================================
# D_DIFUSION_REFERENCIA: de donde sale y por que esta escrito asi.
#
# El valor de la literatura (Budyko-Sellers / North 1975), D ~ 0.55
# W/m2/K, esta calibrado para la Tierra real -- su radio esta metido
# implicitamente en ese numero (ver la nota larga mas abajo sobre por
# que la formula no lleva R^2 explicito). Usar ese 0.55 tal cual en
# P3N seria repetir, de otra forma, el mismo error que ya se cometio
# una vez con este mismo archivo (dividir sin querer por un R^2 que
# no correspondia): tratar una cifra calibrada para la Tierra como si
# fuera valida para cualquier planeta.
#
# Por eso el valor de la Tierra y el radio de referencia se dejan
# EXPLICITOS aqui, y P3N_RADIO_REFERENCIA se calcula a partir de
# ellos con la formula de reescalado (D_planeta = D_tierra *
# (R_planeta/R_tierra)^2 -- ver la nota de correccion mas abajo para
# la derivacion). Si algun dia cambia el radio de P3N en parametros.py,
# o si se reutiliza esta misma logica para otro planeta ficticio, el
# valor correcto sale solo con volver a ejecutar este archivo -- no
# hay un numero ya calculado escondido que se pueda quedar
# desactualizado sin que nadie se de cuenta.
D_DIFUSION_REFERENCIA_TIERRA = 0.55  # W/m2/K -- valor de la literatura, tal cual, SIN adaptar a P3N
RADIO_TIERRA_REFERENCIA = 6_371_000  # m -- el radio para el que esta calibrado ese 0.55
D_DIFUSION_REFERENCIA = D_DIFUSION_REFERENCIA_TIERRA * (P3N_RADIO / RADIO_TIERRA_REFERENCIA) ** 2

DLAMBDA_RAD = math.radians(GRADOS_POR_COL)
DPHI_RAD = math.radians(GRADOS_POR_FILA)
_PHI_GRID = np.radians(LATITUDES_GRADOS).reshape(-1, 1)  # (FILAS, 1), constante para toda la rejilla


def _vecino_norte(T):
    """Vecino hacia el polo norte (fila - 1) de cada celda. La fila 0
    no tiene una fila -1 real; se repite la fila 0 (el valor no importa
    para el resultado, porque su flujo se multiplica por cos(90 deg) =
    0 -- ver cabecera del archivo), asi que basta con no salirse del
    array."""
    vecino = np.empty_like(T)
    vecino[1:, :] = T[:-1, :]
    vecino[0, :] = T[0, :]
    return vecino


def _vecino_sur(T):
    """Igual que _vecino_norte pero hacia el polo sur (fila + 1)."""
    vecino = np.empty_like(T)
    vecino[:-1, :] = T[1:, :]
    vecino[-1, :] = T[-1, :]
    return vecino


def _vecino_este(T):
    return np.roll(T, -1, axis=1)


def _vecino_oeste(T):
    return np.roll(T, 1, axis=1)


def divergencia_difusion(T, D_grid):
    """
    Operador de difusion esferico de T por diferencias finitas (forma
    de flujo, conservativa -- ver cabecera del archivo), escrito
    directamente en radianes de latitud/longitud (SIN dividir por el
    radio del planeta -- ver la nota de correccion en la cabecera del
    archivo sobre por que esta es la forma que corresponde al valor de
    D de la literatura). Devuelve el termino en W/m^2, listo para
    sumarlo a la ecuacion de balance de energia:
        C * dT/dt = absorbido - emitido + divergencia_difusion(...)

    D_grid: array (FILAS, COLUMNAS) con el coeficiente de difusion de
    cada celda (W/m2/K). Puede ser uniforme (mismo valor en todas las
    celdas) o distinto por tipo de celda.

    Esta funcion se mantiene (ya no la usa el bucle principal, que
    ahora usa la matriz dispersa equivalente -- ver
    construir_matriz_difusion()) porque es la definicion de referencia,
    clara y directa, contra la que se verifica que la matriz dispersa
    es exactamente equivalente (prueba 0 mas abajo).
    """
    # ---- termino zonal (este-oeste) ----
    T_este = _vecino_este(T)
    T_oeste = _vecino_oeste(T)
    termino_zonal = (T_este - 2 * T + T_oeste) / DLAMBDA_RAD**2 / np.cos(_PHI_GRID)**2

    # ---- termino meridional (norte-sur), forma de flujo ----
    T_norte = _vecino_norte(T)
    T_sur = _vecino_sur(T)

    phi_frontera_norte = _PHI_GRID + DPHI_RAD / 2  # frontera con la fila norte (latitud mas alta)
    phi_frontera_sur = _PHI_GRID - DPHI_RAD / 2    # frontera con la fila sur (latitud mas baja)

    flujo_norte = np.cos(phi_frontera_norte) * (T_norte - T) / DPHI_RAD
    flujo_sur = np.cos(phi_frontera_sur) * (T - T_sur) / DPHI_RAD

    termino_meridional = (flujo_norte - flujo_sur) / DPHI_RAD / np.cos(_PHI_GRID)

    return D_grid * (termino_zonal + termino_meridional)


def construir_matriz_difusion(D_grid):
    """
    Construye la matriz dispersa (sparse) L de tamaño (FILAS*COLUMNAS,
    FILAS*COLUMNAS) tal que, para cualquier campo de temperatura T:

        (L @ T.flatten()).reshape(FILAS, COLUMNAS) == divergencia_difusion(T, D_grid)

    (verificado numericamente en la prueba 0 del bloque __main__). Es
    la misma formula de diferencias finitas, escrita como sistema
    lineal en vez de aplicada directamente sobre el array, para poder
    resolverla de forma implicita.

    Indexado de la rejilla aplanada: idx(fila, col) = fila*COLUMNAS + col.
    """
    N = FILAS * COLUMNAS

    def idx(fila, col):
        return fila * COLUMNAS + col

    filas_idx = []
    cols_idx = []
    valores = []

    cos_phi = np.cos(_PHI_GRID).flatten()  # (FILAS,)
    cos_frontera_norte = np.cos(_PHI_GRID.flatten() + DPHI_RAD / 2)  # (FILAS,)
    cos_frontera_sur = np.cos(_PHI_GRID.flatten() - DPHI_RAD / 2)   # (FILAS,)

    for i in range(FILAS):
        coef_zonal = 1.0 / (DLAMBDA_RAD**2 * cos_phi[i]**2)
        coef_norte = cos_frontera_norte[i] / (DPHI_RAD**2 * cos_phi[i])
        coef_sur = cos_frontera_sur[i] / (DPHI_RAD**2 * cos_phi[i])

        for j in range(COLUMNAS):
            d = D_grid[i, j]
            fila_out = idx(i, j)

            # ---- termino zonal: T_este - 2T + T_oeste (periodico en columnas) ----
            j_este = (j + 1) % COLUMNAS
            j_oeste = (j - 1) % COLUMNAS
            filas_idx.append(fila_out); cols_idx.append(idx(i, j_este)); valores.append(d * coef_zonal)
            filas_idx.append(fila_out); cols_idx.append(idx(i, j_oeste)); valores.append(d * coef_zonal)
            filas_idx.append(fila_out); cols_idx.append(fila_out); valores.append(d * (-2.0 * coef_zonal))

            # ---- termino meridional: sin envolvido (frontera con flujo -> 0 en los polos) ----
            if i > 0:
                filas_idx.append(fila_out); cols_idx.append(idx(i - 1, j)); valores.append(d * coef_norte)
            if i < FILAS - 1:
                filas_idx.append(fila_out); cols_idx.append(idx(i + 1, j)); valores.append(d * coef_sur)
            filas_idx.append(fila_out); cols_idx.append(fila_out); valores.append(d * (-(coef_norte + coef_sur)))

    L = sp.coo_matrix((valores, (filas_idx, cols_idx)), shape=(N, N)).tocsr()
    return L


def simular_rejilla_difusion(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.01,
):
    """
    Igual que simular_rejilla_geografia() (fase1_geografia.py), con un
    termino de difusion horizontal añadido. Si D_grid es todo ceros,
    debe coincidir EXACTAMENTE con Fase 1 (validacion obligatoria).

    Metodo: "operator splitting" por paso de tiempo -- paso radiativo
    explicito (igual que siempre) + paso difusivo implicito (matriz
    dispersa factorizada UNA VEZ al principio, reutilizada en cada
    paso). Ver la nota larga en la cabecera del archivo.
    """
    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])
    C_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)
    C_flat = C_grid.flatten()

    N = FILAS * COLUMNAS
    tiene_difusion = np.max(D_grid) > 0

    if tiene_difusion:
        L = construir_matriz_difusion(D_grid)
        A = sp.diags(C_flat / paso_tiempo) - L
        A_factorizada = splu(A.tocsc())
    else:
        A_factorizada = None

    T = estimar_T_inicial_equilibrio(datos_orbita, albedo_grid, emisividad, profundidad_optica)

    for ano in range(max_anos):
        T_inicio_ano = T.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            emitido = (1 - emisividad / 2) * CONSTANTE_SB * T**4

            # ---- paso 1: radiativo, explicito ----
            T_estrella = T + (paso_tiempo / C_grid) * (abs_local - emitido)

            # ---- paso 2: difusivo, implicito (si D_grid tiene difusion) ----
            if tiene_difusion:
                lado_derecho = (C_flat / paso_tiempo) * T_estrella.flatten()
                T = A_factorizada.solve(lado_derecho).reshape(FILAS, COLUMNAS)
            else:
                T = T_estrella

        diferencia_maxima = np.max(np.abs(T - T_inicio_ano))
        if diferencia_maxima < tolerancia_convergencia:
            return correccion_altitud(T - 273.15, altitud_metros), ano + 1

    return correccion_altitud(T - 273.15, altitud_metros), max_anos


def simular_rejilla_difusion_con_registro(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.01,
):
    """
    Igual que simular_rejilla_geografia_con_registro() (fase1_geografia.py),
    con difusion horizontal añadida -- mismo metodo que
    simular_rejilla_difusion() (paso radiativo explicito + paso difusivo
    implicito, matriz factorizada una vez). Primero converge al
    equilibrio (igual que simular_rejilla_difusion), y despues hace UN
    año mas registrando el minimo/medio/maximo diario de cada celda --
    exactamente igual que la version de Fase 1, para que
    consulta_punto.py, mapa_calor.py y analisis_latitudes.py puedan usar
    esta version sin cambiar su forma de leer el resultado (mismo
    orden y forma de los 5 valores que devuelve).
    """
    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])
    C_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)
    C_flat = C_grid.flatten()

    tiene_difusion = np.max(D_grid) > 0
    if tiene_difusion:
        L = construir_matriz_difusion(D_grid)
        A = sp.diags(C_flat / paso_tiempo) - L
        A_factorizada = splu(A.tocsc())
    else:
        A_factorizada = None

    def paso_fisica(T, abs_local):
        emitido = (1 - emisividad / 2) * CONSTANTE_SB * T**4
        T_estrella = T + (paso_tiempo / C_grid) * (abs_local - emitido)
        if tiene_difusion:
            lado_derecho = (C_flat / paso_tiempo) * T_estrella.flatten()
            return A_factorizada.solve(lado_derecho).reshape(FILAS, COLUMNAS)
        return T_estrella

    T = estimar_T_inicial_equilibrio(datos_orbita, albedo_grid, emisividad, profundidad_optica)

    for ano in range(max_anos):
        T_inicio_ano = T.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            T = paso_fisica(T, abs_local)
        diferencia_maxima = np.max(np.abs(T - T_inicio_ano))
        if diferencia_maxima < tolerancia_convergencia:
            anos_convergencia = ano + 1
            break
    else:
        anos_convergencia = max_anos

    pasos_por_dia = round(ROTACION_PERIODO / paso_tiempo)

    registro_minima = []
    registro_media = []
    registro_maxima = []

    T_min_dia = None
    T_max_dia = None
    suma_dia = None
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
            T_min_dia = np.full_like(T, np.inf)
            T_max_dia = np.full_like(T, -np.inf)
            suma_dia = np.zeros_like(T)
            contador_dia = 0

        abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
            toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
        )
        T = paso_fisica(T, abs_local)

        T_min_dia = np.minimum(T_min_dia, T)
        T_max_dia = np.maximum(T_max_dia, T)
        suma_dia = suma_dia + T
        contador_dia += 1

    if contador_dia == pasos_por_dia:
        guardar_dia()
    elif contador_dia > 0:
        print(f"Aviso: se descarta el ultimo dia registrado por estar incompleto "
              f"({contador_dia} de {pasos_por_dia} pasos).")

    registro_minima = np.array(registro_minima)
    registro_media = np.array(registro_media)
    registro_maxima = np.array(registro_maxima)

    T_final_corregida = correccion_altitud(T - 273.15, altitud_metros)

    return T_final_corregida, anos_convergencia, registro_minima, registro_media, registro_maxima


if __name__ == "__main__":
    import time
    from parametros import EMISIVIDAD, PROFUNDIDAD_OPTICA, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, mapa_falso_todo_tierra, simular_rejilla_geografia

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo_superficie, altitud_metros = mapa_falso_todo_tierra()

    # ---- Prueba 0: la matriz dispersa reproduce EXACTO divergencia_difusion() ----
    print("\n--- Prueba 0: la matriz dispersa equivale a divergencia_difusion() ---")
    rng = np.random.default_rng(42)
    D_prueba = rng.uniform(0.1, 1.0, size=(FILAS, COLUMNAS))
    L_prueba = construir_matriz_difusion(D_prueba)
    peor_diferencia_relativa = 0.0
    for _ in range(5):
        T_prueba = 250 + rng.normal(0, 30, size=(FILAS, COLUMNAS))
        directo = divergencia_difusion(T_prueba, D_prueba)
        via_matriz = (L_prueba @ T_prueba.flatten()).reshape(FILAS, COLUMNAS)
        escala = np.max(np.abs(directo))
        peor_diferencia_relativa = max(peor_diferencia_relativa, np.max(np.abs(directo - via_matriz)) / escala)
    print(f"Diferencia relativa maxima (5 campos T aleatorios): {peor_diferencia_relativa:.2e}")
    print("OK: la matriz dispersa es exactamente equivalente." if peor_diferencia_relativa < 1e-10 else "AVISO: no coinciden -- revisar.")

    # ---- Prueba 1: D = 0 debe coincidir EXACTO con Fase 1 ----
    print("\n--- Prueba 1: difusion desactivada (D=0) debe coincidir con Fase 1 ---")
    D_cero = np.zeros((FILAS, COLUMNAS))
    T_fase1, anos_fase1 = simular_rejilla_geografia(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )
    T_fase2_sin_dif, anos_fase2_sin_dif = simular_rejilla_difusion(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
        PROFUNDIDAD_OPTICA, D_cero,
    )
    diferencia = np.max(np.abs(T_fase1 - T_fase2_sin_dif))
    print(f"Convergencia: {anos_fase1} (Fase 1) vs {anos_fase2_sin_dif} (Fase 2, D=0)")
    print(f"Diferencia maxima: {diferencia:.2e} C")
    print("OK: coincide EXACTA." if diferencia < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    # ---- Prueba 2: conservacion de energia global, D UNIFORME (caso de referencia) ----
    # OJO: la conservacion exacta de esta formula (D * Laplaciano(T), con D
    # sacado fuera del operador) solo esta garantizada cuando D es uniforme
    # en toda la rejilla. Si D varia por celda (p.ej. un valor distinto para
    # tierra y otro para agua), esta forma deja de ser estrictamente
    # conservativa -- para eso haria falta evaluar D en la FRONTERA entre
    # celdas (forma "D en la cara"), no en el centro. Ver el aviso que
    # imprime este mismo bloque para el caso D variable.
    print("\n--- Prueba 2: conservacion de energia global de divergencia_difusion() ---")
    peso = np.cos(_PHI_GRID)  # ponderacion por area real de la celda (mas peso cerca del ecuador)
    T_aleatoria = 250 + rng.normal(0, 40, size=(FILAS, COLUMNAS))

    D_uniforme_prueba = np.full((FILAS, COLUMNAS), D_DIFUSION_REFERENCIA)
    dif_uniforme = divergencia_difusion(T_aleatoria, D_uniforme_prueba)
    energia_neta_uniforme = np.sum(dif_uniforme * peso)
    escala_uniforme = np.sum(np.abs(dif_uniforme) * peso)
    print(f"D uniforme -- energia neta (deberia ser ~0): {energia_neta_uniforme:.3e}  (escala: {escala_uniforme:.3e})")
    print("OK: conserva energia con D uniforme." if abs(energia_neta_uniforme) < 1e-6 * escala_uniforme else "AVISO: no conserva -- revisar.")

    dif_variable = divergencia_difusion(T_aleatoria, D_prueba)  # D_prueba: aleatorio por celda, de la Prueba 0
    energia_neta_variable = np.sum(dif_variable * peso)
    escala_variable = np.sum(np.abs(dif_variable) * peso)
    print(f"D variable por celda -- energia neta: {energia_neta_variable:.3e}  (escala: {escala_variable:.3e}, "
          f"{100*abs(energia_neta_variable)/escala_variable:.1f}% relativo)")
    print("AVISO ESPERADO: con D variable por celda esta forma NO conserva energia exactamente "
          "(ver nota en el codigo) -- pendiente de decidir con Carlos si hace falta D por tipo.")

    # ---- Prueba 3: mapa mixto agua/tierra real, con difusion activada, D uniforme de referencia ----
    print("\n--- Prueba 3: mapa mixto tierra/agua (valores reales de Fase 1), D uniforme de referencia, cronometrado ---")
    from fase1_mapa_mixto import mapa_mitad_agua_mitad_tierra
    tipo_mixto, altitud_mixta = mapa_mitad_agua_mitad_tierra()

    print("Referencia: mismo mapa, MISMA fisica, SIN difusion (Fase 1 pura) --")
    inicio = time.time()
    T_sin_dif, anos_sin_dif = simular_rejilla_geografia(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )
    duracion_sin_dif = time.time() - inicio
    print(f"  Fase 1 (sin difusion): {anos_sin_dif} año(s) en {duracion_sin_dif:.2f} s")

    D_referencia = np.full((FILAS, COLUMNAS), D_DIFUSION_REFERENCIA)
    inicio = time.time()
    T_mixto, anos_mixto = simular_rejilla_difusion(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
        PROFUNDIDAD_OPTICA, D_referencia,
    )
    duracion = time.time() - inicio
    print(f"  Fase 2 (con difusion implicita, D={D_DIFUSION_REFERENCIA}): {anos_mixto} año(s) en {duracion:.2f} s")
    print(f"  Rango de temperatura resultante (un instante congelado del ciclo, NO el extremo anual real): "
          f"{T_mixto.min():.1f} C a {T_mixto.max():.1f} C "
          f"(sin difusion: {T_sin_dif.min():.1f} C a {T_sin_dif.max():.1f} C)")
    print("OK: el metodo implicito termina sin problema (el explicito con sub-pasos no terminaba ni en 5 minutos).")

    # ---- Prueba 4: version con registro diario (la que usan consulta_punto.py,
    #                mapa_calor.py, analisis_latitudes.py) -- verificar que D=0
    #                coincide exacto con Fase 1, y ver el rango anual REAL
    #                (minimo y maximo de verdad, con ciclo diurno y estacional
    #                incluido, no solo un instante congelado como la Prueba 3) ----
    print("\n--- Prueba 4: version con registro diario (la que usan las herramientas) ---")
    from fase1_geografia import simular_rejilla_geografia_con_registro

    print("D=0 -- debe coincidir EXACTO con la version de registro de Fase 1 --")
    _, anos_f1, reg_min_f1, reg_media_f1, reg_max_f1 = simular_rejilla_geografia_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )
    _, anos_f2, reg_min_f2, reg_media_f2, reg_max_f2 = simular_rejilla_difusion_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
        PROFUNDIDAD_OPTICA, D_cero,
    )
    diferencia_registro = max(
        np.max(np.abs(reg_min_f1 - reg_min_f2)),
        np.max(np.abs(reg_media_f1 - reg_media_f2)),
        np.max(np.abs(reg_max_f1 - reg_max_f2)),
    )
    print(f"Convergencia: {anos_f1} (Fase 1) vs {anos_f2} (Fase 2, D=0) | diferencia maxima en el registro: {diferencia_registro:.2e} C")
    print("OK: coincide EXACTA." if diferencia_registro < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    print("\nCon difusion (D=0.55) -- rango anual REAL (con ciclo diurno + estacional, todo el año) --")
    inicio = time.time()
    _, anos_f2d, reg_min_f2d, reg_media_f2d, reg_max_f2d = simular_rejilla_difusion_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
        PROFUNDIDAD_OPTICA, D_referencia,
    )
    duracion_f2d = time.time() - inicio
    print(f"  SIN difusion -- minimo real del año: {reg_min_f1.min():.1f} C | maximo real del año: {reg_max_f1.max():.1f} C")
    print(f"  CON difusion -- minimo real del año: {reg_min_f2d.min():.1f} C | maximo real del año: {reg_max_f2d.max():.1f} C")
    print(f"  ({anos_f2d} año(s) de convergencia, {duracion_f2d:.1f} s con el registro diario)")

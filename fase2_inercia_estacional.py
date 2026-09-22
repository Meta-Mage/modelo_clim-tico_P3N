# fase2_inercia_estacional.py -- Fase 2: inercia termica estacional
# (arquitectura de dos capas).
#
# Hasta ahora el modelo tiene UNA sola capa por celda, con una inercia
# termica ("inercia", el parametro en INERCIA_POR_TIPO) calibrada para
# la respuesta DIURNA (dia/noche): C = inercia * sqrt(ROTACION_PERIODO/PI).
# Eso significa que, a lo largo del año, la temperatura de cada celda
# sigue el ciclo de estaciones casi sin retraso ni amortiguacion -- como
# si el suelo (o el agua) no tuviera ninguna "memoria" mas alla de un
# dia. En la realidad esto no es asi: el subsuelo profundo, y sobre
# todo el oceano profundo, absorben y liberan calor a lo largo de MESES,
# amortiguando y retrasando el ciclo estacional en la superficie (por
# eso el mar suaviza el clima costero, y el mes mas caluroso del año no
# suele coincidir exactamente con el solsticio de verano).
#
# ================================================================
# LA IDEA: separar cada celda en DOS capas acopladas
#   - T_superficie: la de siempre, responde a la radiacion cada paso de
#     tiempo, con su inercia diurna de siempre.
#   - T_profunda: una segunda capa, con MUCHA mas inercia (calibrada al
#     periodo del año, no del dia), que no ve la radiacion directamente
#     -- solo intercambia calor con la capa superficial.
#
# EL PARAMETRO "inercia" ES UNA PROPIEDAD DEL MATERIAL, NO DEL PERIODO:
# revisando las unidades de la formula ya existente, C = inercia *
# sqrt(P/pi), se ve que "inercia" tiene las unidades exactas de la
# EFUSIVIDAD TERMICA real, e = sqrt(k * rho * c_p) (conductividad *
# densidad * calor especifico) -- una propiedad fija del material
# (tierra o agua), que NO depende de si el periodo P que le metamos en
# la formula es el dia o el año. Por eso la capacidad de la capa
# profunda se puede calcular reutilizando el MISMO valor de "inercia"
# que ya tienes tabulado, cambiando solo el periodo:
#
#     C_profunda = inercia * sqrt(ORBITA_PERIODO / PI)
#
# Con el año de P3N (~270 dias) frente al dia (1 dia), esto da una capa
# profunda con una capacidad termica ~16 veces mayor que la superficial
# (sqrt(270) ~ 16.4) -- del orden de magnitud correcto para una "reserva"
# que amortigua el ciclo estacional sin llegar a anularlo del todo.
#
# EL ACOPLAMIENTO ENTRE LAS DOS CAPAS (k, en W/m2/K) es la parte donde
# NO hay una formula unica evidente: derivarlo con rigor total
# requeriria resolver la ecuacion de conduccion 1D completa en
# profundidad (un perfil continuo de temperatura bajo cada celda, no
# solo dos numeros), que es un cambio de arquitectura mucho mayor. En
# su lugar, aqui se usa una aproximacion estandar en modelos de
# parametros concentrados (el equivalente termico a un circuito RC):
# elegir k de forma que la constante de tiempo del acoplamiento,
# tau = C_profunda / k, sea igual al periodo del año dividido por 2*pi
# (el criterio habitual de "media potencia" de un filtro de primer
# orden -- el punto en el que la señal de esa frecuencia queda
# atenuada de forma significativa pero no arbitraria):
#
#     k = C_profunda * (2*PI / ORBITA_PERIODO)
#
# Esto es una eleccion razonada, no una cifra inventada al azar, pero
# tampoco es EL valor fisicamente exacto -- es una aproximacion de
# ingenieria, igual que lo es el propio hecho de reducir un perfil de
# profundidad continuo a solo dos capas. Lo señalo explicitamente
# porque es una simplificacion real del modelo, no un dato verificado
# como los de la Fase 0/1.
# ================================================================
#
# METODO NUMERICO: igual que en fase2_difusion.py, el paso de tiempo se
# separa en dos ("operator splitting"):
#   1. Paso radiativo, explicito, SOLO sobre la capa superficial (la
#      capa profunda no ve la radiacion directamente):
#        T_superficie* = T_superficie + (paso_tiempo/C_superficie) * (absorbido - emitido)
#   2. Paso de acoplamiento, implicito: como cada celda solo intercambia
#      calor con su propia capa profunda (no con las celdas vecinas,
#      a diferencia de la difusion horizontal), esto es un sistema
#      lineal 2x2 LOCAL por celda -- no hace falta una matriz dispersa
#      ni factorizacion, un despeje algebraico directo (regla de
#      Cramer) basta, y es incondicionalmente estable igual que la
#      difusion implicita.

import math
import numpy as np
from parametros import ROTACION_PERIODO, ORBITA_PERIODO, PI, CONSTANTE_SB
from temperatura import precalcular_orbita, t_eq, estimar_T_inicial_equilibrio_punto, PASO_TIEMPO
from geometria import angulo_cenital
from atmosfera import masa_aire
from radiacion import i_inst, trans, i_atm, i_abs

FRECUENCIA_ANGULAR_ANUAL = 2 * PI / ORBITA_PERIODO


def capacidad_profunda(inercia):
    """C_profunda (J/m2/K) -- mismo 'inercia' (efusividad del material),
    calculado con el periodo del año en vez del periodo del dia."""
    return inercia * math.sqrt(ORBITA_PERIODO / PI)


def acoplamiento_estacional(inercia):
    """k (W/m2/K) -- ver la nota de cabecera sobre el criterio de
    'media potencia' usado para elegir este valor, y sobre por que
    SOLO se aplica a tierra (para agua, ver AGUA_SIN_CAPA_PROFUNDA
    mas abajo)."""
    return capacidad_profunda(inercia) * FRECUENCIA_ANGULAR_ANUAL


# ================================================================
# INVESTIGACION ADICIONAL (tras revisar una conversacion previa sobre
# los valores de INERCIA_POR_TIPO): esta arquitectura de dos capas NO
# se debe aplicar de la misma forma a tierra y a agua. Resultado de la
# investigacion:
#
# AGUA: la conversacion previa establecio que el mecanismo fisico del
# agua no es conduccion (por eso NO tiene sentido escalar su inercia
# por sqrt(periodo) como se hace con la tierra -- la profundidad de
# mezcla turbulenta del oceano no depende del periodo de forzamiento
# de la misma forma que la profundidad de penetracion por conduccion).
# Comprobacion numerica directa: la capacidad calorifica DIURNA que ya
# tiene el agua ahora mismo, calculada con la formula existente,
#     C_diurno_agua = INERCIA_POR_TIPO[AGUA] * sqrt(ROTACION_PERIODO/PI)
#                   = 1.200.000 * sqrt(86400/PI) ~= 1,99e8 J/m2/K
# cae casi exactamente dentro del rango que da la literatura para la
# capacidad calorifica real del oceano (Williams & Kasting 1997 y
# fisica de la capa de mezcla oceanica: 1,2e8 - 2,1e8 J/m2/K). Es
# decir: la capacidad diurna del agua YA representa razonablemente
# bien su capacidad fisica real -- no hace falta (y seria un error,
# el mismo tipo de error que ya se detecto una vez con este mismo
# parametro) darle una capa profunda adicional escalada por
# sqrt(año/dia), porque eso la infla ~16 veces por encima de su valor
# fisico real. Por eso, para AGUA, k = 0 (sin acoplamiento -- el agua
# se queda en el modelo de una sola capa, ya validado).
#
# TIERRA: la busqueda encontro el precedente real y bien establecido
# para esto -- el "metodo force-restore" (Bhumralkar 1975; Blackadar
# 1976), usado en esquemas de superficie terrestre reales (p.ej. ISBA),
# que acopla una capa de superficie rapida con una capa mas profunda
# mediante un termino de relajacion de la forma (tasa) * (T_superficie
# - T_profunda), con tasa del orden de 1/periodo -- estructuralmente
# el mismo enfoque que uso aqui (mi eleccion de "media potencia",
# tasa = 2*PI/periodo, es del mismo orden de magnitud, con un factor
# ~pi de diferencia frente al 2/periodo de la formulacion clasica).
#
# OJO, limite real de esta aproximacion para tierra (encontrado en la
# propia literatura, no una duda mia): los estudios que aplican el
# force-restore clasico señalan que un esquema de SOLO DOS capas
# acumula un sesgo sistematico (~5.5 K en un estudio real) cuando se
# le pide representar a la vez el ciclo diurno Y el estacional -- por
# eso existen extensiones de TRES capas en la literatura ("Inclusion
# of a Third Soil Layer in a Land Surface Scheme Using the
# Force-Restore Method", 1999). Es decir: el enfoque de dos capas que
# he construido es una aproximacion razonable y con respaldo en la
# literatura, pero no es el tratamiento mas riguroso posible -- ese
# seria de tres capas. Lo dejo como esta (dos capas) por ahora, y lo
# señalo para que Carlos decida si en algun momento merece la pena
# ese paso adicional.
# ================================================================

AGUA_SIN_CAPA_PROFUNDA = True  # ver nota de arriba: para agua, k = 0


def _paso_acoplamiento(T_sup_estrella, T_prof, C_sup, C_prof, k, paso_tiempo):
    """Resuelve el sistema lineal 2x2 (implicito) del intercambio de
    calor entre las dos capas, para un paso de tiempo. Vectorizable:
    T_sup_estrella, T_prof, C_sup, C_prof, k pueden ser escalares (un
    punto) o arrays de la misma forma (una rejilla completa)."""
    a = C_sup / paso_tiempo + k
    b = -k
    c = -k
    d = C_prof / paso_tiempo + k
    rhs1 = (C_sup / paso_tiempo) * T_sup_estrella
    rhs2 = (C_prof / paso_tiempo) * T_prof
    det = a * d - b * c
    T_sup_nuevo = (rhs1 * d - b * rhs2) / det
    T_prof_nuevo = (a * rhs2 - c * rhs1) / det
    return T_sup_nuevo, T_prof_nuevo


def simular_dos_capas(
    latitud_rad, datos_orbita, emisividad, inercia, albedo, profundidad_optica,
    paso_tiempo=PASO_TIEMPO, max_anos=50, tolerancia_convergencia=0.01, k_forzado=None,
):
    """
    Version de un unico punto (como temperatura.py:simular(), pero con
    dos capas). k_forzado permite sobreescribir el acoplamiento
    calculado -- se usa en el autotest para comprobar que k=0 (capas
    totalmente desacopladas) reproduce EXACTO el modelo de una capa.

    Devuelve: T_superficie final, T_profunda final, registro_superficie
    (lista de temperaturas, un valor por paso de tiempo, del ultimo
    año simulado), registro_profunda (igual, de la capa profunda).
    """
    C_sup = inercia * math.sqrt(ROTACION_PERIODO / PI)
    C_prof = capacidad_profunda(inercia)
    k = acoplamiento_estacional(inercia) if k_forzado is None else k_forzado

    T_sup = estimar_T_inicial_equilibrio_punto(latitud_rad, datos_orbita, emisividad, albedo, profundidad_optica)
    T_prof = T_sup

    registro_sup = []
    registro_prof = []

    for ano in range(max_anos):
        T_sup_inicio, T_prof_inicio = T_sup, T_prof
        registro_sup = []
        registro_prof = []
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

            emitido = (1 - emisividad / 2) * CONSTANTE_SB * T_sup**4
            T_sup_estrella = T_sup + (paso_tiempo / C_sup) * (abs_local - emitido)

            T_sup, T_prof = _paso_acoplamiento(T_sup_estrella, T_prof, C_sup, C_prof, k, paso_tiempo)

            registro_sup.append(T_sup - 273.15)
            registro_prof.append(T_prof - 273.15)

        if (abs(T_sup - T_sup_inicio) < tolerancia_convergencia
                and abs(T_prof - T_prof_inicio) < tolerancia_convergencia):
            break

    return T_sup - 273.15, T_prof - 273.15, registro_sup, registro_prof


def simular_rejilla_dos_capas(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica,
    paso_tiempo=PASO_TIEMPO, max_anos=50, tolerancia_convergencia=0.01,
):
    """
    Version de rejilla completa de simular_dos_capas(): aplica la capa
    profunda SOLO a las celdas de tierra (AGUA_SIN_CAPA_PROFUNDA -- ver
    la nota de investigacion mas arriba). Las celdas de agua se
    comportan exactamente igual que en el modelo de una capa (Fase 1),
    porque su k_grid es 0 ahi.

    Se apoya en las mismas funciones de radiacion que fase1_geografia.py
    (import local para no crear una dependencia circular al nivel del
    modulo).
    """
    from fase1_geografia import (
        TIERRA, AGUA, correccion_altitud, irradiancia_absorbida_desde_toa_rejilla_con_albedo,
    )

    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])

    C_sup_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)
    C_prof_grid = capacidad_profunda(inercia_grid)

    # k_grid: solo tierra tiene capa profunda activa (ver nota de investigacion).
    k_grid = np.where(
        tipo_superficie == TIERRA,
        acoplamiento_estacional(inercia_por_tipo[TIERRA]),
        0.0,
    )

    from fase1_geografia import estimar_T_inicial_equilibrio
    T_sup = estimar_T_inicial_equilibrio(datos_orbita, albedo_grid, emisividad, profundidad_optica)
    T_prof = T_sup.copy()

    for ano in range(max_anos):
        T_sup_inicio, T_prof_inicio = T_sup.copy(), T_prof.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            emitido = (1 - emisividad / 2) * CONSTANTE_SB * T_sup**4
            T_sup_estrella = T_sup + (paso_tiempo / C_sup_grid) * (abs_local - emitido)

            T_sup, T_prof = _paso_acoplamiento(T_sup_estrella, T_prof, C_sup_grid, C_prof_grid, k_grid, paso_tiempo)

        diferencia_maxima = max(
            np.max(np.abs(T_sup - T_sup_inicio)),
            np.max(np.abs(T_prof - T_prof_inicio)),
        )
        if diferencia_maxima < tolerancia_convergencia:
            return correccion_altitud(T_sup - 273.15, altitud_metros), correccion_altitud(T_prof - 273.15, altitud_metros), ano + 1

    return correccion_altitud(T_sup - 273.15, altitud_metros), correccion_altitud(T_prof - 273.15, altitud_metros), max_anos


if __name__ == "__main__":
    import numpy as np
    from parametros import EMISIVIDAD, PROFUNDIDAD_OPTICA, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR
    from temperatura import simular as simular_una_capa

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    latitud_prueba_rad = math.radians(45)  # latitud media, donde el ciclo estacional es claro
    albedo_prueba = 0.30
    inercia_prueba = 2500  # valor de tierra, INERCIA_POR_TIPO[TIERRA]

    # ---- Prueba 0: k=0 (capas desacopladas) debe reproducir EXACTO el modelo de una capa ----
    print("\n--- Prueba 0: acoplamiento k=0 debe coincidir EXACTO con el modelo de una capa ---")
    T_final_1capa, registro_1capa = simular_una_capa(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA
    )
    T_sup_final, T_prof_final, registro_sup, registro_prof = simular_dos_capas(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA,
        k_forzado=0.0,
    )
    diferencia = max(abs(a - b) for a, b in zip(registro_1capa, registro_sup))
    print(f"Diferencia maxima en el registro de un año: {diferencia:.2e} C")
    print("OK: coincide EXACTA." if diferencia < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    # ---- Prueba 1: conservacion de energia del acoplamiento en aislamiento (sin radiacion) ----
    print("\n--- Prueba 1: el acoplamiento por si solo conserva energia y relaja hacia el equilibrio ---")
    C_sup = inercia_prueba * math.sqrt(ROTACION_PERIODO / PI)
    C_prof = capacidad_profunda(inercia_prueba)
    k = acoplamiento_estacional(inercia_prueba)
    T_sup_prueba, T_prof_prueba = 320.0, 260.0  # arranque muy desequilibrado, a proposito
    energia_inicial = C_sup * T_sup_prueba + C_prof * T_prof_prueba
    for _ in range(20000):  # muchos pasos, sin ninguna radiacion, para ver si relaja
        T_sup_prueba, T_prof_prueba = _paso_acoplamiento(T_sup_prueba, T_prof_prueba, C_sup, C_prof, k, PASO_TIEMPO)
    energia_final = C_sup * T_sup_prueba + C_prof * T_prof_prueba
    diferencia_energia_relativa = abs(energia_final - energia_inicial) / abs(energia_inicial)
    print(f"T_superficie: 320.0 -> {T_sup_prueba:.2f} K | T_profunda: 260.0 -> {T_prof_prueba:.2f} K")
    print(f"Diferencia de temperaturas al final: {abs(T_sup_prueba - T_prof_prueba):.4f} K (deberia ser ~0, ya equilibradas)")
    print(f"Energia (C_sup*T_sup + C_prof*T_prof): diferencia relativa {diferencia_energia_relativa:.2e}")
    print("OK: conserva energia y relaja al equilibrio." if diferencia_energia_relativa < 1e-9 and abs(T_sup_prueba - T_prof_prueba) < 0.01
          else "AVISO: revisar.")

    # ---- Prueba 2: la capa profunda SI amortigua el ciclo estacional en la superficie ----
    print("\n--- Prueba 2: comparacion de la amplitud estacional, una capa vs dos capas ---")
    pasos_por_dia = round(ROTACION_PERIODO / PASO_TIEMPO)
    num_dias = len(registro_1capa) // pasos_por_dia

    def medias_diarias(registro):
        return np.array([
            np.mean(registro[d * pasos_por_dia:(d + 1) * pasos_por_dia])
            for d in range(num_dias)
        ])

    medias_1capa = medias_diarias(registro_1capa)
    T_sup_final2, T_prof_final2, registro_sup2, registro_prof2 = simular_dos_capas(
        latitud_prueba_rad, datos_orbita, EMISIVIDAD, inercia_prueba, albedo_prueba, PROFUNDIDAD_OPTICA,
    )
    medias_2capas = medias_diarias(registro_sup2)

    amplitud_1capa = medias_1capa.max() - medias_1capa.min()
    amplitud_2capas = medias_2capas.max() - medias_2capas.min()
    print(f"Amplitud estacional (max dia medio - min dia medio) -- una capa: {amplitud_1capa:.2f} C")
    print(f"Amplitud estacional -- dos capas (con inercia estacional): {amplitud_2capas:.2f} C")
    print(f"Reduccion: {100 * (1 - amplitud_2capas / amplitud_1capa):.1f}%")
    print(f"(Para referencia, la capa profunda alcanzo un rango de {min(registro_prof2):.1f} C a {max(registro_prof2):.1f} C en el año)")
    print("OK: la capa profunda amortigua el ciclo estacional (amplitud menor)." if amplitud_2capas < amplitud_1capa
          else "AVISO: no se aprecia amortiguacion -- revisar.")

    # ---- Prueba 3: rejilla completa, mapa mixto -- el agua NO debe cambiar
    #                respecto a Fase 1 (k=0 ahi), la tierra SI debe amortiguarse ----
    print("\n--- Prueba 3: rejilla completa (mapa mixto) -- agua sin cambios, tierra amortiguada ---")
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO, TIERRA, AGUA, simular_rejilla_geografia
    from fase1_mapa_mixto import mapa_mitad_agua_mitad_tierra

    tipo_mixto, altitud_mixta = mapa_mitad_agua_mitad_tierra()

    T_f1, anos_f1 = simular_rejilla_geografia(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )
    T_sup_final3, T_prof_final3, anos_f2c = simular_rejilla_dos_capas(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA
    )

    mascara_agua = tipo_mixto == AGUA
    mascara_tierra = tipo_mixto == TIERRA

    # En agua, k=0, asi que esta version deberia coincidir EXACTA con Fase 1
    # (mismo instante congelado del ciclo, mismo criterio de convergencia):
    diferencia_agua = np.max(np.abs(T_sup_final3[mascara_agua] - T_f1[mascara_agua]))
    print(f"Agua -- diferencia maxima vs Fase 1 (deberia ser exacta, k=0 ahi): {diferencia_agua:.2e} C")
    print("OK: coincide EXACTA en agua." if diferencia_agua < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    print(f"Tierra -- capa superficial en este instante: "
          f"{T_sup_final3[mascara_tierra].min():.1f} C a {T_sup_final3[mascara_tierra].max():.1f} C")
    print(f"Tierra -- capa profunda en este instante: "
          f"{T_prof_final3[mascara_tierra].min():.1f} C a {T_prof_final3[mascara_tierra].max():.1f} C")
    variacion_sup = T_sup_final3[mascara_tierra].max() - T_sup_final3[mascara_tierra].min()
    variacion_prof = T_prof_final3[mascara_tierra].max() - T_prof_final3[mascara_tierra].min()
    print(f"(Variacion espacial en tierra -- superficie: {variacion_sup:.1f} C | profunda: {variacion_prof:.1f} C)")

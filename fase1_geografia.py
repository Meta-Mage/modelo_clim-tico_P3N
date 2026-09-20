# fase1_geografia.py -- Fase 1: geografia real (agua/tierra + altitud por celda).
#
# Primer paso: un mapa DE PRUEBA (todo tierra, a nivel del mar) para
# demostrar que el mecanismo nuevo (albedo/inercia por celda + correccion
# por altitud) se reduce exactamente a la Fase 0 cuando no hay ninguna
# diferencia real que aplicar. Los valores FISICOS reales de tierra/agua
# (investigados aparte) se adoptan en un paso posterior, una vez validada
# la mecanica -- mezclarlos aqui impediria saber si un fallo viene de la
# mecanica o de haber cambiado los numeros.

import math
import numpy as np
from parametros import ROTACION_PERIODO, PI, CONSTANTE_SB
from temperatura import precalcular_orbita, t_eq, PASO_TIEMPO as PASO_TIEMPO_POR_DEFECTO
from fase0_radiacion import irradiancia_absorbida_desde_toa_rejilla
from fase0_geometria import angulo_cenital_rejilla
from fase0_atmosfera import masa_aire_rejilla
from fase0_temperatura import simular_rejilla
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS, FILAS, COLUMNAS

TIERRA = 1
AGUA = 0

GRADIENTE_TERMICO = 6.5  # C por km de altitud (correccion estandar simplificada)


def mapa_falso_todo_tierra():
    """Mapa de prueba: todas las celdas son tierra, todas a nivel del mar."""
    tipo_superficie = np.full((FILAS, COLUMNAS), TIERRA)
    altitud_metros = np.zeros((FILAS, COLUMNAS))
    return tipo_superficie, altitud_metros


def correccion_altitud(T_grados_C, altitud_metros):
    return T_grados_C - GRADIENTE_TERMICO * (altitud_metros / 1000)


def simular_rejilla_geografia(datos_orbita, tipo_superficie, altitud_metros, emisividad,
                               albedo_por_tipo, inercia_por_tipo, profundidad_optica,
                               paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.01):
    """
    Version de simular_rejilla() (fase0_temperatura.py) con albedo e
    inercia termica propios de cada celda segun su tipo (agua/tierra), y
    correccion de temperatura por altitud al final.
    """
    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])
    C_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)

    toa_ini, decl_ini, ang_h_ini_lon0 = datos_orbita[len(datos_orbita) // 2]
    ang_h_ini_grid = ang_h_ini_lon0 + np.radians(LONGITUDES_GRADOS)
    cenital_ini = angulo_cenital_rejilla(LATITUDES_GRADOS, decl_ini, ang_h_ini_grid)
    masa_ini = masa_aire_rejilla(cenital_ini)
    es_de_noche_ini = np.isnan(masa_ini)

    # abs_ini necesita el albedo de cada celda -- i_abs() es aritmetica
    # pura, asi que ya acepta un array de albedo sin ningun cambio.
    abs_ini = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
        toa_ini, decl_ini, ang_h_ini_lon0, albedo_grid, profundidad_optica
    )
    T = np.where(es_de_noche_ini, 273.15, t_eq(abs_ini))

    for ano in range(max_anos):
        T_inicio_ano = T.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            dT = (paso_tiempo / C_grid) * (abs_local - (1 - emisividad / 2) * CONSTANTE_SB * T ** 4)
            T = T + dT

        diferencia_maxima = np.max(np.abs(T - T_inicio_ano))
        if diferencia_maxima < tolerancia_convergencia:
            return correccion_altitud(T - 273.15, altitud_metros), ano + 1

    return correccion_altitud(T - 273.15, altitud_metros), max_anos


def irradiancia_absorbida_desde_toa_rejilla_con_albedo(toa, declinacion, ang_h_lon0, albedo_grid, profundidad_optica):
    """
    Igual que irradiancia_absorbida_desde_toa_rejilla() (fase0_radiacion.py)
    pero con albedo por celda en vez de un unico valor -- se reimplementa
    aqui en vez de tocar la version de la Fase 0 (que ya quedo validada tal
    cual, con albedo unico, y no hay que arriesgarla).
    """
    from fase0_geometria import angulo_cenital_rejilla as _acr
    from fase0_atmosfera import trans_rejilla as _tr
    from radiacion import i_atm, i_abs
    from fase0_radiacion import i_inst_rejilla as _iir

    ang_h_grid = ang_h_lon0 + np.radians(LONGITUDES_GRADOS)
    cenital = _acr(LATITUDES_GRADOS, declinacion, ang_h_grid)
    inst = _iir(toa, cenital)
    masa = masa_aire_rejilla(cenital)
    tra = _tr(masa, profundidad_optica)
    atm = i_atm(inst, tra)
    absorbida = i_abs(atm, albedo_grid)  # i_abs es aritmetica pura -- acepta array de albedo sin cambios
    return np.nan_to_num(absorbida, nan=0.0)


if __name__ == "__main__":
    from parametros import EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    print("Simulando con Fase 0 (sin geografia, referencia)...")
    T_fase0, anos_fase0 = simular_rejilla(datos_orbita, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)

    print("Simulando con Fase 1, mapa falso todo tierra a nivel del mar (mismos valores que Fase 0)...")
    tipo_superficie, altitud_metros = mapa_falso_todo_tierra()
    albedo_por_tipo = {TIERRA: ALBEDO, AGUA: ALBEDO}      # AGUA no se usa (mapa sin agua); mismo valor por seguridad
    inercia_por_tipo = {TIERRA: INERCIA_TERMICA, AGUA: INERCIA_TERMICA}
    T_fase1, anos_fase1 = simular_rejilla_geografia(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, albedo_por_tipo, inercia_por_tipo, PROFUNDIDAD_OPTICA
    )

    diferencia = np.abs(T_fase0 - T_fase1)
    print(f"Convergencia: {anos_fase0} año(s) (Fase 0) vs {anos_fase1} año(s) (Fase 1, mapa falso)")
    print(f"Diferencia maxima entre Fase 0 y Fase 1 con mapa todo tierra: {diferencia.max():.2e} C")
    if diferencia.max() < 1e-6:
        print("OK: con mapa falso todo tierra a nivel del mar, la Fase 1 coincide EXACTA con la Fase 0.")
    else:
        print("AVISO: deberian coincidir exactas -- revisar antes de seguir.")


# ============================================================================
# VALORES DE REFERENCIA FIJADOS (Fase 1, 20/09/2026)
# ============================================================================
# Albedo: valores de investigacion academica (NASA MyNASAData, literatura
# MODIS de albedo de vegetacion/suelo). Tierra: valor medio comunmente
# citado. Oceano: punto medio del rango tipico de albedo oceanico real
# (0.06-0.10).
#
# Inercia termica:
#   TIERRA (2500): dentro del rango real de inercia termica de suelo/roca
#     (400-2500 tiu, J m-2 K-1 s-0.5), elegido de forma prudente en el
#     extremo alto de ese rango.
#   AGUA (1200000): NO es una inercia termica de material en sentido
#     estricto -- el oceano no se calienta por conduccion como la tierra,
#     sino por mezcla turbulenta de viento/oleaje de una capa de decenas
#     de metros, cuyo espesor no escala con el periodo de forzamiento
#     igual que la conduccion en tierra. Es un valor EFECTIVO, calculado
#     para que la formula C = inercia * sqrt(periodo_rotacion / pi)
#     reproduzca la capacidad calorifica real de una capa oceanica
#     mezclada (~10^7-10^8 J/m2/K, contrastado con Williams & Kasting 1997
#     y Bhattacharya & Bordoni 2020).
#     PROVISIONAL: pendiente de revision cuando se modele la fisica
#     oceanica (mezcla, corrientes, profundidad variable por latitud y
#     estacion) con mas detalle.
# ============================================================================

ALBEDO_POR_TIPO = {
    TIERRA: 0.20,
    AGUA: 0.08,
}

INERCIA_POR_TIPO = {
    TIERRA: 2500,
    AGUA: 1200000,
}


# ============================================================================
# MOTOR DE REGISTRO ANUAL (minima, media y maxima diaria de la rejilla)
# ============================================================================
# Igual que simular_rejilla_geografia(), pero una vez alcanzada la
# convergencia, simula UN año adicional y para cada dia de ese año
# acumula todos los pasos de tiempo dentro de ese dia para calcular la
# temperatura minima, media y maxima de cada celda -- no una sola foto
# instantanea, sino el ciclo dia/noche completo de cada dia. Estas tres
# "peliculas" (minima, media, maxima) son la base comun para: consultas
# de un punto concreto, tablas por latitud, graficos globales y
# animaciones -- todo se deriva de este mismo registro, sin recalcular
# la simulacion cada vez.
# ============================================================================

def simular_rejilla_geografia_con_registro(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.01,
):
    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])
    inercia_grid = np.where(tipo_superficie == TIERRA, inercia_por_tipo[TIERRA], inercia_por_tipo[AGUA])
    C_grid = inercia_grid * math.sqrt(ROTACION_PERIODO / PI)

    toa_ini, decl_ini, ang_h_ini_lon0 = datos_orbita[len(datos_orbita) // 2]
    ang_h_ini_grid = ang_h_ini_lon0 + np.radians(LONGITUDES_GRADOS)
    cenital_ini = angulo_cenital_rejilla(LATITUDES_GRADOS, decl_ini, ang_h_ini_grid)
    masa_ini = masa_aire_rejilla(cenital_ini)
    es_de_noche_ini = np.isnan(masa_ini)
    abs_ini = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
        toa_ini, decl_ini, ang_h_ini_lon0, albedo_grid, profundidad_optica
    )
    T = np.where(es_de_noche_ini, 273.15, t_eq(abs_ini))

    for ano in range(max_anos):
        T_inicio_ano = T.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            dT = (paso_tiempo / C_grid) * (abs_local - (1 - emisividad / 2) * CONSTANTE_SB * T ** 4)
            T = T + dT
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
        dT = (paso_tiempo / C_grid) * (abs_local - (1 - emisividad / 2) * CONSTANTE_SB * T ** 4)
        T = T + dT

        T_min_dia = np.minimum(T_min_dia, T)
        T_max_dia = np.maximum(T_max_dia, T)
        suma_dia = suma_dia + T
        contador_dia += 1

    if contador_dia > 0:
        guardar_dia()

    registro_minima = np.array(registro_minima)
    registro_media = np.array(registro_media)
    registro_maxima = np.array(registro_maxima)

    T_final_corregida = correccion_altitud(T - 273.15, altitud_metros)

    return T_final_corregida, anos_convergencia, registro_minima, registro_media, registro_maxima


if __name__ == "__main__":
    from temperatura import precalcular_orbita
    from parametros import ORBITA_PERIODO

    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    tipo_tierra, altitud_cero = mapa_falso_todo_tierra()

    T_normal, anos_normal = simular_rejilla_geografia(
        datos_orbita, tipo_tierra, altitud_cero, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
    )

    T_registro, anos_registro, registro_minima, registro_media, registro_maxima = simular_rejilla_geografia_con_registro(
        datos_orbita, tipo_tierra, altitud_cero, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
    )

    print(f"Convergencia (sin registro): {anos_normal} año(s)")
    print(f"Convergencia (con registro): {anos_registro} año(s)")
    print(f"Forma de cada registro (minima/media/maxima): {registro_media.shape} (dias, filas, columnas)")

    dias_esperados_aprox = ORBITA_PERIODO / ROTACION_PERIODO
    print(f"Dias esperados aproximadamente: {dias_esperados_aprox:.1f}")

    orden_correcto = np.all(registro_minima <= registro_media + 1e-9) and np.all(registro_media <= registro_maxima + 1e-9)
    print(f"Orden minima <= media <= maxima en todas las celdas y dias: {orden_correcto}")

    dentro_del_rango = np.all(T_registro >= registro_minima[-1] - 1e-6) and np.all(T_registro <= registro_maxima[-1] + 1e-6)
    print(f"T final dentro del rango [minima, maxima] del ultimo dia registrado: {dentro_del_rango}")

    forma_correcta = abs(registro_media.shape[0] - dias_esperados_aprox) < 2
    if orden_correcto and forma_correcta and dentro_del_rango:
        print("OK: el registro diario (minima/media/maxima) es coherente.")
    else:
        print("AVISO: revisar antes de seguir.")

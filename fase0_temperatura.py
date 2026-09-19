# fase0_temperatura.py -- Fase 0: simular() vectorizado sobre toda la rejilla.
#
# La pieza mas delicada de la Fase 0 (convergencia conjunta, sin registro
# por celda, instante inicial de noche -> 273.15 K -- ver el mensaje que
# acompaño la primera version de este archivo para la explicacion completa).
#
# NUEVO: paso_tiempo ahora es un parametro (por defecto, el PASO_TIEMPO de
# temperatura.py) en vez de estar fijado por dentro -- necesario para poder
# repetir la simulacion con otro paso de tiempo y comprobar sensibilidad,
# sin duplicar la funcion entera para eso.

import math
import numpy as np
from parametros import ROTACION_PERIODO, PI, CONSTANTE_SB
from temperatura import precalcular_orbita, t_eq, PASO_TIEMPO as PASO_TIEMPO_POR_DEFECTO
from fase0_radiacion import irradiancia_absorbida_desde_toa_rejilla
from fase0_geometria import angulo_cenital_rejilla
from fase0_atmosfera import masa_aire_rejilla
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS, FILAS, COLUMNAS


def simular_rejilla(datos_orbita, emisividad, inercia, albedo, profundidad_optica,
                     paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.01):
    """
    Version de simular() (temperatura.py) sobre toda la rejilla lat/lon a
    la vez. Devuelve (T_final_grados_C, anos_hasta_converger).
    """
    C = inercia * math.sqrt(ROTACION_PERIODO / PI)

    toa_ini, decl_ini, ang_h_ini_lon0 = datos_orbita[len(datos_orbita) // 2]
    ang_h_ini_grid = ang_h_ini_lon0 + np.radians(LONGITUDES_GRADOS)
    cenital_ini = angulo_cenital_rejilla(LATITUDES_GRADOS, decl_ini, ang_h_ini_grid)
    masa_ini = masa_aire_rejilla(cenital_ini)
    es_de_noche_ini = np.isnan(masa_ini)

    abs_ini = irradiancia_absorbida_desde_toa_rejilla(
        toa_ini, decl_ini, ang_h_ini_lon0, albedo, profundidad_optica
    )
    T = np.where(es_de_noche_ini, 273.15, t_eq(abs_ini))

    for ano in range(max_anos):
        T_inicio_ano = T.copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla(
                toa, decl, ang_h_lon0, albedo, profundidad_optica
            )
            dT = (paso_tiempo / C) * (abs_local - (1 - emisividad / 2) * CONSTANTE_SB * T ** 4)
            T = T + dT

        diferencia_maxima = np.max(np.abs(T - T_inicio_ano))
        if diferencia_maxima < tolerancia_convergencia:
            return T - 273.15, ano + 1

    return T - 273.15, max_anos


if __name__ == "__main__":
    from temperatura import simular
    from parametros import (EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA,
                             INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR)

    print("Precalculando orbita (un año, decenas de miles de pasos)...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    print(f"Simulando rejilla completa ({FILAS}x{COLUMNAS} celdas) hasta que converja...")
    T_grid, anos = simular_rejilla(datos_orbita, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
    print(f"Convergencia conjunta en {anos} año(s).")

    celdas_de_prueba = [(0, 0), (9, 18), (17, 36), (18, 54), (35, 71)]
    max_diferencia = 0.0
    for i, j in celdas_de_prueba:
        lat = LATITUDES_GRADOS[i]
        lon = LONGITUDES_GRADOS[j]
        lat_rad = math.radians(lat)
        lon_rad = math.radians(lon)
        datos_orbita_desplazados = [(toa, decl, ang_h + lon_rad) for (toa, decl, ang_h) in datos_orbita]
        T_punto, _ = simular(lat_rad, datos_orbita_desplazados, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        diferencia = abs(T_punto - T_grid[i, j])
        max_diferencia = max(max_diferencia, diferencia)
        print(f"  Celda ({i},{j}) lat={lat:.1f} lon={lon:.1f}: rejilla={T_grid[i,j]:.4f} C | punto={T_punto:.4f} C | diff={diferencia:.2e}")

    print(f"Diferencia maxima entre rejilla y punto (temperatura final, C): {max_diferencia:.2e}")
    if max_diferencia < 0.02:
        print("OK: la simulacion vectorizada coincide con el modelo actual en las celdas probadas.")
    else:
        print("AVISO: revisar antes de seguir.")

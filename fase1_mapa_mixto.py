# fase1_mapa_mixto.py -- Fase 1: prueba con un mapa que SI tiene agua.
#
# Valores de este paso: PLACEHOLDER, deliberadamente distintos entre si
# para que la prueba sea inequivoca -- NO son todavia los valores reales
# investigados (eso es un paso aparte, pendiente de resolver antes la
# inercia termica del oceano).
#
# Mapa: mitad occidental (columnas 0-35) agua, mitad oriental (columnas
# 36-71) tierra, todo a nivel del mar. Como las celdas todavia no se
# influyen entre si (sin difusion de calor -- eso es la Fase 2), la forma
# mas rigurosa de validar es comparar cada mitad por separado contra una
# simulacion de un unico valor global (como la Fase 0) con ese mismo valor.

import numpy as np
from temperatura import precalcular_orbita
from fase0_temperatura import simular_rejilla
from fase1_geografia import simular_rejilla_geografia, TIERRA, AGUA
from rejilla import FILAS, COLUMNAS

ALBEDO_TIERRA_PRUEBA = 0.3      # igual que Fase 0, para reusar esa simulacion como referencia
INERCIA_TIERRA_PRUEBA = 3500
ALBEDO_AGUA_PRUEBA = 0.1        # deliberadamente distinto -- NO es el valor real investigado
INERCIA_AGUA_PRUEBA = 1000      # deliberadamente distinto -- NO es el valor real investigado


def mapa_mitad_agua_mitad_tierra():
    tipo_superficie = np.full((FILAS, COLUMNAS), TIERRA)
    tipo_superficie[:, :COLUMNAS // 2] = AGUA
    altitud_metros = np.zeros((FILAS, COLUMNAS))
    return tipo_superficie, altitud_metros


if __name__ == "__main__":
    from parametros import EMISIVIDAD, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR, PROFUNDIDAD_OPTICA

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    print("Simulando mapa mixto (mitad agua, mitad tierra)...")
    tipo_superficie, altitud_metros = mapa_mitad_agua_mitad_tierra()
    albedo_por_tipo = {TIERRA: ALBEDO_TIERRA_PRUEBA, AGUA: ALBEDO_AGUA_PRUEBA}
    inercia_por_tipo = {TIERRA: INERCIA_TIERRA_PRUEBA, AGUA: INERCIA_AGUA_PRUEBA}
    T_mixto, anos_mixto = simular_rejilla_geografia(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD, albedo_por_tipo, inercia_por_tipo, PROFUNDIDAD_OPTICA
    )

    print("Simulando referencia 'todo tierra' con el valor de tierra...")
    T_solo_tierra, _ = simular_rejilla(datos_orbita, EMISIVIDAD, INERCIA_TIERRA_PRUEBA, ALBEDO_TIERRA_PRUEBA, PROFUNDIDAD_OPTICA)

    print("Simulando referencia 'todo agua' con el valor de agua...")
    T_solo_agua, _ = simular_rejilla(datos_orbita, EMISIVIDAD, INERCIA_AGUA_PRUEBA, ALBEDO_AGUA_PRUEBA, PROFUNDIDAD_OPTICA)

    mitad_agua = tipo_superficie == AGUA
    mitad_tierra = tipo_superficie == TIERRA

    diferencia_tierra = np.max(np.abs(T_mixto[mitad_tierra] - T_solo_tierra[mitad_tierra]))
    diferencia_agua = np.max(np.abs(T_mixto[mitad_agua] - T_solo_agua[mitad_agua]))

    print(f"Convergencia del mapa mixto: {anos_mixto} año(s)")
    print(f"Diferencia maxima, celdas de tierra del mapa mixto vs 'todo tierra': {diferencia_tierra:.2e} C")
    print(f"Diferencia maxima, celdas de agua del mapa mixto vs 'todo agua': {diferencia_agua:.2e} C")

    ejemplo_frontera_agua = T_mixto[18, COLUMNAS // 2 - 1]
    ejemplo_frontera_tierra = T_mixto[18, COLUMNAS // 2]
    print(f"Ejemplo en la frontera: celda de agua = {ejemplo_frontera_agua:.2f} C | celda de tierra vecina = {ejemplo_frontera_tierra:.2f} C")

    if diferencia_tierra < 0.01 and diferencia_agua < 0.01:
        print("OK: cada celda del mapa mixto se comporta exactamente segun su propio tipo.")
    else:
        print("AVISO: revisar antes de seguir.")

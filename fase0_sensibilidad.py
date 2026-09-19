# fase0_sensibilidad.py -- Fase 0: sensibilidad al paso de tiempo.
#
# Repite simular_rejilla() con la mitad de PASO_TIEMPO (450 s en vez de
# 900 s) y compara el resultado final con el paso normal. Si la fisica
# esta bien planteada, deberian coincidir dentro de un margen razonable;
# una diferencia grande indicaria que el paso actual es demasiado grueso
# para la inercia termica que se esta usando.
#
# precalcular_orbita() (temperatura.py) usa el PASO_TIEMPO global de ese
# modulo para precalcular los pasos de la orbita -- no esta pensado para
# recibir el paso como argumento. Aqui se cambia ese valor temporalmente,
# se llama, y se restaura -- no se toca el archivo original.

import numpy as np
import temperatura
from temperatura import precalcular_orbita
from fase0_temperatura import simular_rejilla
from parametros import (EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA,
                         INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR)


def precalcular_orbita_con_paso(paso_tiempo):
    paso_original = temperatura.PASO_TIEMPO
    temperatura.PASO_TIEMPO = paso_tiempo
    try:
        return precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    finally:
        temperatura.PASO_TIEMPO = paso_original


if __name__ == "__main__":
    paso_normal = temperatura.PASO_TIEMPO
    paso_mitad = paso_normal / 2

    print(f"Simulando con paso de tiempo normal ({paso_normal:.0f} s)...")
    datos_normal = precalcular_orbita_con_paso(paso_normal)
    T_normal, anos_normal = simular_rejilla(
        datos_normal, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA, paso_tiempo=paso_normal
    )

    print(f"Simulando con la mitad del paso de tiempo ({paso_mitad:.0f} s)...")
    datos_mitad = precalcular_orbita_con_paso(paso_mitad)
    T_mitad, anos_mitad = simular_rejilla(
        datos_mitad, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA, paso_tiempo=paso_mitad
    )

    diferencia = np.abs(T_normal - T_mitad)
    print(f"Convergencia: {anos_normal} año(s) (paso normal) vs {anos_mitad} año(s) (paso mitad)")
    print(f"Diferencia maxima entre ambas resoluciones temporales: {diferencia.max():.4f} C")
    print(f"Diferencia media: {diferencia.mean():.4f} C")

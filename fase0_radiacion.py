# fase0_radiacion.py -- Fase 0: i_inst() vectorizado + cadena completa
# (geometria + atmosfera + radiacion) para calcular la irradiancia
# absorbida de toda la rejilla en un instante dado.
#
# i_toa() no depende de la posicion (solo de la distancia al Sol, igual
# para todo el planeta en un instante dado) -- se llama una vez, sin
# cambios. i_atm() e i_abs() son aritmetica pura, ya funcionan sobre
# arrays sin tocarlas. Solo i_inst() usa math.cos y necesita version nueva.

import math
import numpy as np
from geometria import declinacion_solar, angulo_horario, angulo_cenital
from atmosfera import masa_aire, trans
from radiacion import i_toa, i_inst, i_atm, i_abs
from parametros import S3N_LUMINOSIDAD, ALBEDO, PROFUNDIDAD_OPTICA, SEMIEJE_MAYOR
from orbita import distancia_S3N
from fase0_geometria import angulo_cenital_rejilla
from fase0_atmosfera import masa_aire_rejilla, trans_rejilla
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS


def i_inst_rejilla(irradiancia_toa, cenital_rad):
    return irradiancia_toa * np.cos(cenital_rad)


def irradiancia_absorbida_rejilla(distancia, luminosidad, declinacion, hora, albedo, profundidad_optica):
    """
    Encadena, sobre toda la rejilla a la vez, los mismos pasos que hace
    simular() (temperatura.py) punto a punto en cada paso de tiempo:
    angulo cenital -> irradiancia TOA -> irradiancia instantanea ->
    atenuacion atmosferica -> absorcion. De noche, el resultado es 0 --
    igual que hace temperatura.py con `if masa is None: abs_local = 0`.
    """
    angulos_horarios = angulo_horario(hora, LONGITUDES_GRADOS)
    cenital = angulo_cenital_rejilla(LATITUDES_GRADOS, declinacion, angulos_horarios)

    irradiancia_toa = i_toa(distancia, luminosidad)
    irradiancia_instantanea = i_inst_rejilla(irradiancia_toa, cenital)

    masa = masa_aire_rejilla(cenital)
    transmitancia = trans_rejilla(masa, profundidad_optica)

    irradiancia_atenuada = i_atm(irradiancia_instantanea, transmitancia)
    irradiancia_absorbida = i_abs(irradiancia_atenuada, albedo)

    return np.nan_to_num(irradiancia_absorbida, nan=0.0)


if __name__ == "__main__":
    momento_prueba = 0.7
    inclinacion_prueba = math.radians(1.7)
    hora_prueba = 14.5

    declinacion_prueba = declinacion_solar(momento_prueba, inclinacion_prueba)
    distancia_prueba = distancia_S3N(momento_prueba, SEMIEJE_MAYOR)

    resultado_rejilla = irradiancia_absorbida_rejilla(
        distancia_prueba, S3N_LUMINOSIDAD, declinacion_prueba, hora_prueba, ALBEDO, PROFUNDIDAD_OPTICA
    )

    max_diferencia = 0.0
    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h_punto = angulo_horario(hora_prueba, lon)
            cenital_punto = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h_punto)
            toa_punto = i_toa(distancia_prueba, S3N_LUMINOSIDAD)
            inst_punto = i_inst(toa_punto, cenital_punto)
            masa_punto = masa_aire(cenital_punto)

            if masa_punto is None:
                abs_punto = 0.0
            else:
                tra_punto = trans(masa_punto, PROFUNDIDAD_OPTICA)
                atm_punto = i_atm(inst_punto, tra_punto)
                abs_punto = i_abs(atm_punto, ALBEDO)

            max_diferencia = max(max_diferencia, abs(abs_punto - resultado_rejilla[i, j]))

    print(f"Diferencia maxima (irradiancia absorbida, W/m2): {max_diferencia:.2e}")
    if max_diferencia < 1e-6:
        print("OK: la cadena completa vectorizada coincide con el modelo actual, celda a celda.")
    else:
        print("AVISO: revisar antes de seguir.")

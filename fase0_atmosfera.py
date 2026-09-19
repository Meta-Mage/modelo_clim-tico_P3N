# fase0_atmosfera.py -- Fase 0: version vectorizada de masa_aire() y trans().
#
# Novedad real de este paso: con una rejilla de verdad, columnas distintas
# pueden estar de dia y de noche AL MISMO TIEMPO -- en el modelo de un solo
# punto eso nunca pasaba. masa_aire() marca la noche devolviendo None; en
# un array de NumPy no cabe None, asi que aqui se usa np.nan con el mismo
# significado ("de noche, sin luz que atenuar").

import math
import numpy as np
from atmosfera import masa_aire, trans


def masa_aire_rejilla(cenital_rad):
    """
    Version de masa_aire() sobre una rejilla. Donde es de noche
    (altura_solar <= 0), el resultado es np.nan -- mismo papel que el
    None de la version de un punto.
    """
    cenital_grados = np.degrees(cenital_rad)
    altura_solar = 90 - cenital_grados
    with np.errstate(invalid="ignore"):
        valor = 1 / np.sin(np.radians(altura_solar + 244 / (165 + 47 * altura_solar ** 1.1)))
    return np.where(altura_solar <= 0, np.nan, valor)


def trans_rejilla(masa, profundidad_optica):
    """
    Version de trans() sobre una rejilla. np.exp(nan) da nan de forma
    natural, asi que la "noche" (nan en masa) se propaga sola.
    """
    return np.exp(-masa * profundidad_optica)


if __name__ == "__main__":
    from geometria import declinacion_solar, angulo_horario, angulo_cenital
    from fase0_geometria import angulo_cenital_rejilla
    from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS

    momento_prueba = 0.7
    inclinacion_prueba = math.radians(1.7)
    hora_prueba = 14.5
    profundidad_optica_prueba = 0.3

    declinacion_prueba = declinacion_solar(momento_prueba, inclinacion_prueba)
    angulos_horarios_prueba = angulo_horario(hora_prueba, LONGITUDES_GRADOS)
    cenital_rejilla = angulo_cenital_rejilla(LATITUDES_GRADOS, declinacion_prueba, angulos_horarios_prueba)

    masa_rejilla = masa_aire_rejilla(cenital_rejilla)
    trans_calculada = trans_rejilla(masa_rejilla, profundidad_optica_prueba)

    max_diferencia = 0.0
    celdas_de_noche = 0
    celdas_de_dia = 0
    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h_punto = angulo_horario(hora_prueba, lon)
            cenital_punto = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h_punto)
            masa_punto = masa_aire(cenital_punto)

            if masa_punto is None:
                celdas_de_noche += 1
                if not math.isnan(masa_rejilla[i, j]):
                    print(f"AVISO: celda ({i},{j}) deberia ser de noche (nan) y no lo es")
                continue

            celdas_de_dia += 1
            max_diferencia = max(max_diferencia, abs(masa_punto - masa_rejilla[i, j]))
            trans_punto = trans(masa_punto, profundidad_optica_prueba)
            max_diferencia = max(max_diferencia, abs(trans_punto - trans_calculada[i, j]))

    print(f"Celdas de dia: {celdas_de_dia} | Celdas de noche: {celdas_de_noche}")
    print(f"Diferencia maxima (masa de aire y transmitancia): {max_diferencia:.2e}")
    if max_diferencia < 1e-9 and celdas_de_dia > 0 and celdas_de_noche > 0:
        print("OK: coincide con el modelo actual, celda a celda, con dia y noche presentes a la vez.")
    else:
        print("AVISO: revisar antes de seguir.")

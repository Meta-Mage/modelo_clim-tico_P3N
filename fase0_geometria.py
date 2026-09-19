# fase0_geometria.py — Fase 0: version vectorizada del angulo cenital solar.
#
# declinacion_solar() y angulo_horario() (definidas en geometria.py) ya
# funcionan tal cual sobre arrays de NumPy: ninguna de las dos llama a una
# funcion de `math`, solo hacen aritmetica normal. No se tocan ni se
# duplican aqui.
#
# angulo_cenital() SI usa math.sin/cos/acos, que no aceptan arrays (fallan
# con TypeError). Esta es la unica funcion que necesita una version nueva,
# capaz de trabajar sobre la rejilla entera a la vez.

import math
import numpy as np
from geometria import declinacion_solar, angulo_horario, angulo_cenital
from rejilla import LATITUDES_GRADOS, LONGITUDES_GRADOS


def angulo_cenital_rejilla(latitudes_grados, declinacion, angulos_horarios):
    """
    Version de angulo_cenital() que trabaja sobre una rejilla entera a la
    vez, no sobre un unico punto.

    latitudes_grados:  array de forma (filas,)   -- una latitud por fila.
    declinacion:        escalar en radianes (no depende de la posicion).
    angulos_horarios:  array de forma (columnas,) -- un angulo horario por
                        columna (ya con la longitud real de cada una).

    Devuelve un array (filas, columnas) con el angulo cenital de cada
    celda, en radianes. Misma formula exacta que angulo_cenital(), solo
    que evaluada para todas las celdas a la vez mediante broadcasting.
    """
    lat_rad = np.radians(latitudes_grados).reshape(-1, 1)   # (filas, 1)
    ang_h = np.asarray(angulos_horarios).reshape(1, -1)     # (1, columnas)

    cos_cenital = (
        np.sin(lat_rad) * np.sin(declinacion)
        + np.cos(lat_rad) * np.cos(declinacion) * np.cos(ang_h)
    )
    # Recorte defensivo: por redondeo, cos_cenital puede salirse
    # ligerisimamente de [-1, 1] (p. ej. 1.0000000000000002), lo que
    # rompe arccos. Con un unico punto (math.acos) el margen de error no
    # llega a notarse nunca; con miles de celdas a la vez, mejor curarse
    # en salud. Es una diferencia deliberada frente al original, no un
    # cambio de fisica -- se deja anotado a proposito.
    cos_cenital = np.clip(cos_cenital, -1.0, 1.0)
    return np.arccos(cos_cenital)


if __name__ == "__main__":
    # ---- Validacion: ¿coincide con la version de un solo punto? ----
    momento_prueba = 0.7                 # un AV cualquiera, en radianes
    inclinacion_prueba = math.radians(1.7)
    hora_prueba = 14.5

    declinacion_prueba = declinacion_solar(momento_prueba, inclinacion_prueba)

    # angulo_horario ya acepta el array completo de longitudes tal cual,
    # sin ningun cambio -- se demuestra usandolo asi directamente.
    angulos_horarios_prueba = angulo_horario(hora_prueba, LONGITUDES_GRADOS)

    resultado_rejilla = angulo_cenital_rejilla(
        LATITUDES_GRADOS, declinacion_prueba, angulos_horarios_prueba
    )

    max_diferencia = 0.0
    for i, lat in enumerate(LATITUDES_GRADOS):
        for j, lon in enumerate(LONGITUDES_GRADOS):
            ang_h_punto = angulo_horario(hora_prueba, lon)
            valor_punto = angulo_cenital(math.radians(lat), declinacion_prueba, ang_h_punto)
            diferencia = abs(valor_punto - resultado_rejilla[i, j])
            max_diferencia = max(max_diferencia, diferencia)

    print(f"Diferencia maxima entre rejilla vectorizada y calculo punto a punto: {max_diferencia:.2e} radianes")
    if max_diferencia < 1e-10:
        print("OK: la rejilla vectorizada coincide con el modelo actual, celda a celda.")
    else:
        print("AVISO: hay una diferencia mayor de la esperada -- revisar antes de seguir.")

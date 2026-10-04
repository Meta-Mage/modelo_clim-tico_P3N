#GEOMETRÍA

from parametros import *
from orbita import *
import math

LATITUD_grados = 15
LATITUD_rad = LATITUD_grados * PI / 180

def declinacion_solar(momento, inclinacion_axial_rad):
    declinacion = math.asin(math.sin(inclinacion_axial_rad) * math.sin(momento + DESFASE_SOLSTICIO_RAD))
    return declinacion


# v2.4.2: el angulo horario ya no supone un dia de 24 h. "hora" son horas
# de 3600 s desde la medianoche del meridiano 0, y el sol da una vuelta
# completa (2*PI) en ROTACION_PERIODO (el dia SOLAR de P3N). Con
# ROTACION_PERIODO = 86400, MEDIO_DIA_H = 12.0 y la formula es exactamente
# la de antes, bit a bit: (hora - 12) * PI / 12.
MEDIO_DIA_H = ROTACION_PERIODO / 3600 / 2

def angulo_horario(hora, longitud_grados=0):
    angulo = (hora - MEDIO_DIA_H) * PI / MEDIO_DIA_H + longitud_grados * PI / 180
    return angulo

def angulo_cenital(latitud, declinacion, angulo):
    cos_cenital = math.sin(latitud) * math.sin(declinacion) + math.cos(latitud) * math.cos(declinacion) * math.cos(angulo)
    cenital_rad = math.acos(cos_cenital)
    return cenital_rad


if __name__ == "__main__":
    AV_demo = anomalia_verdadera(anomalia_excentrica(anomalia_media(DIA, HORA)))
    declinacion = declinacion_solar(AV_demo, INCLINACION_AXIAL_RAD)
    angulo = angulo_horario(HORA)
    angulo_cenital_solar = angulo_cenital(LATITUD_rad, declinacion, angulo)

    print(f"Declinación solar: {declinacion}")
    print(f"Ángulo horario: {angulo}")
    print(f"Ángulo cenital solar: {angulo_cenital_solar}")

#GEOMETRÍA

from parametros import *
from orbita import *
import math

LATITUD_grados = 15
LATITUD_rad = LATITUD_grados * PI / 180

def declinacion_solar(momento, inclinacion_axial_rad):
    declinacion = inclinacion_axial_rad * math.sin(momento + DESFASE_SOLSTICIO_RAD)
    return declinacion


def angulo_horario(hora, longitud_grados=0):
    angulo = (hora - 12) * PI / 12 + longitud_grados * PI / 180
    return angulo

def angulo_cenital(latitud, declinacion, angulo):
    cos_cenital = math.sin(latitud) * math.sin(declinacion) + math.cos(latitud) * math.cos(declinacion) * math.cos(angulo)
    cenital_rad = math.acos(cos_cenital)
    return cenital_rad


declinacion = declinacion_solar(AV, INCLINACION_AXIAL_RAD)
angulo = angulo_horario(HORA)
angulo_cenital_solar = angulo_cenital(LATITUD_rad, declinacion, angulo)

print(f"Declinación solar: {declinacion}")
print(f"Ángulo horario: {angulo}")
print(f"Ángulo cenital solar: {angulo_cenital_solar}")
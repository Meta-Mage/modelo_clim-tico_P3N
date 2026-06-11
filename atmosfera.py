#ATMÓSFERA
from parametros import *
from orbita import *
from geometria import * 
import math

def masa_aire(cenital_rad):
    cenital_grados = cenital_rad * 180 / PI
    altura_solar = 90 - cenital_grados
    if altura_solar <= 0:
        return None
    masa_aire_relativa = 1 / math.sin(math.radians(altura_solar + 244 / (165 + 47 * altura_solar ** 1.1)))
    return masa_aire_relativa

def trans(masa, profundidad_optica):
    transmitancia = math.exp(-masa * profundidad_optica)
    return transmitancia


masa_aire_relativa = masa_aire(angulo_cenital_solar)
transmitancia = trans(masa_aire_relativa, PROFUNDIDAD_OPTICA)

print(f"Masa de aire relativa: {masa_aire_relativa}")
print(f"Transmitancia: {transmitancia}")

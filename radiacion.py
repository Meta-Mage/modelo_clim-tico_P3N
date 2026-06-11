#RADIACIÓN
from parametros import *
from orbita import *
from geometria import * 
from atmosfera import *
import math

def i_toa(distancia, luminosidad):
    irradiancia_toa = luminosidad / (4*PI * distancia**2)
    return irradiancia_toa

def i_inst(irradiancia_toa, cenital_rad):
    irradiancia_instantanea = irradiancia_toa * math.cos(cenital_rad)
    return irradiancia_instantanea

def i_atm(irradiancia_instantanea, trans):
    irradiancia_atenuada = irradiancia_instantanea * trans
    return irradiancia_atenuada

def i_abs(irradiancia_atenuada, albedo):
    irradiancia_absorbida = irradiancia_atenuada * (1 - albedo)
    return irradiancia_absorbida

irradiancia_toa = i_toa(distancia, S3N_LUMINOSIDAD)
irradiancia_instantanea = i_inst( irradiancia_toa, angulo_cenital_solar)
irradiancia_atenuada = i_atm(irradiancia_instantanea, transmitancia)
irradiancia_absorbida = i_abs(irradiancia_atenuada, ALBEDO)

print(f"Irradiancia TOA: {irradiancia_toa}")
print(f"Irradiancia instantánea: {irradiancia_instantanea}")
print(f"Irradiancia con atenuación: {irradiancia_atenuada}")
print(f"Irradiancia absorbida: {irradiancia_absorbida}")

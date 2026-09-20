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


if __name__ == "__main__":
    AV_demo = anomalia_verdadera(anomalia_excentrica(anomalia_media(DIA, HORA)))
    distancia = distancia_S3N(AV_demo, SEMIEJE_MAYOR)
    declinacion_demo = declinacion_solar(AV_demo, INCLINACION_AXIAL_RAD)
    angulo_demo = angulo_horario(HORA)
    angulo_cenital_solar = angulo_cenital(LATITUD_rad, declinacion_demo, angulo_demo)
    masa_aire_relativa = masa_aire(angulo_cenital_solar)
    transmitancia = trans(masa_aire_relativa, PROFUNDIDAD_OPTICA)

    irradiancia_toa = i_toa(distancia, S3N_LUMINOSIDAD)
    irradiancia_instantanea = i_inst(irradiancia_toa, angulo_cenital_solar)
    irradiancia_atenuada = i_atm(irradiancia_instantanea, transmitancia)
    irradiancia_absorbida = i_abs(irradiancia_atenuada, ALBEDO)

    print(f"Irradiancia TOA: {irradiancia_toa}")
    print(f"Irradiancia instantánea: {irradiancia_instantanea}")
    print(f"Irradiancia con atenuación: {irradiancia_atenuada}")
    print(f"Irradiancia absorbida: {irradiancia_absorbida}")

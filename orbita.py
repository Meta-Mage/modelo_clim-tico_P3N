#ÓRBITA

from parametros import *
import math

DIA = 98
HORA = 16

def anomalia_media(dia, hora):
    tiempo_transcurrido = (dia -1) * 86400 + hora * 3600
    AM_grados = (19.043217 + tiempo_transcurrido * MOVIMIENTO_MEDIO) % 360
    AM_rad = AM_grados * PI/180
    return AM_rad

def anomalia_excentrica(AM_rad):
    AE_rad = AM_rad
    for i in range(100):
        AE_nuevo = AE_rad - (AE_rad - ORBITA_EXCENTRICIDAD * math.sin(AE_rad) - AM_rad) / (1 - ORBITA_EXCENTRICIDAD * math.cos(AE_rad))
        if abs(AE_nuevo - AE_rad) < 1e-10:
            break
        AE_rad = AE_nuevo
    return AE_rad

def anomalia_verdadera(AE_rad):
    factor = math.sqrt((1 + ORBITA_EXCENTRICIDAD) / (1 - ORBITA_EXCENTRICIDAD))
    AV_rad = 2 * math.atan(factor * math.tan(AE_rad / 2))
    return AV_rad

def distancia_S3N(AV_rad, semieje):
    distancia = semieje * (1 - ORBITA_EXCENTRICIDAD**2) / (1 + ORBITA_EXCENTRICIDAD * math.cos(AV_rad))
    return distancia


# Limites de las estaciones -- se calculan aqui, sin guardia, porque
# estacion() (justo debajo) los necesita como constantes globales, y
# otros modulos (temperatura.py, analisis_latitudes.py) llaman a
# estacion() e info_estaciones() directamente. Esto NO es codigo de
# prueba suelto -- es la base real de esas dos funciones.
AV_INVIERNO = anomalia_verdadera(anomalia_excentrica(anomalia_media(1, 0)))
AV_PRIMAVERA = (AV_INVIERNO + PI/2) % (2*PI)
AV_VERANO = (AV_INVIERNO + PI) % (2*PI)
AV_OTOÑO = (AV_INVIERNO + 3*PI/2) % (2*PI)

def estacion(momento):
    if momento >= AV_INVIERNO and momento < AV_PRIMAVERA:
        return "Invierno"
    elif momento >= AV_PRIMAVERA and momento < AV_VERANO:
        return "Primavera"
    elif momento >= AV_VERANO and momento < AV_OTOÑO:
        return "Verano"
    else:
        return "Otoño"

def info_estaciones():
    estaciones = ['Invierno', 'Primavera', 'Verano', 'Otoño']

    AM_inicio = anomalia_media(1, 0)
    t_perihelio_a_inicio = AM_inicio / (2 * PI) * ORBITA_PERIODO
    AV_inicio = anomalia_verdadera(anomalia_excentrica(AM_inicio))
    limites_AV = [(AV_inicio + i * PI/2) % (2*PI) for i in range(4)]

    def AV_a_dia(AV):
        factor = math.sqrt((1 - ORBITA_EXCENTRICIDAD) / (1 + ORBITA_EXCENTRICIDAD))
        AE = 2 * math.atan(factor * math.tan(AV / 2))
        if AE < 0:
            AE += 2 * PI
        AM = AE - ORBITA_EXCENTRICIDAD * math.sin(AE)
        if AM < 0:
            AM += 2 * PI
        t_desde_perihelio = AM / (2 * PI) * ORBITA_PERIODO
        t_desde_inicio = t_desde_perihelio - t_perihelio_a_inicio
        if t_desde_inicio < 0:
            t_desde_inicio += ORBITA_PERIODO
        return t_desde_inicio / 86400 + 1

    dias = [AV_a_dia(AV) for AV in limites_AV]

    resultados = []
    for i in range(4):
        dia_inicio = dias[i]
        dia_fin = dias[(i + 1) % 4]
        if dia_fin < dia_inicio:
            duracion = (ORBITA_PERIODO/86400 - dia_inicio + 1) + (dia_fin - 1)
        else:
            duracion = dia_fin - dia_inicio
        resultados.append({
            'estacion': estaciones[i],
            'AV_grados': math.degrees(limites_AV[i]),
            'dia_inicio': dia_inicio,
            'duracion': duracion
        })
    return resultados


if __name__ == "__main__":
    AM = anomalia_media(DIA, HORA)
    AE = anomalia_excentrica(AM)
    AV = anomalia_verdadera(AE)
    distancia = distancia_S3N(AV, SEMIEJE_MAYOR)

    print(f"Anomalía media: {AM}")
    print(f"Anomalía excéntrica: {AE}")
    print(f"Anomalía verdadera: {AV}")
    print(f"Distancia al Sol: {distancia}")
    print(f"Estación: {estacion(AV)}")

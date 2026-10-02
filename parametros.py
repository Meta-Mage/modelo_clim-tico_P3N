# CONSTANTES FÍSICAS UNIVERSALES

CONSTANTE_SB = 5.670374419e-8
CONSTANTE_GRAVITACIONAL = 6.6743e-11
PI = 3.141592653589793
MASA_SOL = 1_988_470_000_000_000_000_000_000_000_000
RADIO_SOL = 6.957e8
LUMINOSIDAD_SOL = 3.828e26
TEMPERATURA_SOL = 5778
 
# ESTRELLA S3N
S3N_MASAS_SOLARES = 0.884   # 02/10/2026, Fase 2b: calibrada para aire a 2 m medio global = 15.5 C (antes 0.97). Provisional hasta el barrido final de parametros.
EXP_MASA_RADIO = 0.8
EXP_MASA_TEMPERATURA = 0.475
S3N_MASA =  S3N_MASAS_SOLARES * MASA_SOL
S3N_RADIO = RADIO_SOL * (S3N_MASAS_SOLARES ** EXP_MASA_RADIO)
S3N_TEMPERATURA = TEMPERATURA_SOL * (S3N_MASAS_SOLARES ** EXP_MASA_TEMPERATURA)
S3N_LUMINOSIDAD = 4 * PI * S3N_RADIO**2 * CONSTANTE_SB * S3N_TEMPERATURA**4

# PLANETA P3N
P3N_CIRCUNFERENCIA = 29_384_382     
P3N_MASA = 2.9666081294151e24       

P3N_RADIO = P3N_CIRCUNFERENCIA / (2 * PI)
P3N_VOLUMEN = (4/3) * PI * P3N_RADIO ** 3
P3N_DENSIDAD = P3N_MASA / P3N_VOLUMEN
P3N_GRAVEDAD = CONSTANTE_GRAVITACIONAL * P3N_MASA / P3N_RADIO ** 2

#ÓRBITA
ORBITA_EXCENTRICIDAD = 0.046
ORBITA_PERIODO = 23_337_579.9978  
MOVIMIENTO_MEDIO = 360/ORBITA_PERIODO  
SEMIEJE_MAYOR = (CONSTANTE_GRAVITACIONAL * S3N_MASA * ORBITA_PERIODO**2 / (4 * PI**2)) ** (1/3)
# Anomalia media de P3N en el instante inicial del calendario (dia 1,
# 0 h), en grados: fija EN QUE PUNTO de la orbita empieza el año. Antes
# estaba escrita a mano dentro de orbita.py.
ANOMALIA_MEDIA_INICIO_GRADOS = 19.043217

# DESFASE DEL SOLSTICIO -- ahora se CALCULA. Antes era 242.1084 escrito a
# mano, calibrado para la excentricidad antigua (0.186); con la actual
# (0.046) dejaba el solsticio de invierno del hemisferio norte en el
# dia 5 en vez del dia 1 (corregido el 01/10/2026).
# Criterio del calendario de P3N: el dia 1 a las 0 h es exactamente el
# solsticio de invierno del hemisferio norte (declinacion minima). La
# declinacion es asin(sin(inclinacion) * sin(AV + DESFASE)), minima
# cuando AV + DESFASE = 270 grados, asi que DESFASE = 270 - AV(dia 1, 0 h).
# La ecuacion de Kepler se repite aqui (no se importa orbita.py porque
# orbita.py importa este archivo: seria una importacion circular).
def _anomalia_verdadera_inicio():
    import math
    am = math.radians(ANOMALIA_MEDIA_INICIO_GRADOS)
    ae = am
    for _ in range(100):
        ae_nuevo = ae - (ae - ORBITA_EXCENTRICIDAD * math.sin(ae) - am) / (1 - ORBITA_EXCENTRICIDAD * math.cos(ae))
        if abs(ae_nuevo - ae) < 1e-12:
            ae = ae_nuevo
            break
        ae = ae_nuevo
    factor = math.sqrt((1 + ORBITA_EXCENTRICIDAD) / (1 - ORBITA_EXCENTRICIDAD))
    return 2 * math.atan(factor * math.tan(ae / 2))

DESFASE_SOLSTICIO_RAD = (3 * PI / 2 - _anomalia_verdadera_inicio()) % (2 * PI)
ROTACION_PERIODO = 86400

#ATMÓSFERA Y SUPERFICIE
INCLINACION_AXIAL = 1.7
INCLINACION_AXIAL_RAD = INCLINACION_AXIAL * PI / 180
ALBEDO = 0.3
PROFUNDIDAD_OPTICA = 0.3
EMISIVIDAD = 0.77
INERCIA_TERMICA = 3500


# PRUEBA
if __name__ == "__main__":
    print(f"S3N. Masa: {S3N_MASA}; Radio: {S3N_RADIO}; Temperatura: {S3N_TEMPERATURA}; Luminosidad: {S3N_LUMINOSIDAD}.")
    print(f"P3N. Radio: {P3N_RADIO}; Volumen: {P3N_VOLUMEN}; Densidad: {P3N_DENSIDAD}; Gravedad: {P3N_GRAVEDAD}.")
    print(f"Semieje mayor: {SEMIEJE_MAYOR}")

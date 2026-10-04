# CONSTANTES FÍSICAS UNIVERSALES

CONSTANTE_SB = 5.670374419e-8
CONSTANTE_GRAVITACIONAL = 6.6743e-11
PI = 3.141592653589793
MASA_SOL = 1_988_470_000_000_000_000_000_000_000_000
RADIO_SOL = 6.957e8
LUMINOSIDAD_SOL = 3.828e26
TEMPERATURA_SOL = 5778
 
# ESTRELLA S3N
S3N_MASAS_SOLARES = 0.874   # 02/10/2026, v2.2c: calibrada para aire a 2 m medio global = 15.94 C (objetivo sorteado al azar entre 15.0 y 16.0, decision de Carlos). Con el hielo marino (v2.3) da 15.30 C, dentro del rango: Carlos decidio no recalibrar. Antes 0.97 (v2.2) y 0.884 (v2.2b). Provisional hasta el barrido final de parametros.
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
# Duracion del DIA SOLAR de P3N en segundos (de mediodia a mediodia).
# v2.4.3: 19,84 h (decision de Carlos y Ozan, 04/10/2026; antes 24 h
# provisionales). Es el dia SOLAR, el que viven los habitantes; el dia
# sideral (giro respecto a las estrellas), para un giro en el mismo
# sentido que la orbita, es 1 / (1/ROTACION_PERIODO + 1/ORBITA_PERIODO)
# = 71 206 s (19 h 46 min 46 s). Año: 326,75 dias solares.
# Requisito: multiplo exacto del paso de tiempo (temperatura.PASO_TIEMPO)
# -> lo comprueba la simulacion.
ROTACION_PERIODO = 71_424

# HORA DE P3N (v2.4.3, decision de Carlos 04/10/2026): 1/24 del dia solar,
# 2976 s (unos 49,6 min terrestres). El registro horario y la exportacion
# a H3N van en horas de P3N: 24 por dia.
HORAS_POR_DIA = 24
DURACION_HORA = ROTACION_PERIODO / HORAS_POR_DIA

# Rotacion y transporte de calor (v2.4.3): los coeficientes de difusion de
# la atmosfera estan calibrados en la Tierra; la difusividad de los
# remolinos escala como 1/Omega^2 (Williams y Kasting 1997), con Omega la
# velocidad de giro SIDERAL. FACTOR_ROTACION_D = (Omega_Tierra/Omega_P3N)^2.
OMEGA_TIERRA_SIDERAL = 7.2921159e-5      # rad/s
OMEGA_SIDERAL = 2 * PI * (1 / ROTACION_PERIODO + 1 / ORBITA_PERIODO)
FACTOR_ROTACION_D = (OMEGA_TIERRA_SIDERAL / OMEGA_SIDERAL) ** 2

#ATMÓSFERA Y SUPERFICIE
INCLINACION_AXIAL = 1.7
INCLINACION_AXIAL_RAD = INCLINACION_AXIAL * PI / 180
ALBEDO = 0.3
PROFUNDIDAD_OPTICA = 0.3
EMISIVIDAD = 0.77
INERCIA_TERMICA = 3500


# ================================================================
# MODO TIERRA (v3.0) -- SOLO PARA CALIBRAR Y VALIDAR
# ================================================================
# Con la variable de entorno M3N_MODO=tierra, el modelo simula la Tierra
# en vez de P3N: Sol, orbita, eje, dia y gravedad terrestres. Sirve para
# calibrar contra observaciones (los parametros etiquetados "(c)") y para
# comprobar que la fisica nueva da una Tierra razonable antes de fiarse
# de ella en P3N. Nunca se usa para el clima de P3N. El mapa de la Tierra
# lo da modo_tierra.py. Valores (todos de referencia, J2000):
#   - irradiancia a 1 UA: 1361 W/m2 (Kopp y Lean 2011);
#   - año sideral 365,256363 dias; excentricidad 0,0167086;
#   - longitud del perihelio 282,9373 grados (geocentrica) -> en el
#     solsticio de diciembre (longitud del Sol 270) la anomalia media es
#     347,486 grados: el dia 1 a las 0 h es el solsticio de invierno del
#     norte, mismo criterio que el calendario de P3N;
#   - oblicuidad 23,44 grados; dia solar 86 400 s; g = 9,80665 m/s2;
#     radio medio 6 371 km.
import os as _os
MODO_TIERRA = _os.environ.get("M3N_MODO", "").strip().lower() == "tierra"
if MODO_TIERRA:
    UNIDAD_ASTRONOMICA = 1.495978707e11
    IRRADIANCIA_TIERRA = 1361.0
    SEMIEJE_MAYOR = UNIDAD_ASTRONOMICA
    S3N_LUMINOSIDAD = IRRADIANCIA_TIERRA * 4 * PI * UNIDAD_ASTRONOMICA ** 2
    S3N_TEMPERATURA = 5772
    ORBITA_EXCENTRICIDAD = 0.0167086
    ORBITA_PERIODO = 365.256363 * 86400
    MOVIMIENTO_MEDIO = 360 / ORBITA_PERIODO
    ANOMALIA_MEDIA_INICIO_GRADOS = 347.486184
    DESFASE_SOLSTICIO_RAD = (3 * PI / 2 - _anomalia_verdadera_inicio()) % (2 * PI)
    ROTACION_PERIODO = 86_400
    DURACION_HORA = ROTACION_PERIODO / HORAS_POR_DIA
    OMEGA_SIDERAL = 2 * PI * (1 / ROTACION_PERIODO + 1 / ORBITA_PERIODO)
    FACTOR_ROTACION_D = (OMEGA_TIERRA_SIDERAL / OMEGA_SIDERAL) ** 2     # = 1 (por construccion)
    P3N_RADIO = 6.371e6
    P3N_GRAVEDAD = 9.80665
    INCLINACION_AXIAL = 23.44
    INCLINACION_AXIAL_RAD = INCLINACION_AXIAL * PI / 180


# PRUEBA
if __name__ == "__main__":
    print(f"S3N. Masa: {S3N_MASA}; Radio: {S3N_RADIO}; Temperatura: {S3N_TEMPERATURA}; Luminosidad: {S3N_LUMINOSIDAD}.")
    print(f"P3N. Radio: {P3N_RADIO}; Volumen: {P3N_VOLUMEN}; Densidad: {P3N_DENSIDAD}; Gravedad: {P3N_GRAVEDAD}.")
    print(f"Semieje mayor: {SEMIEJE_MAYOR}")

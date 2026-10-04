# prueba_paso.py -- prueba de sensibilidad del paso de tiempo (v2.4.3).
#
# Repite la simulacion del mapa activo con el paso normal (temperatura.PASO_TIEMPO,
# 992 s) y con la MITAD (496 s), y compara el aire a 2 m medio de cada latitud.
# Si la fisica esta bien resuelta, las diferencias deben ser pequeñas
# (criterio: < 0,1 C de media global y < 0,3 C en cualquier banda de latitud).
# Es la misma prueba que se hizo en la Fase 0 con 900 y 450 s.
#
# Uso:  python prueba_paso.py      (la simulacion con 496 s tarda el doble; ~15 min)

import numpy as np

import temperatura
from parametros import S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR, EMISIVIDAD, PROFUNDIDAD_OPTICA
from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from cache_simulacion import precalcular_orbita_cacheada, simular_fase2b_cacheada
from puente_c3n import cargar_mapa_activo_de_c3n
from rejilla import LATITUDES_GRADOS


def simular(paso):
    original = temperatura.PASO_TIEMPO
    temperatura.PASO_TIEMPO = paso
    try:
        orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
        return simular_fase2b_cacheada(orbita, tipo, altitud, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO,
                                       PROFUNDIDAD_OPTICA, float(D_DIFUSION_REFERENCIA),
                                       nombre_mapa=f"{nombre}_paso{paso}", paso_tiempo=paso)
    finally:
        temperatura.PASO_TIEMPO = original


if __name__ == "__main__":
    tipo, altitud, nombre = cargar_mapa_activo_de_c3n()
    paso = temperatura.PASO_TIEMPO
    print(f"Mapa '{nombre}': simulando con {paso} s y con {paso // 2} s...")
    a = simular(paso)
    b = simular(paso // 2)
    peso = np.cos(np.radians(LATITUDES_GRADOS))
    za = a["reg_media"].mean(axis=(0, 2)); zb = b["reg_media"].mean(axis=(0, 2))
    ga = np.average(za, weights=peso); gb = np.average(zb, weights=peso)
    print(f"Aire a 2 m medio global: {ga:.3f} C ({paso} s)  frente a  {gb:.3f} C ({paso // 2} s)  -> diferencia {ga - gb:+.3f} C")
    print(f"Mayor diferencia en una banda de latitud: {np.max(np.abs(za - zb)):.3f} C")
    ok = abs(ga - gb) < 0.1 and np.max(np.abs(za - zb)) < 0.3
    print("RESULTADO:", "SUPERADA (el paso no influye de forma apreciable)" if ok else "NO SUPERADA: pegalo en el chat")

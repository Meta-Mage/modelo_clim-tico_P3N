import numpy as np
from parametros import *
from temperatura import precalcular_orbita
from rejilla import FILAS, COLUMNAS
from fase1_geografia import (
    TIERRA, AGUA, GRADIENTE_TERMICO,
    mapa_falso_todo_tierra,
    simular_rejilla_geografia,
)

# ============================================================================
# AVISO IMPORTANTE (PROVISIONAL):
# GRADIENTE_TERMICO = 6.5 C/km es el "lapse rate" ambiental estandar de la
# atmosfera TERRESTRE (valor empirico real de la Tierra), aplicado aqui como
# aproximacion provisional porque M3N todavia no modela estructura atmosferica
# vertical (capas, humedad, presion). NO esta derivado de la fisica propia de
# P3N. Debe revisarse cuando exista un modelo atmosferico multicapa.
# Ver NOTAS_DISENO_FASE0.md.
# ============================================================================

ALBEDO_PRUEBA = {TIERRA: 0.3, AGUA: 0.1}
INERCIA_PRUEBA = {TIERRA: 3500, AGUA: 1000}

ALTITUD_PRUEBA_METROS = 2000
FILA_PRUEBA, COLUMNA_PRUEBA = 10, 20  # celda cualquiera, lejos de los polos


def mapa_con_una_celda_alta():
    tipo_superficie, altitud_metros = mapa_falso_todo_tierra()
    altitud_metros = altitud_metros.copy()
    altitud_metros[FILA_PRUEBA, COLUMNA_PRUEBA] = ALTITUD_PRUEBA_METROS
    return tipo_superficie, altitud_metros


if __name__ == "__main__":
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    tipo_base, altitud_base = mapa_falso_todo_tierra()
    T_base, anos_base = simular_rejilla_geografia(
        datos_orbita, tipo_base, altitud_base, EMISIVIDAD,
        ALBEDO_PRUEBA, INERCIA_PRUEBA, PROFUNDIDAD_OPTICA,
    )

    tipo_alta, altitud_alta = mapa_con_una_celda_alta()
    T_alta, anos_alta = simular_rejilla_geografia(
        datos_orbita, tipo_alta, altitud_alta, EMISIVIDAD,
        ALBEDO_PRUEBA, INERCIA_PRUEBA, PROFUNDIDAD_OPTICA,
    )

    print(f"Convergencia mapa base (altitud 0 en todas las celdas): {anos_base} ano(s)")
    print(f"Convergencia mapa con una celda a {ALTITUD_PRUEBA_METROS} m: {anos_alta} ano(s)")

    mascara_sin_cambio = np.ones((FILAS, COLUMNAS), dtype=bool)
    mascara_sin_cambio[FILA_PRUEBA, COLUMNA_PRUEBA] = False
    diferencia_resto = np.max(np.abs(T_alta[mascara_sin_cambio] - T_base[mascara_sin_cambio]))
    print(f"Diferencia maxima en celdas NO modificadas: {diferencia_resto:.2e} C")

    diferencia_celda = T_base[FILA_PRUEBA, COLUMNA_PRUEBA] - T_alta[FILA_PRUEBA, COLUMNA_PRUEBA]
    esperado = GRADIENTE_TERMICO * (ALTITUD_PRUEBA_METROS / 1000)
    print(f"Celda de prueba: nivel del mar = {T_base[FILA_PRUEBA, COLUMNA_PRUEBA]:.4f} C | "
          f"a {ALTITUD_PRUEBA_METROS} m = {T_alta[FILA_PRUEBA, COLUMNA_PRUEBA]:.4f} C | "
          f"diferencia = {diferencia_celda:.4f} C (esperado {esperado:.4f} C)")

    tolerancia = 1e-9
    ok_resto = diferencia_resto < tolerancia
    ok_celda = abs(diferencia_celda - esperado) < tolerancia

    if ok_resto and ok_celda:
        print("OK: la correccion por altitud actua exactamente como se espera.")
        print("RECORDATORIO: el gradiente de 6.5 C/km es un valor terrestre PROVISIONAL,")
        print("pendiente de revision cuando exista un modelo atmosferico multicapa.")
    else:
        print("AVISO: revisar antes de seguir.")

# fase2_combinado.py -- Fase 2: transporte horizontal (difusion,
# fase2_difusion.py) + conduccion vertical por columna (fase2_inercia_
# multicapa.py), COMBINADOS en el mismo paso de tiempo.
#
# ================================================================
# POR QUE SE COMBINAN (analisis resumido, la version completa esta en
# la conversacion, no repetida aqui para no duplicar texto largo)
# ================================================================
# Las dos fisicas actuan sobre grados de libertad distintos y no
# compiten entre si:
#   - Difusion horizontal: transporte de calor ENTRE columnas vecinas,
#     representando circulacion atmosferica/oceanica (no conduccion
#     literal por el terreno).
#   - Conduccion vertical: transporte DENTRO de una misma columna,
#     hacia el subsuelo (conduccion solida real).
#
# La conduccion horizontal directa entre columnas de roca (calor
# viajando por el subsuelo de una celda a la vecina) es fisicamente
# despreciable a las distancias de la rejilla (cientos de km): el
# tiempo caracteristico de esa conduccion, distancia^2/difusividad,
# da ~10^16 segundos con nuestra difusividad (K=1e-6 m2/s) -- cientos
# de millones de años. Por eso, igual que en los esquemas de
# superficie terrestre reales (ISBA, Noah, CLM), la difusion
# horizontal se aplica SOLO a la capa expuesta a la atmosfera (la
# capa 0 de la columna), nunca a las capas profundas -- estas solo
# intercambian calor verticalmente, dentro de su propia columna.
# ================================================================
#
# METODO NUMERICO: "operator splitting" con TRES sub-pasos por paso de
# tiempo (en vez de los dos habituales en el resto de Fase 2):
#   1. Radiativo, explicito -- solo en la capa 0 (como siempre).
#   2. Difusion horizontal, implicita -- solo en la capa 0 (matriz
#      dispersa de fase2_difusion.py, pero construida con la
#      capacidad de la capa 0 de la columna, NO la capacidad de una
#      unica capa como en fase2_difusion.py puro -- ver mas abajo).
#   3. Conduccion vertical, implicita -- toda la columna, celda por
#      celda de la rejilla (fase2_inercia_multicapa.py).
#
# Cada sub-paso es incondicionalmente estable (todos implicitos, salvo
# el radiativo). El orden (difusion horizontal antes que conduccion
# vertical) es una eleccion razonable pero no la unica posible -- con
# un paso de tiempo pequeño frente a las escalas fisicas relevantes
# (que es el caso aqui, 900 s), el error de "particionar" el problema
# en sub-pasos en vez de resolverlo como un unico sistema 3D es
# pequeño, del mismo tipo que ya se acepta en el resto de Fase 2.
#
# DETALLE TECNICO IMPORTANTE: la matriz de difusion horizontal necesita
# la capacidad calorifica de la capa que esta difundiendo (aparece en
# la ecuacion como C/paso_tiempo). Antes (fase2_difusion.py puro) esa
# capacidad era la de la UNICA capa del modelo (inercia * sqrt(dia/pi),
# ~415.000 J/m2/K en tierra). Ahora que existe una columna, la capa 0
# tiene mucha menos masa (solo el primer centimetro y pico, ~104.000
# J/m2/K en tierra con 8 capas) -- es FISICAMENTE coherente (esa capa
# 0 representa de verdad solo una fina piel superficial, no toda la
# profundidad diurna), pero significa que la superficie reacciona mas
# rapido tanto a la radiacion como a la difusion horizontal que en la
# version anterior. Este archivo usa la capacidad de la capa 0 de la
# columna en la matriz de difusion, correctamente.

import math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from parametros import ROTACION_PERIODO, PI, CONSTANTE_SB
from temperatura import precalcular_orbita, PASO_TIEMPO as PASO_TIEMPO_POR_DEFECTO
from rejilla import FILAS, COLUMNAS
from fase1_geografia import (
    TIERRA, AGUA, correccion_altitud, estimar_T_inicial_equilibrio,
    irradiancia_absorbida_desde_toa_rejilla_con_albedo,
)
from fase2_difusion import construir_matriz_difusion, D_DIFUSION_REFERENCIA
from fase2_inercia_multicapa import (
    construir_columna_rejilla, paso_conduccion_implicito,
    K_DIFUSIVIDAD_TIERRA, N_CAPAS_DEFECTO,
)


def simular_rejilla_combinada_con_registro(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
    K_difusividad=K_DIFUSIVIDAD_TIERRA, n_capas=N_CAPAS_DEFECTO,
    paso_tiempo=PASO_TIEMPO_POR_DEFECTO, max_anos=50, tolerancia_convergencia=0.015,
    # 0.015 en vez de 0.01 (eleccion de Carlos, 21/09): relajar "muy
    # poco" la tolerancia para converger en algun año menos, sin
    # renunciar a una precision fina (sigue siendo un cambio de 0.015 C
    # entre un año y el siguiente, imperceptible en la practica).
    verificar_energia=False,
):
    """
    Version de rejilla completa con registro diario (minima/media/
    maxima) de la CAPA SUPERFICIAL -- misma forma de devolver el
    resultado que el resto de Fase 2 (consulta_punto.py, mapa_calor.py,
    analisis_latitudes.py pueden usarla sin cambiar como leen el
    resultado). Combina difusion horizontal (solo capa 0) y conduccion
    vertical por columna (toda la columna) -- ver cabecera del archivo.

    verificar_energia=True añade un SEXTO valor de retorno: un dict con
    la comprobacion de conservacion de energia global pedida como
    criterio de validacion de Fase 2 en la tabla de fases ("suma
    absorbida = suma emitida en equilibrio") -- acumulada, ponderada por
    area real (cos(latitud)), durante el año final que ya se simula de
    todas formas para construir el registro. En un sistema en equilibrio
    periodico verdadero, la energia total al principio y al final de ese
    año coincide, asi que la suma de (absorbido - emitido) de todo el
    año, en todas las celdas, debe ser ~0 comparada con la escala de lo
    absorbido o emitido (la difusion y la conduccion SOLO redistribuyen
    energia, no la crean ni la destruyen).
    """
    albedo_grid = np.where(tipo_superficie == TIERRA, albedo_por_tipo[TIERRA], albedo_por_tipo[AGUA])

    capacidades, conductancias = construir_columna_rejilla(tipo_superficie, inercia_por_tipo, K_difusividad, n_capas)
    C0_flat = capacidades[..., 0].flatten()

    tiene_difusion = np.max(D_grid) > 0
    if tiene_difusion:
        L = construir_matriz_difusion(D_grid)
        A = sp.diags(C0_flat / paso_tiempo) - L
        A_factorizada = splu(A.tocsc())
    else:
        A_factorizada = None

    def paso_fisica(T_columna, abs_local):
        # ---- 1. radiativo, explicito, solo capa 0 ----
        emitido = (1 - emisividad / 2) * CONSTANTE_SB * T_columna[..., 0] ** 4
        T_columna = T_columna.copy()
        T_columna[..., 0] = T_columna[..., 0] + (paso_tiempo / capacidades[..., 0]) * (abs_local - emitido)

        # ---- 2. difusion horizontal, implicita, solo capa 0 ----
        if tiene_difusion:
            lado_derecho = (C0_flat / paso_tiempo) * T_columna[..., 0].flatten()
            T_columna[..., 0] = A_factorizada.solve(lado_derecho).reshape(FILAS, COLUMNAS)

        # ---- 3. conduccion vertical, implicita, toda la columna ----
        T_columna = paso_conduccion_implicito(T_columna, capacidades, conductancias, paso_tiempo)

        return T_columna

    T_inicial = estimar_T_inicial_equilibrio(datos_orbita, albedo_grid, emisividad, profundidad_optica)
    T_columna = np.zeros((FILAS, COLUMNAS, n_capas))
    for i in range(n_capas):
        T_columna[..., i] = T_inicial

    for ano in range(max_anos):
        T_inicio_ano = T_columna[..., 0].copy()
        for toa, decl, ang_h_lon0 in datos_orbita:
            abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
                toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
            )
            T_columna = paso_fisica(T_columna, abs_local)
        diferencia_maxima = np.max(np.abs(T_columna[..., 0] - T_inicio_ano))
        if diferencia_maxima < tolerancia_convergencia:
            anos_convergencia = ano + 1
            break
    else:
        anos_convergencia = max_anos

    pasos_por_dia = round(ROTACION_PERIODO / paso_tiempo)

    registro_minima, registro_media, registro_maxima = [], [], []
    T_min_dia = T_max_dia = suma_dia = None
    contador_dia = 0

    if verificar_energia:
        peso_area = np.cos(np.radians(90 - (np.arange(FILAS) + 0.5) * (180 / FILAS))).reshape(-1, 1)
        suma_absorbido_ponderado = 0.0
        suma_emitido_ponderado = 0.0

    def guardar_dia():
        T_media_dia = suma_dia / contador_dia
        registro_minima.append(correccion_altitud(T_min_dia - 273.15, altitud_metros))
        registro_media.append(correccion_altitud(T_media_dia - 273.15, altitud_metros))
        registro_maxima.append(correccion_altitud(T_max_dia - 273.15, altitud_metros))

    for paso, (toa, decl, ang_h_lon0) in enumerate(datos_orbita):
        if paso % pasos_por_dia == 0:
            if contador_dia > 0:
                guardar_dia()
            T_sup = T_columna[..., 0]
            T_min_dia = np.full_like(T_sup, np.inf)
            T_max_dia = np.full_like(T_sup, -np.inf)
            suma_dia = np.zeros_like(T_sup)
            contador_dia = 0

        abs_local = irradiancia_absorbida_desde_toa_rejilla_con_albedo(
            toa, decl, ang_h_lon0, albedo_grid, profundidad_optica
        )

        if verificar_energia:
            emitido_local = (1 - emisividad / 2) * CONSTANTE_SB * T_columna[..., 0] ** 4
            suma_absorbido_ponderado += np.sum(abs_local * peso_area)
            suma_emitido_ponderado += np.sum(emitido_local * peso_area)

        T_columna = paso_fisica(T_columna, abs_local)
        T_sup = T_columna[..., 0]

        T_min_dia = np.minimum(T_min_dia, T_sup)
        T_max_dia = np.maximum(T_max_dia, T_sup)
        suma_dia = suma_dia + T_sup
        contador_dia += 1

    if contador_dia == pasos_por_dia:
        guardar_dia()

    registro_minima = np.array(registro_minima)
    registro_media = np.array(registro_media)
    registro_maxima = np.array(registro_maxima)

    T_final_corregida = correccion_altitud(T_columna[..., 0] - 273.15, altitud_metros)

    if verificar_energia:
        diferencia_relativa = abs(suma_absorbido_ponderado - suma_emitido_ponderado) / suma_absorbido_ponderado
        diagnostico_energia = {
            "suma_absorbido": suma_absorbido_ponderado,
            "suma_emitido": suma_emitido_ponderado,
            "diferencia_relativa": diferencia_relativa,
        }
        return T_final_corregida, anos_convergencia, registro_minima, registro_media, registro_maxima, diagnostico_energia

    return T_final_corregida, anos_convergencia, registro_minima, registro_media, registro_maxima


if __name__ == "__main__":
    import time
    from parametros import EMISIVIDAD, PROFUNDIDAD_OPTICA, INCLINACION_AXIAL_RAD, S3N_LUMINOSIDAD, SEMIEJE_MAYOR
    from fase1_geografia import ALBEDO_POR_TIPO, INERCIA_POR_TIPO
    from fase1_mapa_mixto import mapa_mitad_agua_mitad_tierra
    from fase2_inercia_multicapa import simular_rejilla_multicapa_con_registro
    from fase2_difusion import simular_rejilla_difusion_con_registro

    print("Precalculando orbita...")
    datos_orbita = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    tipo_mixto, altitud_mixta = mapa_mitad_agua_mitad_tierra()

    # ---- Prueba 0: D=0 debe coincidir EXACTO con la version multicapa pura (sin difusion) ----
    print("\n--- Prueba 0: D=0 debe coincidir EXACTO con fase2_inercia_multicapa.py (sin difusion) ---")
    D_cero = np.zeros((FILAS, COLUMNAS))
    _, anos_multi, reg_min_multi, reg_media_multi, reg_max_multi = simular_rejilla_multicapa_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
    )
    _, anos_comb0, reg_min_comb0, reg_media_comb0, reg_max_comb0 = simular_rejilla_combinada_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
        D_cero,
    )
    diferencia0 = max(
        np.max(np.abs(reg_min_multi - reg_min_comb0)),
        np.max(np.abs(reg_media_multi - reg_media_comb0)),
        np.max(np.abs(reg_max_multi - reg_max_comb0)),
    )
    print(f"Convergencia: {anos_multi} (multicapa pura) vs {anos_comb0} (combinado, D=0)")
    print(f"Diferencia maxima en el registro: {diferencia0:.2e} C")
    print("OK: coincide EXACTA." if diferencia0 < 1e-6 else "AVISO: deberian coincidir exactas -- revisar.")

    # ---- Prueba 1 (RETIRADA): se intento comparar, celda a celda de
    # agua, el combinado con n_capas=1 contra fase2_difusion.py puro,
    # esperando una coincidencia EXACTA. No es una comparacion valida,
    # y por partida doble:
    #   1. Con n_capas=1, la "columna" de TIERRA no reduce a su
    #      capacidad diurna de siempre -- pasa a cubrir toda la
    #      profundidad de la columna (metros, no cm), algo distinto a
    #      proposito (ya se documento igual en fase2_inercia_multicapa.py).
    #   2. Mas importante: la difusion horizontal es un UNICO sistema
    #      que acopla TODAS las celdas vecinas a la vez, agua y tierra
    #      mezcladas -- no se puede aislar el comportamiento del agua
    #      sin que le afecte la capacidad de sus vecinas de tierra. Como
    #      la tierra en el modelo combinado tiene una capa superficial
    #      mas fina (solo la "piel", no toda la profundidad diurna),
    #      eso cambia ligeramente el resultado tambien en las celdas de
    #      agua vecinas a tierra -- correctamente, es el acoplamiento
    #      haciendo su trabajo, no un error. Una comparacion "exacta"
    #      solo tendria sentido en un mapa sin ninguna celda de tierra,
    #      caso trivial que no aporta mas confianza que la Prueba 0.
    # Se deja esta nota en vez de la prueba (que en un intento real dio
    # una diferencia de 1.18 C, pequeña y razonable para un efecto de
    # frontera, no un fallo) para no repetir una simulacion de rejilla
    # completa (~10 min) por una comparacion que no es correcta de partida.

    # ---- Prueba 2: combinado de verdad (D=0.55 escalado a P3N + n_capas
    # por defecto, ahora 6 en vez de 8 -- 21/09, eleccion de Carlos para
    # optimizar tiempo de calculo) -- comportamiento y tiempo ----
    print(f"\n--- Prueba 2: combinado real (difusion + {N_CAPAS_DEFECTO} capas, por defecto) -- rango anual y tiempo ---")
    D_referencia = np.full((FILAS, COLUMNAS), D_DIFUSION_REFERENCIA)
    _, anos_dif, reg_min_dif, reg_media_dif, reg_max_dif = simular_rejilla_difusion_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
        D_referencia,
    )
    inicio = time.time()
    T_comb, anos_comb, reg_min_comb, reg_media_comb, reg_max_comb, diag_energia = simular_rejilla_combinada_con_registro(
        datos_orbita, tipo_mixto, altitud_mixta, EMISIVIDAD, ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA,
        D_referencia, verificar_energia=True,
    )
    duracion = time.time() - inicio
    print(f"Convergencia: {anos_comb} año(s) en {duracion:.1f} s")
    print(f"Rango anual real (minimo/maximo de todo el año, con ciclo diurno+estacional): "
          f"{reg_min_comb.min():.1f} C a {reg_max_comb.max():.1f} C")
    print(f"(Para referencia -- solo difusion, sin columna: {reg_min_dif.min():.1f} C a {reg_max_dif.max():.1f} C | "
          f"solo columna, sin difusion: {reg_min_multi.min():.1f} C a {reg_max_multi.max():.1f} C)")

    # ---- Prueba 3: conservacion de energia global (criterio de
    # validacion de Fase 2 segun la tabla de fases: "suma absorbida =
    # suma emitida en equilibrio"). Reusa el diagnostico calculado
    # durante la Prueba 2 (mismo año final ya simulado), evitando otra
    # simulacion de rejilla completa (~10 min) solo para esto.
    print("\n--- Prueba 3: conservacion de energia global (año final, ponderado por area) ---")
    print(f"Suma absorbida (ponderada): {diag_energia['suma_absorbido']:.6e}")
    print(f"Suma emitida (ponderada):   {diag_energia['suma_emitido']:.6e}")
    print(f"Diferencia relativa: {diag_energia['diferencia_relativa']:.2e}")
    print("OK: conserva energia (diferencia < 1%)." if diag_energia["diferencia_relativa"] < 0.01
          else "AVISO: diferencia mayor de lo esperado -- revisar.")

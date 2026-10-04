# cache_simulacion.py -- cache de resultados costosos (orbita, simulacion
# de rejilla), invalidada AUTOMATICAMENTE si cambia cualquier parametro de
# entrada. Peticion de Carlos (21/09): "hay cosas como los calculos
# orbitales que no deberian cambiar nunca, pero es importante que si
# cambiamos algun parametro se pueda detectar para que reaccione bien".
#
# ================================================================
# COMO FUNCIONA
# ================================================================
# Cada resultado costoso se guarda en disco con un nombre de archivo que
# es un hash (SHA-256, truncado) de TODOS los valores de entrada que
# afectan a ese resultado -- constantes fisicas, el mapa activo, n_capas,
# la tolerancia de convergencia, etc. Si vuelves a pedir el mismo
# resultado con los MISMOS parametros, el hash coincide, se encuentra el
# archivo, y se carga en vez de recalcular. Si cambias CUALQUIER
# parametro (una constante en parametros.py, el mapa activo en C3N,
# n_capas, la tolerancia...), el hash cambia automaticamente, no se
# encuentra archivo, y se recalcula -- no hace falta borrar nada a mano
# ni acordarse de invalidar la cache. Es un mecanismo generico: sirve
# igual para el precalculo de la orbita (barato, pero se beneficia igual)
# y para la simulacion completa de rejilla (el coste real, minutos).
#
# Importante: esto NO es "inteligente" -- no sabe que dos conjuntos de
# parametros parecidos dan resultados parecidos, es una comparacion
# EXACTA. Un cambio de una coma flotante en la ultima cifra decimal
# produce un hash distinto y fuerza un recalculo. Es intencionado:
# preferible recalcular de mas (unos minutos perdidos) a servir sin
# darse cuenta un resultado que ya no corresponde a los parametros
# actuales.
#
# La cache vive en outputs/cache/ (se crea sola). Si en algun momento
# quieres forzar que TODO se recalcule desde cero (por ejemplo, para
# comprobar que no hay ningun resultado obsoleto sirviendose por error),
# basta con borrar esa carpeta entera -- no rompe nada, simplemente se
# volvera a rellenar la primera vez que haga falta cada resultado.

import hashlib
import os
import pickle

import numpy as np

# v2.4.1: ruta absoluta (junto a este archivo), para que la cache sea la
# misma aunque M3N se use desde otra carpeta (por ejemplo, desde H3N).
CARPETA_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "cache")


def _alimentar_hash(hasher, obj):
    if isinstance(obj, np.ndarray):
        hasher.update(b"ndarray")
        hasher.update(str(obj.shape).encode())
        hasher.update(str(obj.dtype).encode())
        hasher.update(np.ascontiguousarray(obj).tobytes())
    elif isinstance(obj, (list, tuple)):
        hasher.update(f"seq{len(obj)}".encode())
        for elemento in obj:
            _alimentar_hash(hasher, elemento)
    elif isinstance(obj, dict):
        hasher.update(f"dict{len(obj)}".encode())
        for clave in sorted(obj.keys()):
            hasher.update(str(clave).encode())
            _alimentar_hash(hasher, obj[clave])
    elif isinstance(obj, float):
        # repr() da la representacion EXACTA de un float de Python (no
        # redondeada como str()), asi que un cambio real en la ultima
        # cifra decimal de un parametro se detecta siempre.
        hasher.update(repr(obj).encode())
    else:
        hasher.update(repr(obj).encode())


def _hash_de_objeto(obj):
    hasher = hashlib.sha256()
    _alimentar_hash(hasher, obj)
    return hasher.hexdigest()[:24]  # 24 caracteres hex bastan para este uso (nada criptografico)


def cargar_o_calcular(nombre_funcion, parametros_clave, funcion_calculo, etiqueta=""):
    """
    nombre_funcion: identifica QUE se cachea (p.ej. "orbita",
        "simulacion_combinada") -- separa archivos de cosas distintas.
    parametros_clave: lista/tupla con TODOS los valores de entrada que
        afectan al resultado, y SOLO esos (incluir algo que no influye en
        el resultado invalidaria la cache sin necesidad cada vez que
        cambie, aunque el resultado real fuera identico).
    funcion_calculo: funcion sin argumentos que calcula el resultado real
        si no esta en cache (normalmente un lambda que cierra sobre los
        argumentos reales).
    etiqueta: texto opcional SOLO para que el nombre de archivo sea mas
        legible para un humano mirando la carpeta (p.ej. el nombre del
        mapa activo) -- no participa en el hash, cambiarla no invalida
        nada.
    """
    os.makedirs(CARPETA_CACHE, exist_ok=True)
    hash_corto = _hash_de_objeto(parametros_clave)
    sufijo_legible = f"_{etiqueta}" if etiqueta else ""
    nombre_archivo = f"{nombre_funcion}{sufijo_legible}_{hash_corto}.pkl"
    ruta = os.path.join(CARPETA_CACHE, nombre_archivo)

    if os.path.exists(ruta):
        with open(ruta, "rb") as f:
            resultado = pickle.load(f)
        print(f"[cache] '{nombre_funcion}': encontrado en cache ({nombre_archivo}), no se recalcula.")
        return resultado

    print(f"[cache] '{nombre_funcion}': no encontrado (parametros nuevos o cambiados) -- calculando...")
    resultado = funcion_calculo()

    with open(ruta, "wb") as f:
        pickle.dump(resultado, f)
    print(f"[cache] '{nombre_funcion}': guardado para la proxima vez ({nombre_archivo}).")
    return resultado


# ================================================================
# ENVOLTORIOS CONCRETOS -- para precalcular_orbita(), para el modelo de
# la Fase 2 (simular_rejilla_combinada_con_registro, se conserva como
# referencia), para el modelo actual (simular_fase2b) y el punto de
# entrada de las herramientas (simular_modelo_cacheado). Mismos
# argumentos y mismo resultado que las funciones reales; la unica
# diferencia es que pasan primero por la cache.
# ================================================================

def precalcular_orbita_cacheada(luminosidad, inclinacion_axial_rad, semieje):
    from parametros import (
        ORBITA_PERIODO, ORBITA_EXCENTRICIDAD,
        ANOMALIA_MEDIA_INICIO_GRADOS, DESFASE_SOLSTICIO_RAD, ROTACION_PERIODO,
    )
    from temperatura import precalcular_orbita, PASO_TIEMPO

    # luminosidad, inclinacion_axial_rad y semieje ya llegan como
    # argumentos, pero el resto se lee como globales dentro de
    # precalcular_orbita() (y de anomalia_media, anomalia_excentrica,
    # declinacion_solar...) -- hay que incluirlos aqui a mano o un cambio
    # en cualquiera de ellos no se detectaria. v2 (01/10/2026): se añaden
    # excentricidad, anomalia media inicial y desfase del solsticio, que
    # faltaban en v1.
    # v3 (v2.4.3): CORRIGE que la clave no incluyera la duracion del dia
    # (ROTACION_PERIODO): al cambiarla sin cambiar el paso de tiempo se habria
    # reutilizado una orbita con el sol en la posicion del dia antiguo. Se
    # añade tambien la huella del codigo de la fisica, como en la simulacion.
    clave = (
        "precalcular_orbita_v3", luminosidad, inclinacion_axial_rad, semieje,
        ORBITA_PERIODO, PASO_TIEMPO, ROTACION_PERIODO,
        ORBITA_EXCENTRICIDAD, ANOMALIA_MEDIA_INICIO_GRADOS, DESFASE_SOLSTICIO_RAD,
        huella_codigo(),
    )
    return cargar_o_calcular(
        "orbita", clave,
        lambda: precalcular_orbita(luminosidad, inclinacion_axial_rad, semieje),
    )


def simular_rejilla_combinada_cacheada(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
    K_difusividad=None, n_capas=None, paso_tiempo=None, max_anos=50,
    tolerancia_convergencia=None, verificar_energia=False, nombre_mapa="",
    acelerar_convergencia=True,
):
    from fase2_combinado import simular_rejilla_combinada_con_registro
    from fase2_inercia_multicapa import K_DIFUSIVIDAD_TIERRA, N_CAPAS_DEFECTO
    from temperatura import PASO_TIEMPO as PASO_TIEMPO_DEFECTO

    # Los valores por defecto se resuelven aqui explicitamente (en vez de
    # dejarlos como None en la clave) para que, si el valor por defecto
    # de la funcion real cambia en el futuro, ese cambio SI invalide la
    # cache -- si se dejara "None" en la clave, dos llamadas con defaults
    # distintos compartirian hash por error.
    if K_difusividad is None:
        K_difusividad = K_DIFUSIVIDAD_TIERRA
    if n_capas is None:
        n_capas = N_CAPAS_DEFECTO
    if paso_tiempo is None:
        paso_tiempo = PASO_TIEMPO_DEFECTO
    if tolerancia_convergencia is None:
        tolerancia_convergencia = 0.015

    albedo_items = tuple(sorted(albedo_por_tipo.items()))
    inercia_items = tuple(sorted(inercia_por_tipo.items()))

    clave = (
        # v2 (02/10/2026): añade acelerar_convergencia (arranque cercano
        # al equilibrio, ver fase2_combinado.py), que cambia ligeramente
        # el resultado (centesimas de grado) frente a v1.
        "simular_combinada_v2", datos_orbita, tipo_superficie, altitud_metros,
        emisividad, albedo_items, inercia_items, profundidad_optica, D_grid,
        K_difusividad, n_capas, paso_tiempo, max_anos, tolerancia_convergencia,
        verificar_energia, acelerar_convergencia,
    )
    return cargar_o_calcular(
        "simulacion_combinada", clave,
        lambda: simular_rejilla_combinada_con_registro(
            datos_orbita, tipo_superficie, altitud_metros, emisividad,
            albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
            K_difusividad=K_difusividad, n_capas=n_capas, paso_tiempo=paso_tiempo,
            max_anos=max_anos, tolerancia_convergencia=tolerancia_convergencia,
            verificar_energia=verificar_energia,
            acelerar_convergencia=acelerar_convergencia,
        ),
        etiqueta=nombre_mapa,
    )



# ================================================================
# FASE 2b (02/10/2026): simulacion con atmosfera de dos capas
# (fase2b_atmosfera.py). La clave incluye:
#   - todas las constantes en mayusculas de ese modulo (interruptores,
#     capacidades, emisividades, viento, rugosidades...), menos DEPURAR,
#     que solo cambia lo que se imprime;
#   - v2.4.1: la HUELLA DEL CODIGO FUENTE de los modulos que hacen la
#     fisica. Antes, una correccion de la fisica que no cambiara ninguna
#     constante seguia sirviendo resultados antiguos de la cache;
#   - los argumentos de la simulacion, tambien los que tienen valor por
#     defecto (max_anos, acelerar, estado_inicial).
# Asi, cambiar CUALQUIER cosa que afecte al resultado invalida la cache
# sola, sin tener que acordarse de nada.
# ================================================================

MODULOS_FISICA = (
    "fase2b_atmosfera.py", "fase2_inercia_multicapa.py", "fase2_difusion.py",
    "fase1_geografia.py", "rejilla.py", "parametros.py", "temperatura.py",
    "orbita.py", "geometria.py", "fase30_multicapa.py",
)


def huella_codigo(modulos=MODULOS_FISICA):
    """SHA-256 del codigo fuente de los modulos de la fisica."""
    carpeta = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    for nombre in modulos:
        with open(os.path.join(carpeta, nombre), "rb") as f:
            h.update(nombre.encode())
            h.update(f.read())
    return h.hexdigest()[:24]


def clave_fase2b(datos_orbita, tipo_superficie, altitud_metros, emisividad,
                 albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
                 paso_tiempo, tolerancia_convergencia, max_anos, acelerar, estado_inicial):
    import fase2b_atmosfera as F
    from fase2_inercia_multicapa import K_DIFUSIVIDAD_TIERRA, N_CAPAS_DEFECTO
    constantes = tuple(sorted(
        (nombre, repr(valor)) for nombre, valor in vars(F).items()
        if nombre.isupper() and nombre != "DEPURAR" and isinstance(valor, (int, float, dict, tuple))
    ))
    return (
        # v3 (v2.4.1): la clave incluye la huella del codigo fuente
        "simular_fase2b_v3", datos_orbita, tipo_superficie, altitud_metros, emisividad,
        tuple(sorted(albedo_por_tipo.items())), tuple(sorted(inercia_por_tipo.items())),
        profundidad_optica, D, K_DIFUSIVIDAD_TIERRA, N_CAPAS_DEFECTO, paso_tiempo,
        tolerancia_convergencia, max_anos, acelerar, estado_inicial, constantes, huella_codigo(),
    )


def simular_fase2b_cacheada(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
    nombre_mapa="", paso_tiempo=None, tolerancia_convergencia=0.015,
    max_anos=50, acelerar=True, estado_inicial="libre",
):
    """Devuelve el dict completo de simular_fase2b()."""
    import fase2b_atmosfera as F
    from temperatura import PASO_TIEMPO as PASO_TIEMPO_DEFECTO
    if paso_tiempo is None:
        paso_tiempo = PASO_TIEMPO_DEFECTO
    clave = clave_fase2b(datos_orbita, tipo_superficie, altitud_metros, emisividad,
                         albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
                         paso_tiempo, tolerancia_convergencia, max_anos, acelerar, estado_inicial)
    return cargar_o_calcular(
        "simulacion_fase2b", clave,
        lambda: F.simular_fase2b(
            datos_orbita, tipo_superficie, altitud_metros, emisividad,
            albedo_por_tipo, inercia_por_tipo, profundidad_optica, D,
            paso_tiempo=paso_tiempo, max_anos=max_anos, tolerancia_convergencia=tolerancia_convergencia,
            acelerar=acelerar, estado_inicial=estado_inicial,
        ),
        etiqueta=nombre_mapa,
    )


def simular_modelo_cacheado(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid, nombre_mapa="",
):
    """
    Punto de entrada de las herramientas (mapa_calor, consulta_punto,
    analisis_latitudes, analisis_global): mismos argumentos y misma forma
    de devolver el resultado que simular_rejilla_combinada_cacheada()
    -- (T_final, anos, registro_minima, registro_media, registro_maxima)
    -- pero con el modelo de la Fase 2b. Los registros son la
    temperatura del AIRE A 2 m (la de un parte meteorologico).
    """
    D = float(np.max(D_grid))
    r = simular_fase2b_cacheada(
        datos_orbita, tipo_superficie, altitud_metros, emisividad,
        albedo_por_tipo, inercia_por_tipo, profundidad_optica, D, nombre_mapa=nombre_mapa,
    )
    return r["T_final"], r["anos"], r["reg_min"], r["reg_media"], r["reg_max"]

if __name__ == "__main__":
    # Auto-prueba rapida y barata: NO ejecuta la simulacion de rejilla
    # completa (seria repetir minutos de calculo solo para probar la
    # cache) -- prueba el mecanismo generico cargar_o_calcular() con un
    # calculo trivial, y por separado la cache real de la orbita (que si
    # es barata, unos segundos).
    print("--- Prueba 1: cargar_o_calcular() generico, con una funcion trivial ---")
    contador_llamadas = {"n": 0}

    def calculo_trivial():
        contador_llamadas["n"] += 1
        return np.array([1.0, 2.0, 3.0]) * contador_llamadas["n"]

    r1 = cargar_o_calcular("prueba", (1, 2.5, "x", np.array([1, 2, 3])), calculo_trivial)
    r2 = cargar_o_calcular("prueba", (1, 2.5, "x", np.array([1, 2, 3])), calculo_trivial)
    r3 = cargar_o_calcular("prueba", (1, 2.5, "y", np.array([1, 2, 3])), calculo_trivial)  # clave distinta
    print(f"Llamadas reales a calculo_trivial(): {contador_llamadas['n']} (deberia ser 2: r1/r2 comparten cache, r3 no)")
    print("OK: cache funciona." if contador_llamadas["n"] == 2 and np.array_equal(r1, r2) and not np.array_equal(r1, r3)
          else "AVISO: revisar cargar_o_calcular().")

    print("\n--- Prueba 2: precalcular_orbita_cacheada(), primera vez vs segunda vez ---")
    import time
    from parametros import S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR

    inicio = time.time()
    datos_a = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    duracion_a = time.time() - inicio

    inicio = time.time()
    datos_b = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    duracion_b = time.time() - inicio

    print(f"Primera llamada: {duracion_a:.3f} s | Segunda llamada (deberia ser cache): {duracion_b:.3f} s")
    print(f"Mismos datos: {datos_a == datos_b}")
    print("OK: segunda llamada mucho mas rapida y datos identicos." if duracion_b < duracion_a / 3 and datos_a == datos_b
          else "AVISO: revisar -- la segunda llamada deberia ser casi instantanea.")

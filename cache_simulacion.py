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

CARPETA_CACHE = os.path.join("outputs", "cache")


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
# ENVOLTORIOS CONCRETOS -- uno para precalcular_orbita(), otro para
# simular_rejilla_combinada_con_registro(). Mismos argumentos y mismo
# resultado que las funciones reales; la unica diferencia es que pasan
# primero por la cache.
# ================================================================

def precalcular_orbita_cacheada(luminosidad, inclinacion_axial_rad, semieje):
    from parametros import ORBITA_PERIODO
    from temperatura import precalcular_orbita, PASO_TIEMPO

    # luminosidad, inclinacion_axial_rad y semieje ya llegan como
    # argumentos, pero ORBITA_PERIODO y PASO_TIEMPO se leen como globales
    # dentro de precalcular_orbita() -- hay que incluirlos aqui a mano
    # o un cambio en cualquiera de los dos no se detectaria.
    clave = ("precalcular_orbita_v1", luminosidad, inclinacion_axial_rad, semieje, ORBITA_PERIODO, PASO_TIEMPO)
    return cargar_o_calcular(
        "orbita", clave,
        lambda: precalcular_orbita(luminosidad, inclinacion_axial_rad, semieje),
    )


def simular_rejilla_combinada_cacheada(
    datos_orbita, tipo_superficie, altitud_metros, emisividad,
    albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
    K_difusividad=None, n_capas=None, paso_tiempo=None, max_anos=50,
    tolerancia_convergencia=None, verificar_energia=False, nombre_mapa="",
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
        "simular_combinada_v1", datos_orbita, tipo_superficie, altitud_metros,
        emisividad, albedo_items, inercia_items, profundidad_optica, D_grid,
        K_difusividad, n_capas, paso_tiempo, max_anos, tolerancia_convergencia,
        verificar_energia,
    )
    return cargar_o_calcular(
        "simulacion_combinada", clave,
        lambda: simular_rejilla_combinada_con_registro(
            datos_orbita, tipo_superficie, altitud_metros, emisividad,
            albedo_por_tipo, inercia_por_tipo, profundidad_optica, D_grid,
            K_difusividad=K_difusividad, n_capas=n_capas, paso_tiempo=paso_tiempo,
            max_anos=max_anos, tolerancia_convergencia=tolerancia_convergencia,
            verificar_energia=verificar_energia,
        ),
        etiqueta=nombre_mapa,
    )


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

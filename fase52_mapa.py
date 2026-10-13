# fase52_mapa.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §7, decision tecnica T9): las tablas k del
# infrarrojo con un MAPA ESPECTRAL FIJO por banda, como ecCKD (Hogan y Matricardi 2022, JAMES 14, e2022MS003033 ✅t;
# codigo Apache 2.0 leido el 13/10/2026). Lo usa herramientas/generar_tablas_k.py. No se usa todavia en el modelo.
#
# Idea: en cada banda, cada numero de onda se asigna UNA vez a un grupo ("punto g") segun lo hondo que llega su luz en
# una columna de referencia de P3N (fase52_referencia: el nivel donde el espesor optico desde arriba llega a 1;
# fase52_infrarrojo.tabla_mapa). Para cada gas, presion, temperatura (y fraccion de vapor), su coeficiente en cada
# grupo es el que da la misma transmision difusa media del grupo, pesada con Planck a esa temperatura, para la cantidad
# de ese gas que tiene una capa de la columna de referencia a esa presion (el promedio "transmission" de ecCKD); para el
# vapor, que cambia mucho de un sitio a otro, la cantidad es la de la fraccion molar de esa casilla de la tabla. En el
# modelo, el espesor optico de cada grupo es la SUMA de los gases (u_gas * k_gas), y la emision de cada grupo, su
# fraccion de Planck, que con grupos fijos es exacta.
#
# Por que (T9, medido con herramientas/medir_solapamiento.py; DISENO_FASE_5.2.md §7): con el mismo numero de puntos g,
# el error frente al linea a linea es tan pequeño como el de la mejor alternativa (la mezcla ordenada capa a capa), pero
#   - sirve para cualquier composicion (el CO2 se puede cambiar por worldbuilding sin rehacer las tablas: se suman los
#     gases), al contrario que la mezcla ordenada;
#   - cuesta lo mismo que un solo gas, al contrario que el solapamiento aleatorio con reordenacion (que habria que
#     ordenar en cada capa y cada llamada);
#   - la fraccion de Planck de cada grupo es exacta, lo que en las bandas anchas reduce el error del flujo al espacio
#     de -0,75 W/m2 a -0,04 W/m2 (banda de 10-350 cm^-1, columna de referencia).

import numpy as np

import fase52_generador as G
import fase52_infrarrojo as F
import fase52_k as K
import fase52_referencia as R


def mapa_banda(datos, banda, mtckd, n_capas=30):
    """El mapa de la banda con la columna de referencia: dict con nu, dnu (la rejilla comun, T12 con las anchuras de la
    capa de arriba de todos los gases), grupo de cada nu, pesos de cada grupo (fraccion de la banda), y la presion y la
    cantidad de cada gas de cada capa de la columna (para el promedio)."""
    col = R.columna(n_capas, Ts=R.COLUMNAS["ref"])
    nu, dnu = R.rejilla_comun(datos, banda, col)
    tau = np.zeros((n_capas, len(nu)))                 # la suma de los gases, gas a gas (asi cabe en memoria)
    for gas, par in datos.items():
        tau += sum(R.espesores({gas: par}, nu, col, mtckd).values(), np.zeros((n_capas, len(nu))))
    _, pesos = K.nodos_g()
    acum = np.vstack([np.zeros(len(nu)), np.cumsum(tau, axis=0)])
    grupo = F.tabla_mapa(acum, dnu, pesos)
    pesos_mapa = np.bincount(grupo, weights=dnu, minlength=len(pesos)) / dnu.sum()
    u_ref = {gas: col["vmr"][gas] * col["n_aire"] for gas in R.GASES}
    u_ref["aire"] = col["n_aire"]
    return {"banda": banda, "nu": nu, "dnu": dnu, "grupo": grupo.astype(np.int16), "pesos": pesos_mapa,
            "pm_ref": col["pm"], "u_ref": u_ref}


def u_referencia(mapa, gas, p_hPa):
    """Cantidad del gas (moleculas/cm^2) de una capa de la columna de referencia a la presion p (interpolada en ln p y
    en ln u; fuera de la columna, la de la capa del extremo)."""
    lp = np.log(mapa["pm_ref"])
    lu = np.log(np.maximum(mapa["u_ref"][gas], 1e-300))
    return float(np.exp(np.interp(np.log(p_hPa), lp, lu)))


def k_grupos(lineas, isotopos, mapa, gas, p_hPa, T, vmr_self, u_max, umbral, mtckd=None):
    """k (cm^2/molecula) de un gas en cada grupo del mapa a (p, T, vmr_self), y la cota de lo descartado (T14).
    Vapor: lineas sin pedestal (T4) mas el continuo MT_CKD."""
    nu, dnu, grupo = mapa["nu"], mapa["dnu"], mapa["grupo"]
    n_g = len(mapa["pesos"])
    p = p_hPa / 1013.25
    k = np.zeros_like(nu)
    perdida = 0.0
    es_vapor = gas == "H2O"
    if lineas is not None:
        a, b = mapa["banda"]
        for iso, inf in isotopos.items():
            m = (lineas["iso"] == iso) & (lineas["nu0"] > a - 25.0) & (lineas["nu0"] < b + 25.0)
            if not m.any():
                continue
            sub = {kk: v[m] for kk, v in lineas.items()}
            Q_T, Q_ref = G.q_de(inf["T"], inf["Q"], T), G.q_de(inf["T"], inf["Q"], G.T_REF)
            S, c, gL, aD = G.parametros(sub, T, p, vmr_self * p, Q_T, Q_ref, inf["masa_u"])
            q, per = G.descartar(S, u_max, umbral)
            perdida += per
            if len(q):
                k += G.k_rapido(nu, S[q], c[q], gL[q], aD[q], quitar_pedestal=es_vapor)
    if es_vapor:
        k = np.maximum(k, 0.0)
        if mtckd is not None:
            from fase52_continuo import k_continuo
            k += k_continuo(nu, T, p_hPa, vmr_self, mtckd)
    if es_vapor:        # el vapor varia mucho: la cantidad del promedio es la de su propia fraccion molar (la de la tabla)
        u = max(vmr_self, 1e-12) * u_referencia(mapa, "aire", p_hPa)
    else:
        u = u_referencia(mapa, gas, p_hPa)
    peso = nu ** 3 / np.expm1(F.C2 * nu / T) * dnu
    kg = F.media_transmision_por_grupo((k * u)[None, :], peso, grupo, n_g)[0] / u
    return kg, perdida


def fracciones_planck(mapa, T):
    """(grupos, len(T)): fraccion de la emision de la banda que cae en cada grupo, a cada temperatura (exacta: los
    grupos son conjuntos fijos de numeros de onda; la suma de los grupos es 1)."""
    nu, dnu, grupo = mapa["nu"], mapa["dnu"], mapa["grupo"]
    T = np.atleast_1d(np.asarray(T, dtype=float))
    out = np.zeros((len(mapa["pesos"]), len(T)))
    for i, t in enumerate(T):
        B = nu ** 3 / np.expm1(F.C2 * nu / t) * dnu
        out[:, i] = np.bincount(grupo, weights=B, minlength=len(mapa["pesos"])) / B.sum()
    return out

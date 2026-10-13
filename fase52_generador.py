# fase52_generador.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §7): el GENERADOR de las tablas k a partir
# de las lineas de HITRAN. Calcula k(nu) de un gas en una banda a una presion y temperatura dadas, rapido y con el error
# medido frente al calculo de referencia (fase52_lineas.coeficiente_absorcion), y lo convierte en k(g) (fase52_k).
# Lo usa herramientas/generar_tablas_k.py en el PC de Carlos. No se usa todavia en el modelo.
#
# Como se acelera (decisiones tecnicas provisionales, cada una con su prueba de error en test_fase52_generador.py):
#   T12 REJILLA ADAPTADA (no uniforme): una base uniforme de paso 0,01 cm^-1 mas, alrededor de cada linea estrecha
#       (semianchura de Voigt w < 0,04 cm^-1; Olivero y Longbothum 1977 ⚠️), puntos con paso w/4 hasta 2 w y despues
#       cada vez mas separados (cada uno un 10 % mas lejos) hasta que el paso alcanza el de la base. Cada punto pesa lo que mide su celda (la mitad de la
#       distancia a cada vecino), y las tablas k usan esos pesos (fase52_k.distribucion_k(k, dnu)). Asi se resuelve el
#       centro de cada linea aunque sea estrechisima (en la alta atmosfera y a numero de onda bajo, el Doppler da
#       anchuras de 1e-5 cm^-1, que con una rejilla uniforme serian cientos de millones de puntos), sin pagar esa
#       resolucion entre lineas, donde k varia despacio.
#   T13 NUCLEO Y ALAS: cada linea, en su NUCLEO, |nu - nu0| < max(50 semianchuras, 0,25 cm^-1), sobre la rejilla
#       fina: con el perfil de Voigt a menos de 50 anchuras Doppler del centro y con el de Lorentz en el resto (alli
#       coinciden: el error relativo va como (alfa_D/dnu)^2 < 4e-4); y en las ALAS, hasta el corte (25 cm^-1, el convenio de los continuos del
#       vapor, T4), con el perfil de Lorentz, que alli coincide con el de Voigt (el error relativo va como
#       (alfa_D / dnu)^2 < 4e-4 a 50 semianchuras), sobre una rejilla GRUESA (paso 0,02 cm^-1) que despues se
#       interpola a la fina (error de la interpolacion lineal ~ (h/dnu)^2/8 < 1e-3 a 0,25 cm^-1, del ala de esa linea). Asi el coste por linea
#       no depende del paso fino en las alas, que son la mayor parte de la ventana.
#   T14 DESCARTE DE LINEAS DEBILES: se descartan las lineas cuya absorcion integrada con la columna MAXIMA del gas
#       (S(T) * u_max, en cm^-1: su anchura equivalente en el limite debil, la mayor posible) es menor que un umbral
#       por linea; y se devuelve la SUMA de lo descartado, que es una cota superior de la absorcion perdida en la banda
#       (en el limite debil la anchura equivalente es aditiva y nunca mayor que S u): queda registrada en la tabla.

import math

import numpy as np

import fase52_lineas as L

PASO_GRUESO = 0.02          # cm^-1 (T13)
NUCLEO_MIN = 0.25           # cm^-1 (T13)
VOIGT_DOPPLERS = 50.0       # (T13): Voigt solo a menos de 50 anchuras Doppler del centro; fuera, Lorentz
NUCLEO_SEMIANCHURAS = 50.0  # (T13)
FRACCION_PASO = 0.25        # paso fino / menor semianchura de Voigt (T12)


def semianchura_voigt(alfa_D, gamma_L):
    """Semianchura de Voigt (HWHM), aproximacion de Olivero y Longbothum (1977), error < 0,02 % (lo comprueba
    test_fase52_lineas)."""
    return 0.5346 * gamma_L + np.sqrt(0.2166 * gamma_L ** 2 + alfa_D ** 2)


def parametros(lineas, T, p_atm, p_self_atm, Q_T, Q_ref, masa_u):
    """S(T), centros desplazados, gamma_L y alfa_D de cada linea (las mismas formulas que fase52_lineas)."""
    S = L.intensidad(lineas["S"], lineas["nu0"], lineas["E"], T, Q_T, Q_ref)
    gL = L.hwhm_lorentz(lineas["g_air"], lineas["g_self"], lineas["n_air"], p_atm, p_self_atm, T)
    c = lineas["nu0"] + lineas["d_air"] * p_atm
    aD = L.hwhm_doppler(lineas["nu0"], masa_u, T)
    return S, c, gL, aD


def descartar(S, u_max, umbral):
    """T14: indices de las lineas que se calculan y cota superior (cm^-1) de la absorcion de las descartadas."""
    W = S * u_max
    quedan = W >= umbral
    return np.nonzero(quedan)[0], float(W[~quedan].sum())


PASO_BASE = 0.01            # cm^-1 (T12)
CRECIMIENTO = 1.1           # (T12)


def rejilla_fina(nu_a, nu_b, gL, aD):
    """Rejilla uniforme de paso 1/4 de la menor semianchura (la de referencia para las pruebas; cara en la alta
    atmosfera)."""
    paso = FRACCION_PASO * float(np.min(semianchura_voigt(aD, gL))) if len(gL) else (nu_b - nu_a) / 1000
    n = int(math.ceil((nu_b - nu_a) / paso))
    return nu_a + (np.arange(n) + 0.5) * (nu_b - nu_a) / n          # centros de intervalos iguales


def rejilla_adaptada(nu_a, nu_b, c, gL, aD):
    """T12: (nu, peso) en [nu_a, nu_b]: base uniforme mas puntos alrededor de cada linea; peso = anchura de su celda."""
    n_base = max(int(math.ceil((nu_b - nu_a) / PASO_BASE)), 2)
    partes = [nu_a + (np.arange(n_base) + 0.5) * (nu_b - nu_a) / n_base]
    w = semianchura_voigt(aD, gL)
    for ci, wi in zip(c, w):
        if wi >= PASO_BASE / FRACCION_PASO:        # linea ancha: la base ya la resuelve (paso <= w/4)
            continue
        fino = np.arange(0.0, 2.0 * wi, FRACCION_PASO * wi)
        x = [fino]
        # despues, cada punto un 10 % mas lejos, hasta que el paso alcanza el de la base
        ultimo = fino[-1]
        n_geo = max(int(math.ceil(math.log(max(PASO_BASE / (ultimo * (CRECIMIENTO - 1)), 1.0)) / math.log(CRECIMIENTO))), 0)
        if n_geo:
            x.append(ultimo * CRECIMIENTO ** np.arange(1, n_geo + 1))
        x = np.concatenate(x)
        partes.append(ci + x)
        partes.append(ci - x[1:])
    nu = np.concatenate(partes)
    nu = np.unique(nu[(nu > nu_a) & (nu < nu_b)])
    bordes = np.concatenate([[nu_a], 0.5 * (nu[1:] + nu[:-1]), [nu_b]])
    return nu, np.diff(bordes)


def k_rapido(nu, S, c, gL, aD, corte=25.0, quitar_pedestal=False):
    """k(nu) en la rejilla fina nu (centros de intervalos iguales) con nucleo de Voigt y alas de Lorentz (T13).
    Las lineas fuera de la banda tambien contribuyen con sus alas si su centro esta a menos de 'corte' de ella.
    quitar_pedestal (T4, para el vapor con MT_CKD): a cada linea se le resta, en toda su ventana, su valor en el corte
    (alli Voigt = Lorentz: el Doppler es despreciable a 25 cm^-1). Es la definicion de MT_CKD: "the Lorentz line shape
    minus a pedestal ... defined as the Lorentz line shape value at +-25 cm-1 from the line center", y el continuo se
    suma a esas contribuciones locales (hitran.org/mtckd, comprobado el 12/10/2026 ✅)."""
    nu = np.asarray(nu, dtype=float)
    k = np.zeros_like(nu)
    a, b = nu[0] - 0.5 * (nu[1] - nu[0]), nu[-1] + 0.5 * (nu[-1] - nu[-2])
    g0 = a - corte
    n_g = int(math.ceil((b - a + 2 * corte) / PASO_GRUESO)) + 1
    gruesa = g0 + np.arange(n_g) * PASO_GRUESO
    k_alas = np.zeros(n_g)
    wV = semianchura_voigt(aD, gL)
    for ci, si, gl, ad, w in zip(c, S, gL, aD, wV):
        if ci < a - corte or ci > b + corte:
            continue
        nucleo = max(NUCLEO_SEMIANCHURAS * w, NUCLEO_MIN)
        ped = si * gl / math.pi / (corte * corte + gl * gl) if quitar_pedestal else 0.0
        # nucleo, en la rejilla fina: Voigt a menos de 50 anchuras Doppler del centro (fuera, Voigt y Lorentz coinciden:
        # el error relativo va como (alfa_D/dnu)^2 < 4e-4) y Lorentz en el resto del nucleo
        i0, i1 = np.searchsorted(nu, (ci - nucleo, ci + nucleo))
        if i1 > i0:
            cv = min(VOIGT_DOPPLERS * ad, nucleo)
            v0, v1 = np.searchsorted(nu, (ci - cv, ci + cv))
            x = nu[i0:i1] - ci
            k[i0:i1] += si * gl / math.pi / (x * x + gl * gl) - ped
            if v1 > v0:
                xv = nu[v0:v1] - ci
                k[v0:v1] += si * (L.voigt(nu[v0:v1], ci, ad, gl) - gl / math.pi / (xv * xv + gl * gl))
        # alas: Lorentz en la rejilla gruesa, de nucleo a corte a cada lado
        j0, j1 = np.searchsorted(gruesa, (ci - corte, ci + corte))
        k_alas[j0:j1] += _lorentz_ala(gruesa[j0:j1] - ci, si, gl, nucleo, corte, ped)
    # las alas se pasan a la rejilla fina interpolando; junto al borde de cada nucleo la gruesa mezcla "dentro" (0) y
    # "fuera" (Lorentz) de esa linea: ahi el ala se pone exacta (_corregir_bordes)
    out = k + np.interp(nu, gruesa, k_alas)
    _corregir_bordes(out, nu, gruesa, c, S, gL, wV, a, b, corte, quitar_pedestal)
    return out


def _lorentz_ala(x, si, gl, nucleo, corte, ped=0.0):
    return np.where((np.abs(x) >= nucleo) & (np.abs(x) < corte), si * gl / math.pi / (x * x + gl * gl) - ped, 0.0)


def _corregir_bordes(out, nu, gruesa, c, S, gL, wV, a, b, corte, quitar_pedestal=False):
    """En los dos intervalos gruesos que rodean cada borde del nucleo de cada linea, quita lo interpolado de ESA linea
    y pone su ala exacta en la rejilla fina (y tambien en el borde del corte)."""
    for ci, si, gl, w in zip(c, S, gL, wV):
        if ci < a - corte or ci > b + corte:
            continue
        nucleo = max(NUCLEO_SEMIANCHURAS * w, NUCLEO_MIN)
        ped = si * gl / math.pi / (corte * corte + gl * gl) if quitar_pedestal else 0.0
        for borde in (ci - nucleo, ci + nucleo, ci - corte, ci + corte):
            jb = int(math.floor((borde - gruesa[0]) / PASO_GRUESO))
            js = np.arange(max(jb - 1, 0), min(jb + 3, len(gruesa)))
            if len(js) < 2:
                continue
            i0, i1 = np.searchsorted(nu, (gruesa[js[0]], gruesa[js[-1]]))
            if i1 <= i0:
                continue
            interp = np.interp(nu[i0:i1], gruesa[js], _lorentz_ala(gruesa[js] - ci, si, gl, nucleo, corte, ped))
            out[i0:i1] += _lorentz_ala(nu[i0:i1] - ci, si, gl, nucleo, corte, ped) - interp


# ---------------------------------------------------------------------------------------------- una casilla de la tabla
T_REF = 296.0


def q_de(tabla_T, tabla_Q, T):
    """Funcion de particion a T interpolando la tabla TIPS (cada 1 K: el error de la interpolacion lineal es < 1e-5)."""
    return float(np.interp(T, tabla_T, tabla_Q))


def k_gas_banda(lineas, isotopos, banda, p_atm, T, vmr_self, u_max, umbral, quitar_pedestal=False, corte=25.0,
                k_extra=None):
    """k(g) de un gas en una banda (a, b) a la presion p_atm (atm), temperatura T (K) y fraccion molar propia
    vmr_self (ensancha por choques consigo mismo). lineas: las del gas (fase52_hitran.leer_par). isotopos: dict
    iso -> {"masa_u", "T": tabla de T, "Q": tabla de Q}. Devuelve (k(g) en cm^2/molecula en los puntos de
    fase52_k.nodos_g(), cota de la absorcion descartada en cm^-1 (T14), lineas usadas, puntos de la rejilla).
    k_extra: funcion nu -> k (cm^2/molecula del mismo gas) que se suma antes de ordenar (el continuo del vapor)."""
    import fase52_k as K
    a, b = banda
    sel = (lineas["nu0"] > a - corte) & (lineas["nu0"] < b + corte)
    S_, c_, gL_, aD_ = [], [], [], []
    for iso, inf in isotopos.items():
        m = sel & (lineas["iso"] == iso)
        if not m.any():
            continue
        sub = {k: v[m] for k, v in lineas.items()}
        Q_T, Q_ref = q_de(inf["T"], inf["Q"], T), q_de(inf["T"], inf["Q"], T_REF)
        S, c, gL, aD = parametros(sub, T, p_atm, vmr_self * p_atm, Q_T, Q_ref, inf["masa_u"])
        S_.append(S); c_.append(c); gL_.append(gL); aD_.append(aD)
    g, _ = K.nodos_g()
    perdida, quedan = 0.0, np.zeros(0, dtype=int)
    if S_:
        S, c, gL, aD = (np.concatenate(x) for x in (S_, c_, gL_, aD_))
        quedan, perdida = descartar(S, u_max, umbral)
    if len(quedan) == 0:
        if k_extra is None:
            return np.zeros(len(g)), perdida, 0, None
        nu = np.linspace(a, b, int(math.ceil((b - a) / PASO_BASE)) + 1)          # solo el continuo: varia despacio
        nu = 0.5 * (nu[1:] + nu[:-1])
        _, _, kg = K.tabla_k(np.maximum(k_extra(nu), 0.0))
        return kg, perdida, 0, len(nu)
    S, c, gL, aD = S[quedan], c[quedan], gL[quedan], aD[quedan]
    dentro = (c > a - NUCLEO_MIN * 50) & (c < b + NUCLEO_MIN * 50)   # las que refinan la rejilla: cerca de la banda
    nu, peso = rejilla_adaptada(a, b, c[dentro], gL[dentro], aD[dentro])
    k = k_rapido(nu, S, c, gL, aD, corte=corte, quitar_pedestal=quitar_pedestal)
    if quitar_pedestal:
        k = np.maximum(k, 0.0)       # con el pedestal restado puede quedar un negativo por redondeo lejos de las lineas
    if k_extra is not None:
        k = k + np.maximum(k_extra(nu), 0.0)
    _, _, kg = K.tabla_k(k, dnu=peso)
    return kg, perdida, len(quedan), len(nu)

# fase52_lineas.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §3.1 y §7): calculo LINEA A LINEA propio de
# la absorcion de un gas a partir de los parametros de sus lineas (formato y definiciones de HITRAN). Es la base de
# las tablas k correlacionadas de la opcion (B). No se usa todavia en el modelo.
#
# Unidades (las de HITRAN): numero de onda nu en cm^-1; intensidad de linea S en cm^-1/(molecula cm^-2); anchuras en
# cm^-1/atm (semianchura a media altura, HWHM); coeficiente de absorcion k en cm^2/molecula.
#
# Fisica (a), con las definiciones de HITRAN (Simeckova et al. 2006, J. Quant. Spectrosc. Radiat. Transf. 98, 130 ⚠️;
# hitran.org/docs/definitions-and-units ⚠️):
#   S(T) = S(Tref) * Q(Tref)/Q(T) * exp(-c2 E''/T)/exp(-c2 E''/Tref) * (1 - exp(-c2 nu0/T))/(1 - exp(-c2 nu0/Tref))
#   gamma_L(p, T) = (Tref/T)^n_air * (gamma_air (p - p_self) + gamma_self p_self)            [p en atm]
#   nu0* = nu0 + delta_air p                                                                   (desplazamiento)
#   alfa_D (HWHM Doppler) = nu0/c * sqrt(2 k_B T ln2 / m)
#   perfil de Voigt normalizado (area 1): V = Re[w(z)] / (sigma sqrt(2 pi)), z = (nu - nu0* + i gamma_L)/(sigma sqrt 2),
#   sigma = alfa_D / sqrt(2 ln2), w = funcion de Faddeeva (scipy.special.wofz).
# c2 = h c / k_B = 1,438776877 cm K (CODATA 2018, exacto a partir de h, c y k_B definidos en el SI de 2019).

import math

import numpy as np
from scipy.special import wofz

C2 = 1.438776877                 # cm K (CODATA 2018)
K_B = 1.380649e-23               # J/K (exacta, SI 2019)
C_LUZ = 2.99792458e8             # m/s (exacta)
U_MASA = 1.66053906660e-27       # kg (CODATA 2018)
T_REF = 296.0                    # K (referencia de HITRAN)
P_REF_ATM = 1.0                  # atm (referencia de HITRAN)
ATM_PA = 101325.0                # Pa (exacta)


def intensidad(S_ref, nu0, E_baja, T, Q_T, Q_ref, T_ref=T_REF):
    """Intensidad de la linea a la temperatura T (definicion de HITRAN). Q_T, Q_ref: funcion de particion total
    interna del isotopologo a T y a T_ref (de las tablas TIPS)."""
    return (S_ref * (Q_ref / Q_T) * np.exp(-C2 * E_baja / T) / np.exp(-C2 * E_baja / T_ref)
            * (-np.expm1(-C2 * nu0 / T)) / (-np.expm1(-C2 * nu0 / T_ref)))


def hwhm_lorentz(gamma_air, gamma_self, n_air, p_atm, p_self_atm, T, T_ref=T_REF):
    """Semianchura de presion (HWHM, cm^-1). Convenio de HITRAN: el exponente n_air se aplica a las dos anchuras."""
    return (T_ref / T) ** n_air * (gamma_air * (p_atm - p_self_atm) + gamma_self * p_self_atm)


def hwhm_doppler(nu0, masa_u, T):
    """Semianchura Doppler (HWHM, cm^-1); masa_u en unidades de masa atomica."""
    return nu0 / C_LUZ * np.sqrt(2.0 * K_B * T * math.log(2.0) / (masa_u * U_MASA))


def voigt(nu, nu0, alfa_D, gamma_L):
    """Perfil de Voigt normalizado (area 1, en 1/cm^-1). nu puede ser un vector; nu0, alfa_D, gamma_L escalares."""
    sigma = alfa_D / math.sqrt(2.0 * math.log(2.0))
    z = ((nu - nu0) + 1j * gamma_L) / (sigma * math.sqrt(2.0))
    return wofz(z).real / (sigma * math.sqrt(2.0 * math.pi))


def coeficiente_absorcion(nu, lineas, T, p_atm, p_self_atm, Q_T, Q_ref, masa_u, corte=25.0, quitar_pedestal=False):
    """k(nu) en cm^2/molecula de un isotopologo: suma de lineas de Voigt cortadas a +-corte cm^-1 del centro.
    lineas: dict de arrays HITRAN: nu0, S, E (E''), g_air, g_self, n_air, d_air.
    quitar_pedestal: restar a cada linea su valor en el corte (convenio de los continuos tipo CKD/MT_CKD para el
    vapor, ⚠️ por comprobar en la documentacion de MT_CKD antes de usarlo)."""
    nu = np.asarray(nu, dtype=float)
    k = np.zeros_like(nu)
    S = intensidad(lineas["S"], lineas["nu0"], lineas["E"], T, Q_T, Q_ref)
    gL = hwhm_lorentz(lineas["g_air"], lineas["g_self"], lineas["n_air"], p_atm, p_self_atm, T)
    centros = lineas["nu0"] + lineas["d_air"] * p_atm
    aD = hwhm_doppler(lineas["nu0"], masa_u, T)
    for c, s, gl, ad in zip(centros, S, gL, aD):
        i0, i1 = np.searchsorted(nu, (c - corte, c + corte))
        if i1 <= i0:
            continue
        perfil = voigt(nu[i0:i1], c, ad, gl)
        if quitar_pedestal:
            perfil = perfil - voigt(np.array([c + corte]), c, ad, gl)[0]
        k[i0:i1] += s * perfil
    return k

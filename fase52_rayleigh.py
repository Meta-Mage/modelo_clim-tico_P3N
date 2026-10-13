# fase52_rayleigh.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §3.2 y §7, decision tecnica T22): la
# DISPERSION DE RAYLEIGH del aire (la que hace azul el cielo), calculada GAS A GAS para que sirva para el aire de P3N
# (menos oxigeno que la Tierra). No se usa todavia en el modelo.
#
# Seccion eficaz de un gas (por molecula; Bodhaine et al. 1999, J. Atmos. Oceanic Technol. 16, 1854, ec. 22, y Bates
# 1984, Planet. Space Sci. 32, 785 ⚠️):
#     sigma = 24 pi^3 / (lambda^4 N^2) * ((n^2 - 1)/(n^2 + 2))^2 * F(lambda)
# con n el indice de refraccion del gas medido a la densidad N (la combinacion ((n^2-1)/(n^2+2))/N no depende de la
# densidad: Lorentz-Lorenz) y F el factor de King (la correccion por la anisotropia de la molecula). En una mezcla, las
# secciones se suman pesadas por la fraccion molar de cada gas (la dispersion de moleculas distintas no interfiere;
# Eberhard 2010, Appl. Opt. 49, 1116 ⚠️).
#
# T22 (provisional). Indices de refraccion (formulas y datos de la base de refractiveindex.info, dominio publico,
# commit c5c2f18; lambda en µm; forma "formula 6": n - 1 = C0 + Σ Ci / (Di - lambda^-2)):
#   - N2: Peck y Khanna 1966 (J. Opt. Soc. Am. 56, 1059), 0 °C; medida de 0,468 a 2,06 µm, y en el ultravioleta
#     coincide con Griesmann y Burnett 1999 (0,145-0,270 µm) dentro del 0,03 % (comprobado el 13/10/2026);
#   - O2: Zhang et al. 2008 (Appl. Opt. 47, 3143), 20 °C, de 0,4 a 1,8 µm; por debajo de 0,288 µm, la tabla de Smith
#     et al. 1976 (Phys. Rev. A 13, 1422; 0 °C), porque la formula de Zhang extrapolada se queda corta en el
#     ultravioleta (−3,5 % a 0,2 µm); en 0,288 µm las dos difieren un 0,12 %;
#   - Ar: Peck y Fisher 1964 (J. Opt. Soc. Am. 54, 1362), 15 °C;
#   - CO2: Bideau-Mehu et al. 1973 (Opt. Commun. 9, 432), 0 °C.
#   Comprobado: la mezcla de la Tierra (78,084 % N2, 20,946 % O2, 0,934 % Ar y 300 ppm de CO2) reproduce el indice
#   del aire de Peck y Reeder 1972 (la formula de Bodhaine et al., ec. 18) dentro del 0,1 % de 0,25 a 1,5 µm.
# Factores de King (Bates 1984, como los usan Bodhaine et al., ec. 5, 6 y 23): N2 1,034 + 3,17e-4 lambda^-2;
# O2 1,096 + 1,385e-3 lambda^-2 + 1,448e-4 lambda^-4; Ar 1 (un atomo); CO2 1,15.
# Densidad de cada medida: gas ideal, N = p / (k T) (la correccion por gas real es < 0,05 %).
# El vapor de agua no se incluye (cambiaria la seccion del aire < 1 % con un 4 % de vapor: limitacion declarada ⚠️).

import numpy as np

K_B = 1.380649e-23            # J/K (exacta, SI 2019)
N_A = 6.02214076e23           # 1/mol (exacta)

# (coeficientes de la formula 6, temperatura K, presion Pa)
FORMULAS = {
    "N2": ((6.8552e-5, 3.243157e-2, 144.0), 273.15, 101325.0),
    "O2": ((1.181494e-4, 9.708931e-3, 75.4), 293.15, 101325.0),
    "Ar": ((6.432135e-5, 2.8606021e-2, 144.0), 288.15, 101325.0),
    "CO2": ((0.0, 6.99100e-2, 166.175, 1.44720e-3, 79.609, 6.42941e-5, 56.3064, 5.21306e-5, 46.0196, 1.46847e-6,
             0.0584738), 273.15, 101325.0),
}
# O2 en el ultravioleta: Smith et al. 1976 (lambda µm, n), 0 °C y 760 torr
SMITH_O2 = np.array([(0.1843, 1.0003792), (0.19185, 1.0003612), (0.19776, 1.0003482), (0.19792, 1.0003479),
                     (0.19806, 1.0003476), (0.19832, 1.0003471), (0.19864, 1.0003467), (0.1989, 1.0003462),
                     (0.20588, 1.0003361), (0.21248, 1.0003285), (0.22087, 1.0003208), (0.22116, 1.0003206),
                     (0.22174, 1.0003201), (0.24359, 1.0003065), (0.25077, 1.0003033), (0.25151, 1.000303),
                     (0.25169, 1.0003029), (0.252, 1.0003027), (0.25249, 1.0003027), (0.25293, 1.0003024),
                     (0.26321, 1.0002986), (0.28824, 1.0002916)])
LAMBDA_SMITH_MAX = 0.28824
# composicion de referencia de la Tierra (Bodhaine et al. 1999, ec. 23; CO2 aparte)
TIERRA = {"N2": 0.78084, "O2": 0.20946, "Ar": 0.00934}


def _formula6(c, lam):
    lam = np.asarray(lam, dtype=float)
    out = np.full(lam.shape, c[0], dtype=float)
    for i in range(1, len(c), 2):
        out = out + c[i] / (c[i + 1] - lam ** -2.0)
    return out


def refractividad(gas, lam_um):
    """(n - 1, N en m^-3) del gas a la longitud de onda lam_um (µm), en las condiciones de su medida."""
    lam = np.asarray(lam_um, dtype=float)
    c, T, p = FORMULAS[gas]
    nm1 = _formula6(c, lam)
    N = np.full(lam.shape, p / (K_B * T))
    if gas == "O2":
        uv = lam < LAMBDA_SMITH_MAX
        if np.any(uv):
            nm1 = np.where(uv, np.interp(lam, SMITH_O2[:, 0], SMITH_O2[:, 1] - 1.0), nm1)
            N = np.where(uv, 101325.0 / (K_B * 273.15), N)
    return nm1, N


def factor_king(gas, lam_um):
    lam = np.asarray(lam_um, dtype=float)
    if gas == "N2":
        return 1.034 + 3.17e-4 * lam ** -2
    if gas == "O2":
        return 1.096 + 1.385e-3 * lam ** -2 + 1.448e-4 * lam ** -4
    if gas == "Ar":
        return np.ones_like(lam)
    if gas == "CO2":
        return np.full(lam.shape, 1.15)
    raise ValueError(f"gas sin datos de Rayleigh: {gas}")


def seccion_gas(gas, lam_um):
    """Seccion eficaz de Rayleigh de una molecula del gas (cm^2)."""
    lam_cm = np.asarray(lam_um, dtype=float) * 1e-4
    nm1, N = refractividad(gas, lam_um)
    n2 = (1.0 + nm1) ** 2
    N_cm = N * 1e-6
    return 24.0 * np.pi ** 3 / (lam_cm ** 4 * N_cm ** 2) * ((n2 - 1.0) / (n2 + 2.0)) ** 2 * factor_king(gas, lam_um)


def seccion_aire(lam_um, fracciones):
    """Seccion eficaz media por molecula de la mezcla (cm^2): suma de cada gas pesada por su fraccion molar
    (fracciones: {gas: fraccion}; se normalizan a 1, asi que los gases sin Rayleigh propio, como el CH4, se reparten)."""
    tot = sum(fracciones.values())
    return sum(x / tot * seccion_gas(g, lam_um) for g, x in fracciones.items())


def composicion(o2, co2, ar=TIERRA["Ar"]):
    """Fracciones molares del aire seco con O2, CO2 y Ar dados y el resto N2 (los gases traza no cuentan)."""
    return {"N2": 1.0 - o2 - co2 - ar, "O2": o2, "Ar": ar, "CO2": co2}


def espesor_optico(lam_um, p_Pa, g, fracciones, masa_molar_kg=None):
    """Espesor optico de Rayleigh de una columna de aire de presion p (Pa) en superficie y gravedad g (m/s^2):
    sigma * (moleculas por cm^2) = sigma * p N_A / (M g) / 1e4."""
    if masa_molar_kg is None:
        masas = {"N2": 28.0134e-3, "O2": 31.9988e-3, "Ar": 39.948e-3, "CO2": 44.0095e-3}
        tot = sum(fracciones.values())
        masa_molar_kg = sum(x / tot * masas[g_] for g_, x in fracciones.items())
    columna = p_Pa * N_A / (masa_molar_kg * g) / 1e4
    return seccion_aire(lam_um, fracciones) * columna

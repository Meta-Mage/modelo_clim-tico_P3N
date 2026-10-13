# fase52_optica_nubes.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §3.4 y §7, decision tecnica T18): las
# propiedades opticas de las GOTAS DE AGUA de las nubes POR BANDA, a partir de la teoria de Mie propia
# (fase52_mie.py) y de los indices de refraccion medidos (fase52_indices.py, T17). No se usa todavia en el modelo.
# El hielo no es esferico y tendra su propio tratamiento (T16).
#
# Para una poblacion de gotas (gamma de Hansen y Travis, radio efectivo r_e y varianza efectiva v_e) a la temperatura
# T, en la banda [nu1, nu2] (cm^-1), devuelve:
#   - el coeficiente de extincion MASICO beta (m^2 por kg de agua liquida): el espesor optico de una capa con LWP kg/m^2
#     de agua es tau = beta * LWP, con beta = 3 Q_ext / (4 rho_agua r_e) (Stephens 1978 ⚠️; fase52_mie.optica_poblacion);
#   - el albedo de dispersion simple omega y el factor de asimetria g.
#
# T18 (provisional): promedio "fino" en la banda (thin averaging, Edwards y Slingo 1996, Q. J. R. Meteorol. Soc. 122,
# 689 ⚠️), pesado con la fuente de la banda S(nu): la luz de la estrella en la onda corta, Planck a la temperatura de
# la nube en el infrarrojo:
#     beta = Σ S beta(nu) / Σ S;   omega = Σ S beta omega / Σ S beta;   g = Σ S beta omega g / Σ S beta omega.
# Es exacto en el limite de capa fina; para nubes gruesas el error de cualquier promedio de banda se medira en la
# validacion (§4) frente al calculo espectral, no se supone.

import numpy as np

import fase52_indices as IND
import fase52_mie as MIE

RHO_AGUA = 1000.0            # kg/m^3
C2 = 1.438776877             # cm K (CODATA 2018)


def planck_nu(nu, T):
    """Forma espectral de Planck por numero de onda (sin constantes: solo sirve de peso)."""
    nu = np.asarray(nu, dtype=float)
    return nu ** 3 / np.expm1(C2 * nu / T)


def puntos_banda(nu1, nu2, n):
    """n numeros de onda en el centro de n tramos iguales de la banda (regla del punto medio)."""
    d = (nu2 - nu1) / n
    return nu1 + d * (np.arange(n) + 0.5)


def optica_monocromatica_agua(nu, r_e_um, T, v_e=0.1, n_radios=400):
    """(beta en m^2/kg, omega, g) de la poblacion de gotas al numero de onda nu (cm^-1)."""
    lam_um = 1e4 / nu
    m = complex(IND.indice_agua(lam_um, T))
    q_ext, omega, g = MIE.optica_poblacion(m, lam_um, r_e_um, v_e=v_e, n_radios=n_radios)
    beta = 3.0 * q_ext / (4.0 * RHO_AGUA * r_e_um * 1e-6)
    return beta, omega, g


def optica_banda_agua(nu1, nu2, r_e_um, T, peso=None, n_nu=20, v_e=0.1, n_radios=400):
    """(beta, omega, g) de la banda [nu1, nu2] (cm^-1) con el promedio de T18. peso: funcion S(nu) de la fuente
    (por defecto, Planck a la temperatura T de la nube)."""
    nu = puntos_banda(nu1, nu2, n_nu)
    S = planck_nu(nu, T) if peso is None else np.asarray(peso(nu), dtype=float)
    if np.any(S < 0) or S.sum() <= 0:
        raise ValueError("el peso de la banda tiene que ser positivo")
    res = np.array([optica_monocromatica_agua(x, r_e_um, T, v_e, n_radios) for x in nu])
    b, w, g = res[:, 0], res[:, 1], res[:, 2]
    beta = np.sum(S * b) / np.sum(S)
    omega = np.sum(S * b * w) / np.sum(S * b)
    gg = np.sum(S * b * w * g) / np.sum(S * b * w)
    return beta, omega, gg


class TablaGotas:
    """La tabla de herramientas/generar_optica_gotas.py (datos/optica_gotas_infrarrojo.npz) y su interpolacion:
    lineal en ln r_e y en T, de beta * r_e (casi constante: Q_ext / rho es lo que varia), de omega y de g. Fuera de la
    rejilla, el extremo (limitacion: por debajo de 2 µm o por encima de 40 µm de radio efectivo)."""

    def __init__(self, ruta=None):
        import os
        ruta = ruta or os.path.join(os.path.dirname(os.path.abspath(__file__)), "datos", "optica_gotas_infrarrojo.npz")
        with np.load(ruta) as d:
            self.beta, self.omega, self.g = d["beta_m2_kg"], d["omega"], d["g"]
            self.r_e, self.T = d["r_e_um"], d["T_K"]
            self.bandas = [tuple(b) for b in d["bandas"]]

    @staticmethod
    def _peso(x, nodos):
        x = float(np.clip(x, nodos[0], nodos[-1]))
        i = int(np.clip(np.searchsorted(nodos, x) - 1, 0, len(nodos) - 2))
        return i, (x - nodos[i]) / (nodos[i + 1] - nodos[i])

    def propiedades(self, r_e_um, T):
        """(beta m^2/kg, omega, g), cada uno (bandas,), para un radio efectivo y una temperatura."""
        ir, fr = self._peso(np.log(r_e_um), np.log(self.r_e))
        it, ft = self._peso(T, self.T)
        r_eff = float(np.clip(r_e_um, self.r_e[0], self.r_e[-1]))

        def bil(a):
            return ((1 - fr) * (1 - ft) * a[:, ir, it] + fr * (1 - ft) * a[:, ir + 1, it]
                    + (1 - fr) * ft * a[:, ir, it + 1] + fr * ft * a[:, ir + 1, it + 1])
        br = bil(self.beta * self.r_e[None, :, None])
        return br / r_eff, bil(self.omega), bil(self.g)


# ------------------------------------------------------------------------------------------- radio efectivo (T24)
# T24 (provisional): el radio efectivo de las gotas segun Martin, Johnson y Spice 1994 (J. Atmos. Sci. 51, 1823 ✅t:
# r_e = (3 LWC / (4 pi rho_w N k))^(1/3), con k = (r_vol/r_e)^3), con las constantes que usa el IFS del ECMWF
# (ecRad, ifs/liquid_effective_radius.F90 y yoerad.F90, Apache 2.0, leidos el 13/10/2026 ✅):
#   - gotas por cm^3 a partir de los nucleos de condensacion: oceano N = -1,15e-3 CCN^2 + 0,963 CCN + 5,30, con
#     CCN = 50 cm^-3 (N ≈ 51) y k = 0,77; tierra N = -2,10e-4 CCN^2 + 0,568 CCN - 27,9, con CCN = 900 cm^-3
#     (N ≈ 313) y k = 0,69;
#   - la correccion por la lluvia de Wood 2000 (ec. 19), como el IFS;
#   - r_e entre 4 y 30 µm.
# ⚠️ Son valores de la Tierra ACTUAL (calibracion terrestre, c). P3N tiene aire natural preindustrial (ficha §6.4):
# sobre tierra seguramente hay menos gotas que 313 cm^-3 (gotas mas grandes); la cifra preindustrial esta por buscar.
CCN_OCEANO, CCN_TIERRA = 50.0, 900.0
RE_MIN_UM, RE_MAX_UM = 4.0, 30.0


def gotas_por_cm3(tierra, ccn=None):
    """Numero de gotas (cm^-3) a partir de los nucleos de condensacion (formulas del IFS)."""
    if tierra:
        c = CCN_TIERRA if ccn is None else ccn
        return -2.10e-4 * c * c + 0.568 * c - 27.9
    c = CCN_OCEANO if ccn is None else ccn
    return -1.15e-3 * c * c + 0.963 * c + 5.30


def radio_efectivo_um(lwc_gm3, tierra, rwc_gm3=0.0, ccn=None):
    """Radio efectivo (µm) de las gotas de una nube con agua liquida lwc (g/m^3, dentro de la nube) y lluvia rwc."""
    k = 0.69 if tierra else 0.77
    N = gotas_por_cm3(tierra, ccn)
    lwc = np.asarray(lwc_gm3, dtype=float)
    rwc = np.asarray(rwc_gm3, dtype=float)
    ratio = (0.222 / k) ** 0.333
    lluvia = np.where(lwc > 1e-6, rwc / np.maximum(lwc, 1e-30), 0.0)
    wood = (1.0 + lluvia) ** 0.666 / (1.0 + 0.2 * ratio * lluvia)
    r3 = 3.0 * (lwc + rwc) / (4.0 * np.pi * N * k)          # g/m^3 y cm^-3 con rho_w = 1e6 g/m^3: en (100 µm)^3
    re = wood * 100.0 * np.cbrt(np.maximum(r3, 0.0))
    return np.clip(np.where(r3 > 1e-30, re, RE_MIN_UM), RE_MIN_UM, RE_MAX_UM)

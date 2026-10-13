# fase52_mie.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §3.4): la teoria de MIE propia, para calcular
# como dispersan y absorben la luz las gotas de las nubes (y, mas adelante, otras particulas esfericas): eficiencias de
# extincion y de dispersion y factor de asimetria g de una esfera de radio r e indice de refraccion complejo m a la
# longitud de onda lambda. Es la base de las propiedades opticas de las nubes (decision tecnica T16, provisional:
# calcularlas con Mie y los indices de refraccion medidos del agua, en vez de copiar una parametrizacion; el hielo no es
# esferico y tendra su propio tratamiento). No se usa todavia en el modelo.
#
# Algoritmo: BHMIE de Bohren y Huffman (1983), "Absorption and Scattering of Light by Small Particles", apendice A
# (convenio m = n + i k con k >= 0 para un medio que absorbe), con:
#   - numero de terminos nstop = x + 4 x^(1/3) + 2 (Wiscombe 1980, NCAR/TN-140+STR ⚠️);
#   - derivada logaritmica D_n(m x) por recurrencia hacia abajo desde nmx = 1,1 max(nstop, |m x|) + 50 (estable; el
#     arranque de Bohren y Huffman, max + 15, no converge del todo con m casi real y x grande: medido, 5e-5 en Q_ext a
#     x = 100 y 3e-4 a x = 1000);
#   - psi_n y chi_n (Riccati-Bessel) por recurrencia hacia arriba.
# Se comprueba (test_fase52_mie.py) frente a miepython 3.3 (Prahl, licencia MIT ⚠️; convenio m = n - i k) con valores
# tabulados en la prueba, y frente a los limites de Rayleigh (x << 1) y de la paradoja de la extincion (Q_ext -> 2).

import numpy as np


def mie(m, x):
    """(Q_ext, Q_sca, g) de una esfera con indice complejo m (m = n + i k, k >= 0) y parametro de tamaño
    x = 2 pi r / lambda (en el medio)."""
    m = complex(m)
    x = float(x)
    if x <= 0:
        return 0.0, 0.0, 0.0
    mx = m * x
    nstop = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    nmx = int(1.1 * max(nstop, abs(mx))) + 50          # arranque holgado: con m casi real, +15 no basta (medido)
    # derivada logaritmica D_n(mx), hacia abajo
    D = np.zeros(nmx + 1, dtype=complex)
    for n in range(nmx, 0, -1):
        D[n - 1] = n / mx - 1.0 / (D[n] + n / mx)
    # psi y chi de Riccati-Bessel hacia arriba: psi_{-1} = cos x, psi_0 = sin x; chi_{-1} = -sin x, chi_0 = cos x
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = complex(psi1, -chi1)
    qext = qsca = gq = 0.0
    an_ant = bn_ant = None
    for n in range(1, nstop + 1):
        fn = (2.0 * n + 1.0) / (n * (n + 1.0))
        psi = (2.0 * n - 1.0) * psi1 / x - psi0
        chi = (2.0 * n - 1.0) * chi1 / x - chi0
        xi = complex(psi, -chi)
        an = ((D[n] / m + n / x) * psi - psi1) / ((D[n] / m + n / x) * xi - xi1)
        bn = ((m * D[n] + n / x) * psi - psi1) / ((m * D[n] + n / x) * xi - xi1)
        qext += (2.0 * n + 1.0) * (an.real + bn.real)
        qsca += (2.0 * n + 1.0) * (abs(an) ** 2 + abs(bn) ** 2)
        gq += fn * (an.real * bn.real + an.imag * bn.imag)
        if an_ant is not None:
            gq += ((n - 1.0) * (n + 1.0) / n) * (an_ant.real * an.real + an_ant.imag * an.imag
                                                 + bn_ant.real * bn.real + bn_ant.imag * bn.imag)
        an_ant, bn_ant = an, bn
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = complex(psi1, -chi1)
    qext *= 2.0 / x ** 2
    qsca *= 2.0 / x ** 2
    g = 4.0 * gq / (qsca * x ** 2) if qsca > 0 else 0.0
    return qext, qsca, g


# ---------------------------------------------------------------------------------------------- poblacion de gotas
def distribucion_gamma(r, r_e, v_e):
    """Distribucion gamma de radios de Hansen y Travis (1974, Space Sci. Rev. 16, 527 ⚠️): n(r) ∝ r^((1-3v)/v) e^(-r/(a v))
    con radio efectivo a = r_e y varianza efectiva v = v_e (sin normalizar)."""
    r = np.asarray(r, dtype=float)
    return r ** ((1.0 - 3.0 * v_e) / v_e) * np.exp(-r / (r_e * v_e))


def optica_poblacion(m, lam, r_e, v_e=0.1, n_radios=400):
    """Propiedades medias de una poblacion de gotas (gamma de Hansen y Travis) a la longitud de onda lam (misma unidad
    que r_e): Q_ext efectiva (media pesada por la seccion geometrica), albedo de dispersion simple omega y g (media
    pesada por la seccion de dispersion). Con Q_ext efectiva, el espesor optico de una nube de agua liquida LWP es
    tau = 3 LWP Q_ext / (4 rho_agua r_e) (Stephens 1978 ⚠️)."""
    r_max = r_e * (1.0 + 12.0 * v_e ** 0.5) + 6.0 * r_e * v_e
    r = np.linspace(r_e * 1e-3, r_max, n_radios)
    n = distribucion_gamma(r, r_e, v_e)
    area = np.pi * r ** 2 * n
    q = np.array([mie(m, 2.0 * np.pi * ri / lam) for ri in r])
    sec_ext, sec_sca = area * q[:, 0], area * q[:, 1]
    q_ext = np.trapezoid(sec_ext, r) / np.trapezoid(area, r)
    omega = np.trapezoid(sec_sca, r) / np.trapezoid(sec_ext, r)
    g = np.trapezoid(sec_sca * q[:, 2], r) / np.trapezoid(sec_sca, r)
    return q_ext, omega, g

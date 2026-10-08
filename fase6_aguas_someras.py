# fase6_aguas_someras.py -- Fase 6, etapa 6.1: ecuaciones de aguas someras en la esfera.
#
# Primer ladrillo del nucleo dinamico PROPIO de M3N (DISENO_FASE_6.md, DISENO_FASE_6.1.md).
# Una sola capa de fluido sobre la esfera en rotacion; sirve para construir y validar, con soluciones
# conocidas, lo que despues usara el nucleo de 20 capas: la rejilla escalonada, el termino de
# Coriolis y vorticidad, el gradiente de presion, la conservacion de la masa y de la energia, el
# filtro polar y el avance en el tiempo.
#
# Discretizacion (eleccion provisional de la IA para la 6.1, a confirmar por Carlos):
#   - Rejilla: la de M3N (por defecto 72 x 36 celdas de 5 grados, fila 0 = norte), con las variables
#     escalonadas segun la rejilla C de Arakawa (Arakawa y Lamb 1977 ✅): h en el centro de cada
#     celda, u en las caras oeste/este, v en las caras norte/sur. En los polos (caras de longitud
#     cero) no hay v ni flujo. La resolucion es un parametro solo para medir la convergencia.
#   - Forma "vectorial invariante":
#         du/dt =  q (h v) - d/dx (g (h + hs) + K)
#         dv/dt = -q (h u) - d/dy (g (h + hs) + K)
#         dh/dt = - div(h v)
#     con q = (zeta + f)/h la vorticidad potencial en las esquinas (zeta por circulacion, teorema de
#     Stokes) y el producto q(hv) en la forma de Sadourny que conserva la energia (Sadourny 1975 ⚠️
#     cita por verificar). DEMOSTRACION PROPIA (DISENO_FASE_6.1.md §3): con h en las caras = media
#     aritmetica, K = (1/2A) sum_caras (A_cara u^2)/2 y estos pesos, la energia total se conserva
#     EXACTAMENTE en el sistema semidiscreto; solo el paso de tiempo y el filtro polar la alteran.
#   - Masa: forma de flujo -> la masa total se conserva exactamente (redondeo).
#   - Filtro polar de Fourier sobre las tendencias, en las filas por encima de lat_filtro: amortigua
#     las ondas zonales cortas que, con las celdas estrechas de cerca del polo, obligarian a un paso
#     de tiempo minusculo (remedio clasico de las rejillas latitud-longitud; planetWRF ✅). No toca la
#     media zonal (onda 0), asi que no cambia la masa.
#   - Tiempo: leapfrog con el filtro RAW de Williams (2009) ✅ (alfa = 0,53, nu = 0,2); el primer paso,
#     Euler hacia delante.
#
# Unidades SI. Longitudes: lambda en [-pi, pi), centros en (j + 0,5)*dlambda - pi (como rejilla.py).

import math
import numpy as np

LAT_FILTRO_DEFECTO = 60.0          # grados: el filtro polar actua por encima de esta latitud
RAW_NU, RAW_ALFA = 0.2, 0.53       # Williams (2009) ✅


class Rejilla:
    def __init__(self, radio, filas=36):
        a = self.a = radio
        self.filas, self.columnas = filas, 2 * filas
        nf, nc = self.filas, self.columnas
        self.dphi = dphi = math.pi / nf
        self.dlam = dlam = 2 * math.pi / nc
        self.phi_c = math.pi / 2 - (np.arange(nf) + 0.5) * dphi               # centros
        self.phi_f = math.pi / 2 - np.arange(nf + 1) * dphi                    # caras, de +90 a -90
        self.lam_c = (np.arange(nc) + 0.5) * dlam - math.pi
        self.lam_u = np.arange(nc) * dlam - math.pi                            # cara oeste de cada celda
        self.cos_c = np.cos(self.phi_c)
        self.cos_v = np.cos(self.phi_f[1:-1])                                  # caras interiores
        unos = np.ones((1, nc))
        self.area = (a * a * dlam * (np.sin(self.phi_f[:-1]) - np.sin(self.phi_f[1:])))[:, None] * unos
        self.dx_u = (a * self.cos_c * dlam)[:, None]                           # distancia entre centros (zonal)
        self.dy = a * dphi                                                     # distancia entre centros (meridiana)
        self.L_u = a * dphi                                                    # longitud de las caras de u
        self.L_v = (a * self.cos_v * dlam)[:, None]                            # longitud de las caras de v
        self.L_n = (a * self.cos_c * dlam)[:, None]                            # tramos zonales del contorno dual
        self.A_u = self.dx_u * self.dy * unos                                  # area asociada a cada u
        self.A_v = self.L_v * self.dy * unos                                   # area asociada a cada v
        # celda dual alrededor de cada esquina interior (entre filas i, i+1 y columnas j-1, j)
        self.area_dual = (a * a * dlam * (np.sin(self.phi_c[:-1]) - np.sin(self.phi_c[1:])))[:, None] * unos

    def malla(self, donde):
        """Coordenadas (lambda, phi) de los puntos de h ('h'), u ('u') o v ('v')."""
        if donde == "h":
            return np.meshgrid(self.lam_c, self.phi_c)
        if donde == "u":
            return np.meshgrid(self.lam_u, self.phi_c)
        return np.meshgrid(self.lam_c, self.phi_f[1:-1])


def _filtro_fourier(campo, cos_filas, cos_ref, dlam):
    """Filtro polar: en cada fila con cos(lat) < cos_ref multiplica la onda zonal k por
    S(k) = min(1, cos(lat) / (cos_ref * |sin(k*dlambda/2)|)). La onda 0 (media zonal) no se toca."""
    filtrar = cos_filas < cos_ref
    if not filtrar.any():
        return campo
    n = campo.shape[1]
    k = np.arange(n // 2 + 1)
    s = np.abs(np.sin(k * dlam / 2))
    s[0] = 1.0
    S = np.minimum(1.0, cos_filas[filtrar][:, None] / (cos_ref * s[None, :]))
    S[:, 0] = 1.0
    out = campo.copy()
    out[filtrar] = np.fft.irfft(np.fft.rfft(campo[filtrar], axis=1) * S, n=n, axis=1)
    return out


class AguasSomeras:
    """Ecuaciones de aguas someras en la rejilla C de M3N. Estado: h (F, C) [m], u (F, C) [m/s, cara
    oeste de cada celda], v (F-1, C) [m/s, cara sur de las filas 0..F-2, positivo hacia el norte]."""

    def __init__(self, radio, omega, g, hs=None, lat_filtro=LAT_FILTRO_DEFECTO, alfa_coriolis=0.0,
                 filas=36, viento_fijo=False):
        self.R = R = Rejilla(radio, filas)
        self.omega, self.g = omega, g
        self.hs = np.zeros((R.filas, R.columnas)) if hs is None else hs
        # Coriolis en las esquinas: f = 2 Omega (-cos(lambda) cos(phi) sin(alfa) + sin(phi) cos(alfa)).
        # alfa = 0 es el planeta real; alfa != 0 solo para el caso 2 de Williamson con el eje girado,
        # donde el eje de rotacion se gira junto con el flujo (Williamson et al. 1992 ✅ via SWEET).
        phi_e = R.phi_f[1:-1][:, None]
        lam_e = R.lam_u[None, :]
        self.f_esq = 2 * omega * (-np.cos(lam_e) * np.cos(phi_e) * math.sin(alfa_coriolis)
                                  + np.sin(phi_e) * math.cos(alfa_coriolis))
        self.cos_ref = math.cos(math.radians(lat_filtro))
        self.viento_fijo = viento_fijo       # caso 1: solo se transporta h con el viento dado

    # ---------------------------------------------------------------
    def flujos(self, h, u, v):
        """Flujos de volumen por cara [m3/s]: F_u por la cara oeste, F_v por la cara sur (hacia el norte +)."""
        R = self.R
        h_u = 0.5 * (h + np.roll(h, 1, axis=1))
        h_v = 0.5 * (h[:-1] + h[1:])
        return u * h_u * R.L_u, v * h_v * R.L_v

    def vorticidad(self, u, v):
        """Vorticidad relativa en las esquinas interiores (circulacion / area de la celda dual)."""
        R = self.R
        circ = (R.dy * (v - np.roll(v, 1, axis=1)) - u[:-1] * R.L_n[:-1] + u[1:] * R.L_n[1:])
        return circ / R.area_dual

    def h_esquina(self, h):
        return 0.25 * (h[:-1] + h[1:] + np.roll(h[:-1], 1, axis=1) + np.roll(h[1:], 1, axis=1))

    def energia_cinetica(self, u, v):
        """K en el centro de cada celda, ponderada por areas (consistente con la conservacion de la energia)."""
        R = self.R
        ku = R.A_u * u * u
        kv = np.zeros((R.filas + 1, R.columnas)); kv[1:-1] = R.A_v * v * v
        return (0.5 * (ku + np.roll(ku, -1, axis=1)) + 0.5 * (kv[:-1] + kv[1:])) / (2 * R.area)

    # ---------------------------------------------------------------
    def tendencias(self, h, u, v):
        R = self.R
        F_u, F_v = self.flujos(h, u, v)
        # masa: forma de flujo (conserva exactamente)
        div = np.roll(F_u, -1, axis=1) - F_u
        div[:-1] -= F_v                    # la fila i gana lo que entra por su cara sur
        div[1:] += F_v                     # la fila i+1 pierde lo que sale por su cara norte
        dh = _filtro_fourier(-div / R.area, R.cos_c, self.cos_ref, R.dlam)
        if self.viento_fijo:
            return dh, np.zeros_like(u), np.zeros_like(v)
        q = (self.vorticidad(u, v) + self.f_esq) / self.h_esquina(h)          # vorticidad potencial
        # termino de vorticidad de Sadourny (conserva la energia; ver la cabecera)
        qV = q * 0.5 * (F_v + np.roll(F_v, 1, axis=1))
        Fu = np.zeros_like(u)
        Fu[1:] += 0.5 * qV                 # esquina norte de las filas 1..F-1
        Fu[:-1] += 0.5 * qV                # esquina sur de las filas 0..F-2
        Fu /= R.dx_u
        qU = q * 0.5 * (F_u[:-1] + F_u[1:])
        Fv = -0.5 * (qU + np.roll(qU, -1, axis=1)) / R.dy
        # funcion de Bernoulli
        B = self.g * (h + self.hs) + self.energia_cinetica(u, v)
        du = Fu - (B - np.roll(B, 1, axis=1)) / R.dx_u
        dv = Fv - (B[:-1] - B[1:]) / R.dy
        du = _filtro_fourier(du, R.cos_c, self.cos_ref, R.dlam)
        dv = _filtro_fourier(dv, R.cos_v, self.cos_ref, R.dlam)
        return dh, du, dv

    # ---------------------------------------------------------------
    def integrar(self, h, u, v, dt, pasos, cada=None, al_registrar=None):
        """Leapfrog con el filtro RAW (Williams 2009). Devuelve (h, u, v) al final."""
        ant = (h.copy(), u.copy(), v.copy())
        d = self.tendencias(*ant)
        act = tuple(x + dt * dx for x, dx in zip(ant, d))                       # primer paso: Euler
        if al_registrar is not None and cada and 1 % cada == 0:
            al_registrar(1, *act)
        for n in range(1, pasos):
            d = self.tendencias(*act)
            nue = tuple(x + 2 * dt * dx for x, dx in zip(ant, d))
            corr = tuple(0.5 * RAW_NU * (a - 2 * b + c) for a, b, c in zip(ant, act, nue))
            ant = tuple(b + RAW_ALFA * c for b, c in zip(act, corr))
            act = tuple(x + (RAW_ALFA - 1) * c for x, c in zip(nue, corr))
            if al_registrar is not None and cada and (n + 1) % cada == 0:
                al_registrar(n + 1, *act)
        return act

    # ---------------------------------------------------------------
    def integrales(self, h, u, v):
        """Masa, energia total y enstrofia potencial (integrales discretas sobre la esfera)."""
        R = self.R
        masa = float((h * R.area).sum())
        energia = float(((h * self.energia_cinetica(u, v) + 0.5 * self.g * ((h + self.hs) ** 2 - self.hs ** 2))
                         * R.area).sum())
        he = self.h_esquina(h)
        q = (self.vorticidad(u, v) + self.f_esq) / he
        enstrofia = float((0.5 * q * q * he * R.area_dual).sum())
        return {"masa": masa, "energia": energia, "enstrofia_potencial": enstrofia}


# ================================================================
# CASOS DE PRUEBA DE WILLIAMSON ET AL. (1992) ✅
# ================================================================
A_TIERRA = 6.37122e6          # constantes del conjunto de pruebas (Williamson et al. 1992 ✅)
OMEGA_TIERRA = 7.292e-5
G_TIERRA = 9.80616
DIA = 86400.0


def distancia_gc(lam, phi, lam_c, phi_c):
    return np.arccos(np.clip(np.sin(phi_c) * np.sin(phi) + np.cos(phi_c) * np.cos(phi) * np.cos(lam - lam_c), -1, 1))


def viento_rotado(lam, phi, u0, alfa):
    """Rotacion de solido rigido con el eje inclinado alfa (casos 1 y 2)."""
    u = u0 * (np.cos(phi) * np.cos(alfa) + np.cos(lam) * np.sin(phi) * np.sin(alfa))
    v = -u0 * np.sin(lam) * np.sin(alfa)
    return u, v


def _vientos(R, funcion):
    LU, PU = R.malla("u")
    LV, PV = R.malla("v")
    return funcion(LU, PU)[0], funcion(LV, PV)[1]


def rotar_punto(lam_c, phi_c, alfa, t, u0, a):
    """Posicion en el instante t de un punto arrastrado por la rotacion de solido rigido (caso 1).
    El eje de giro es (-sin alfa, 0, cos alfa) en cartesianas; velocidad angular u0/a."""
    p = np.array([math.cos(phi_c) * math.cos(lam_c), math.cos(phi_c) * math.sin(lam_c), math.sin(phi_c)])
    k = np.array([-math.sin(alfa), 0.0, math.cos(alfa)])
    th = u0 / a * t
    p2 = p * math.cos(th) + np.cross(k, p) * math.sin(th) + k * np.dot(k, p) * (1 - math.cos(th))
    return math.atan2(p2[1], p2[0]), math.asin(max(-1.0, min(1.0, p2[2])))


def campana(R, lam_c, phi_c, h0=1000.0):
    """Campana coseno del caso 1: h = h0/2 (1 + cos(pi r / Rb)) si r < Rb = a/3 (en angulo, 1/3)."""
    LAM, PHI = R.malla("h")
    r = distancia_gc(LAM, PHI, lam_c, phi_c)
    Rb = 1.0 / 3.0
    return np.where(r < Rb, 0.5 * h0 * (1 + np.cos(math.pi * r / Rb)), 0.0)


def caso1(alfa=0.0, a=A_TIERRA, filas=36):
    """Caso 1: advection de una campana coseno por un viento de solido rigido, una vuelta en 12 dias.
    Centro inicial (3pi/2, 0) -> en nuestras longitudes, -pi/2 (Williamson 1992 ✅ via SWEET)."""
    R = Rejilla(a, filas)
    u0 = 2 * math.pi * a / (12 * DIA)
    h = campana(R, -math.pi / 2, 0.0)
    u, v = _vientos(R, lambda L, P: viento_rotado(L, P, u0, alfa))
    return h, u, v


def caso1_exacta(t, alfa=0.0, a=A_TIERRA, filas=36):
    R = Rejilla(a, filas)
    u0 = 2 * math.pi * a / (12 * DIA)
    lc, pc = rotar_punto(-math.pi / 2, 0.0, alfa, t, u0, a)
    return campana(R, lc, pc)


def caso2(alfa=0.0, a=A_TIERRA, omega=OMEGA_TIERRA, g=G_TIERRA, filas=36):
    """Caso 2: flujo zonal estacionario en equilibrio geostrofico no lineal (solucion EXACTA: no cambia).
    u0 = 2 pi a / 12 dias, g h0 = 2,94e4 m2/s2; g h = g h0 - (a Omega u0 + u0^2/2) mu^2, con
    mu = -cos(lambda) cos(phi) sin(alfa) + sin(phi) cos(alfa) (Williamson eq. 95 ✅ via SWEET).
    Usar AguasSomeras(..., alfa_coriolis=alfa)."""
    R = Rejilla(a, filas)
    u0 = 2 * math.pi * a / (12 * DIA)
    LAM, PHI = R.malla("h")
    mu = -np.cos(LAM) * np.cos(PHI) * math.sin(alfa) + np.sin(PHI) * math.cos(alfa)
    h = (2.94e4 - (a * omega * u0 + 0.5 * u0 * u0) * mu * mu) / g
    u, v = _vientos(R, lambda L, P: viento_rotado(L, P, u0, alfa))
    return h, u, v


def caso5(a=A_TIERRA, omega=OMEGA_TIERRA, g=G_TIERRA, filas=36):
    """Caso 5: flujo zonal sobre una montana aislada: hs = 2000 m (1 - r/Rm), Rm = pi/9, centro
    (3pi/2, pi/6) = (-pi/2, pi/6) aqui; h0 = 5960 m, u0 = 20 m/s (Williamson 1992 ✅ via SWEET).
    Sin solucion exacta: se comprueban la conservacion y el comportamiento."""
    R = Rejilla(a, filas)
    u0, h0, hs0, Rm = 20.0, 5960.0, 2000.0, math.pi / 9
    LAM, PHI = R.malla("h")
    r = distancia_gc(LAM, PHI, -math.pi / 2, math.pi / 6)
    hs = np.where(r < Rm, hs0 * (1 - r / Rm), 0.0)
    h = h0 - (a * omega * u0 + 0.5 * u0 * u0) * np.sin(PHI) ** 2 / g - hs
    u, v = _vientos(R, lambda L, P: (u0 * np.cos(P), 0 * P))
    return h, u, v, hs


RH_W = RH_K = 7.848e-6
RH_R, RH_H0 = 4, 8000.0


def caso6(a=A_TIERRA, omega=OMEGA_TIERRA, g=G_TIERRA, filas=36):
    """Caso 6: onda de Rossby-Haurwitz de numero 4: w = K = 7,848e-6 1/s, R = 4, h0 = 8000 m
    (Williamson 1992 ✅). Altura con las funciones A, B, C (Williamson ecs. 146-148, transcritas del
    escaneo; su equilibrio se comprueba numericamente en test_fase6_1.py)."""
    Rg = Rejilla(a, filas)
    w = K = RH_W
    Rn = RH_R

    def vel(lam, phi):
        c = np.cos(phi); s = np.sin(phi)
        u = a * w * c + a * K * c ** (Rn - 1) * (Rn * s * s - c * c) * np.cos(Rn * lam)
        v = -a * K * Rn * c ** (Rn - 1) * s * np.sin(Rn * lam)
        return u, v
    LAM, PHI = Rg.malla("h")
    h = altura_rh(LAM, PHI, a, omega, g)
    u, v = _vientos(Rg, vel)
    return h, u, v


def altura_rh(lam, phi, a=A_TIERRA, omega=OMEGA_TIERRA, g=G_TIERRA):
    w = K = RH_W
    Rn = RH_R
    c = np.cos(phi)
    A = 0.5 * w * (2 * omega + w) * c ** 2 + 0.25 * K * K * c ** (2 * Rn) * (
        (Rn + 1) * c ** 2 + (2 * Rn * Rn - Rn - 2) - 2 * Rn * Rn * c ** -2)
    B = 2 * (omega + w) * K / ((Rn + 1) * (Rn + 2)) * c ** Rn * ((Rn * Rn + 2 * Rn + 2) - (Rn + 1) ** 2 * c ** 2)
    C = 0.25 * K * K * c ** (2 * Rn) * ((Rn + 1) * c ** 2 - (Rn + 2))
    return RH_H0 + a * a * (A + B * np.cos(Rn * lam) + C * np.cos(2 * Rn * lam)) / g


def velocidad_fase_rh(omega=OMEGA_TIERRA):
    """Velocidad angular de la onda de Rossby-Haurwitz en la ecuacion de vorticidad barotropica no
    divergente: nu = (R (3 + R) w - 2 Omega) / ((1 + R)(2 + R)) (Haurwitz 1940; Williamson 1992 ✅).
    En aguas someras es solo aproximada (la divergencia la frena un poco)."""
    return (RH_R * (3 + RH_R) * RH_W - 2 * omega) / ((1 + RH_R) * (2 + RH_R))


def errores_normalizados(x, x_exacto, area):
    """Medidas de error de Williamson et al. (1992) ✅: l1, l2 e l_infinito normalizados."""
    d = x - x_exacto
    l1 = (np.abs(d) * area).sum() / (np.abs(x_exacto) * area).sum()
    l2 = math.sqrt((d * d * area).sum() / (x_exacto * x_exacto * area).sum())
    li = np.abs(d).max() / np.abs(x_exacto).max()
    return float(l1), float(l2), float(li)

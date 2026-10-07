# fase30_multicapa.py -- v3.0: atmosfera de varias capas (en seco).
#
# Diseño completo, fuentes y decisiones: DISENO_V3.0.md. Resumen:
#   - N capas en coordenada sigma (sigma = p / p_superficie), que siguen el
#     relieve: sobre una montaña hay menos aire encima (presion local real).
#   - Infrarrojo "gris" de dos flujos con fuente de Planck lineal en tau;
#     espesor optico tau(p) = TAU0*[F*(p/p0) + (1-F)*(p/p0)^4] (forma de
#     Frierson, Held y Zurita-Gotor 2006, verificada en Isca), calibrado
#     con los flujos sin nubes de Wild et al. (2019) y la Atmosfera Estandar.
#   - La luz absorbida por la atmosfera se reparte en altura con la forma
#     de absorcion del vapor de Isca, (p/p0)^4.
#   - Ajuste convectivo seco que conserva la entalpia (Manabe y Strickler
#     1964) a un gradiente critico de 6,5 K/km (PROVISIONAL hasta la v3.1).
#   - Transporte horizontal: difusion de la energia estatica seca cerca de
#     la superficie, s = cp*T + g*z (opcion C' de DISENO_FASE5A.md), y lo
#     que llega se reparte en la troposfera de forma proporcional a la masa.
#
# Este modulo solo contiene la fisica de la columna, vectorizada sobre todas
# las celdas. La usa simular_fase2b() (fase2b_atmosfera.py) cuando el
# interruptor I10 "atmosfera_multicapa" esta encendido.

import math
import numpy as np

try:                                   # numba (opcional): compila el ajuste convectivo; ~100 veces mas rapido
    from numba import njit
    HAY_NUMBA = True
except ImportError:                    # sin numba se usa la forma cerrada en numpy (exacta, mas lenta)
    HAY_NUMBA = False

from parametros import CONSTANTE_SB, P3N_GRAVEDAD

R_AIRE = 287.05
CP_AIRE = 1004.0
KAPPA = R_AIRE / CP_AIRE
P0 = 1.0e5                     # Pa, presion de referencia al nivel del mar (1 bar, decision de Carlos 01/10/2026)
T_ALTURA_ESCALA = 255.0        # K, temperatura del aire para la presion en superficie de las montañas
                               # (la misma que usa ALTURA_ESCALA en fase2b_atmosfera.py; provisional)

# ---- capas (DISENO_V3.0.md, seccion 1) ----
N_CAPAS_ATM = 20               # PROVISIONAL: se fija con la prueba de convergencia del modelo completo
Z_TOPE_RECETA = 40.0e3         # m: receta de espaciado (no es un dato fisico)
H_RECETA = 7.5e3               # m
ESTIRAMIENTO = 1.6             # capas mas finas cerca del suelo

# ---- infrarrojo (seccion 2): calibrado 04/10/2026, convergido en N ----
TAU0_LW = 2.452
F_LW = 0.240
EXP_VAPOR_LW = 4.0             # (c) forma del vapor, provisional hasta la Fase 5b

# ---- luz absorbida por la atmosfera (seccion 3) ----
EXP_VAPOR_SW = 4.0             # (c) forma de Isca (solar_exponent = 4)

# ---- conveccion (seccion 4) ----
GRADIENTE_CRITICO = 0.0065     # K/m (c) PROVISIONAL: en la v3.1, gradiente adiabatico humedo de P3N

# ---- transporte horizontal (seccion 6) ----
SIGMA_TOPE_TRANSPORTE = 0.25

# v3.1 (optimizacion): infrarrojo con kernels compilados (numba) cuando esta disponible; mismo
# resultado que la version numpy salvo redondeo (diferencia maxima medida ~1e-13 W/m2)
USAR_KERNELS_IR = True   # PROVISIONAL: lo que llega se reparte en las capas con sigma >= este valor


def sigma_seminiveles(n=N_CAPAS_ATM):
    """Seminiveles sigma, de arriba (0) a abajo (1). Longitud n+1."""
    k = np.arange(n + 1) / n
    z = Z_TOPE_RECETA * k ** ESTIRAMIENTO
    s = np.exp(-z / H_RECETA)
    s[-1] = 0.0
    return s[::-1].copy()


def presion_superficie(altitud_metros, gravedad=P3N_GRAVEDAD):
    """Presion en superficie de cada celda (Pa), en equilibrio hidrostatico."""
    return P0 * np.exp(-gravedad * np.maximum(altitud_metros, 0.0) / (R_AIRE * T_ALTURA_ESCALA))


def _pav_columnas(th, w):
    """Regresion isotonica ponderada (theta* no decreciente desde abajo) de cada
    columna, por el algoritmo clasico de "mezclar vecinos inestables" (pool
    adjacent violators). th, w: (n, M), indice 0 = capa mas baja. Modifica th."""
    n, m = th.shape
    val = np.empty(n); pes = np.empty(n); lon = np.empty(n, dtype=np.int64)
    for c in range(m):
        nb = 0
        for k in range(n):
            val[nb] = th[k, c]; pes[nb] = w[k, c]; lon[nb] = 1; nb += 1
            while nb > 1 and val[nb - 2] > val[nb - 1]:
                pt = pes[nb - 2] + pes[nb - 1]
                val[nb - 2] = (val[nb - 2] * pes[nb - 2] + val[nb - 1] * pes[nb - 1]) / pt
                pes[nb - 2] = pt; lon[nb - 2] += lon[nb - 1]; nb -= 1
        k = 0
        for b in range(nb):
            for _ in range(lon[b]):
                th[k, c] = val[b]; k += 1
    return th


def _pav_columnas_vapor(th, w, qv, wq):
    """Como _pav_columnas, y ademas mezcla el vapor qv (ponderado con la masa wq) dentro de cada
    tramo de capas que el ajuste ha mezclado: la conveccion seca mezcla el aire entero, con su vapor
    (v3.1). Modifica th y qv."""
    n, m = th.shape
    val = np.empty(n); pes = np.empty(n); lon = np.empty(n, dtype=np.int64)
    for c in range(m):
        nb = 0
        for k in range(n):
            val[nb] = th[k, c]; pes[nb] = w[k, c]; lon[nb] = 1; nb += 1
            while nb > 1 and val[nb - 2] > val[nb - 1]:
                pt = pes[nb - 2] + pes[nb - 1]
                val[nb - 2] = (val[nb - 2] * pes[nb - 2] + val[nb - 1] * pes[nb - 1]) / pt
                pes[nb - 2] = pt; lon[nb - 2] += lon[nb - 1]; nb -= 1
        k = 0
        for b in range(nb):
            if lon[b] > 1:
                sq = 0.0; sw = 0.0
                for j in range(k, k + lon[b]):
                    sq += qv[j, c] * wq[j, c]; sw += wq[j, c]
                for j in range(k, k + lon[b]):
                    qv[j, c] = sq / sw
            for _ in range(lon[b]):
                th[k, c] = val[b]; k += 1
    return th, qv


def _ir_bajada(T, T_aire, peso, tr, g1, sb):
    """Kernel del infrarrojo descendente (v3.1, optimizacion): mismas operaciones, en el mismo orden,
    que la version numpy de Columna.infrarrojo_bajada, columna a columna. Arrays (n, M)."""
    n, m = T.shape
    Bh = np.empty((n + 1, m)); D = np.empty((n + 1, m))
    for c in range(m):
        t = T[0, c]
        Bh[0, c] = sb * t ** 4
        for k in range(1, n):
            t = T[k - 1, c] + peso[k - 1, c] * (T[k, c] - T[k - 1, c])
            Bh[k, c] = sb * t ** 4
        t = T_aire[c]
        Bh[n, c] = sb * t ** 4
        D[0, c] = 0.0
        for k in range(n):
            D[k + 1, c] = D[k, c] * tr[k, c] + Bh[k + 1, c] * (1 - tr[k, c]) - (Bh[k + 1, c] - Bh[k, c]) * g1[k, c]
    return D, Bh


def _ir_subida(Bh, emision, tr, g1):
    n1, m = Bh.shape
    n = n1 - 1
    U = np.empty((n1, m))
    for c in range(m):
        U[n, c] = emision[c]
        for k in range(n - 1, -1, -1):
            U[k, c] = U[k + 1, c] * tr[k, c] + Bh[k, c] * (1 - tr[k, c]) + (Bh[k + 1, c] - Bh[k, c]) * g1[k, c]
    return U


if HAY_NUMBA:
    _pav_columnas = njit(cache=True)(_pav_columnas)
    _ir_bajada = njit(cache=True)(_ir_bajada)
    _ir_subida = njit(cache=True)(_ir_subida)
    _pav_columnas_vapor = njit(cache=True)(_pav_columnas_vapor)


def alfa_simmons_burridge(sh):
    """alfa_k de Simmons y Burridge (1981) en sigma pura, igual que fase6_nucleo.NucleoSeco: alfa_0 = ln 2
    (tope a p = 0); alfa_k = 1 - sigma_{k-1/2}/dsigma_k ln(sigma_{k+1/2}/sigma_{k-1/2})."""
    n = len(sh) - 1
    dsig = np.diff(sh)
    dln = np.zeros(n)
    dln[1:] = np.log(sh[2:] / sh[1:-1])
    alfa = np.empty(n)
    alfa[0] = math.log(2.0)
    alfa[1:] = 1.0 - sh[1:-1] / dsig[1:] * dln[1:]
    return alfa


class Columna:
    """Geometria vertical de todas las celdas. Fija en el tiempo salvo que se llame a actualizar_ps
    (v3.1-pre9: con el nucleo dinamico, la presion en superficie cambia en cada paso)."""

    def __init__(self, altitud_metros, gravedad=P3N_GRAVEDAD, n=N_CAPAS_ATM, gradiente=GRADIENTE_CRITICO,
                 presion_capa="media"):
        # gradiente (K/m): el del ajuste convectivo SECO. v3.0: 6,5 K/km provisional; v3.1 con la
        # conveccion humeda (I13): el adiabatico seco g/cp de cada planeta (fisica, a)
        # presion_capa (v3.1-pre9, DISENO_FASE6_3.md §6.2, decision 1.2): "media" = media aritmetica de los
        # seminiveles (la de la v3.0/v3.1, por defecto: sin cambios); "sb81" = la de Simmons y Burridge
        # (1981), la del nucleo dinamico: ln p_k = ln p_{k+1/2} - alfa_k.
        if presion_capa not in ("media", "sb81"):
            raise ValueError("presion_capa debe ser 'media' o 'sb81'")
        self.n = n
        self.gradiente = gradiente
        self.g = gravedad
        self.forma = altitud_metros.shape
        self.presion_capa = presion_capa
        self.actualizar_ps(presion_superficie(altitud_metros, gravedad))

    def actualizar_ps(self, ps):
        """Recalcula TODO lo que depende de la presion en superficie ps (F, C). Con el nucleo dinamico
        (I16) se llama en cada paso de la fisica con la p_s del nucleo."""
        n, gravedad, gradiente = self.n, self.g, self.gradiente
        self.ps = ps                                                             # (F, C)
        sh = sigma_seminiveles(n)
        self.ph = sh[:, None, None] * self.ps[None]                              # (n+1, F, C)
        if self.presion_capa == "media":
            self.pm = 0.5 * (self.ph[1:] + self.ph[:-1])                         # (n, F, C)
        else:
            self.pm = (sh[1:] * np.exp(-alfa_simmons_burridge(sh)))[:, None, None] * self.ps[None]
        self.dp = np.diff(self.ph, axis=0)                                       # (n, F, C)
        self.cap = self.dp / gravedad * CP_AIRE                                  # J/m2/K de cada capa
        self.sigma_media = 0.5 * (sh[1:] + sh[:-1])
        # infrarrojo: espesor optico en los seminiveles (presion ABSOLUTA)
        x = self.ph / P0
        self.tau = TAU0_LW * (F_LW * x + (1 - F_LW) * x ** EXP_VAPOR_LW)
        dt = np.maximum(np.diff(self.tau, axis=0), 1e-12)
        self.trans = np.exp(-dt)
        self.g1 = (1 - self.trans * (1 + dt)) / dt
        # interpolacion de la temperatura a los seminiveles interiores (en ln p)
        lp = np.log(self.pm)
        lh = np.log(np.maximum(self.ph[1:-1], 1e-6))
        self.peso_interp = (lh - lp[:-1]) / (lp[1:] - lp[:-1])                    # (n-1, F, C)
        # luz absorbida por la atmosfera: reparto vertical (suma 1 en cada celda)
        w = np.diff((self.ph / P0) ** EXP_VAPOR_SW, axis=0)
        self.peso_sw = w / w.sum(axis=0, keepdims=True)
        # ajuste convectivo: T_arriba >= T_abajo * (p_arriba/p_abajo)^(R*Gamma/g)
        # PI_k: temperatura de un perfil critico relativa a la capa mas baja (PI = 1 abajo)
        self.pi_ajuste = (self.pm / self.pm[-1][None]) ** (R_AIRE * gradiente / gravedad)
        # extrapolacion adiabatica seca de la capa mas baja a la superficie
        self.factor_superficie = (self.ps / self.pm[-1]) ** KAPPA
        # transporte horizontal: capas que reciben lo que llega
        self.capas_transporte = self.sigma_media >= SIGMA_TOPE_TRANSPORTE
        self.cap_transporte = self.cap[self.capas_transporte].sum(axis=0)          # (F, C)
        # copias contiguas (n, celdas) para los kernels compilados del infrarrojo
        self._peso_2d = np.ascontiguousarray(self.peso_interp.reshape(n - 1, -1)) if n > 1 else np.zeros((1, 1))
        self._tr_2d = np.ascontiguousarray(self.trans.reshape(n, -1))
        self._g1_2d = np.ascontiguousarray(self.g1.reshape(n, -1))

    # ------------------------------------------------------------------
    def actualizar_tau_vapor(self, q, a, b):
        """v3.1 (I15, prototipo): espesor optico infrarrojo de cada capa segun el vapor de la propia
        capa, d(tau) = (a + b*q) * dp / P0 (forma de Byrne y O'Gorman 2013, la del codigo de Isca),
        con a (aire seco: CO2 y demas) y b (vapor) calibrados para M3N (fase31_agua). Recalcula las
        transmisiones de cada capa; el resto del esquema (dos flujos, fuente lineal) no cambia."""
        dt = np.maximum((a + b * np.maximum(q, 0.0)) * self.dp / P0, 1e-12)
        self.trans = np.exp(-dt)
        self.g1 = (1 - self.trans * (1 + dt)) / dt
        n = self.n
        self._tr_2d = np.ascontiguousarray(self.trans.reshape(n, -1))
        self._g1_2d = np.ascontiguousarray(self.g1.reshape(n, -1))

    def temperatura_seminiveles(self, T, T_aire_superficie):
        Th = np.empty((self.n + 1,) + self.forma)
        Th[0] = T[0]
        Th[1:-1] = T[:-1] + self.peso_interp * (T[1:] - T[:-1])
        Th[-1] = T_aire_superficie
        return Th

    def infrarrojo_bajada(self, T, T_aire_superficie):
        """Flujo descendente en los seminiveles. D[-1] = infrarrojo que llega al suelo."""
        if HAY_NUMBA and USAR_KERNELS_IR:
            n = self.n
            D, Bh = _ir_bajada(np.ascontiguousarray(T.reshape(n, -1)), np.ascontiguousarray(T_aire_superficie.reshape(-1)),
                               self._peso_2d, self._tr_2d, self._g1_2d, CONSTANTE_SB)
            return D.reshape((n + 1,) + self.forma), Bh.reshape((n + 1,) + self.forma)
        Bh = CONSTANTE_SB * self.temperatura_seminiveles(T, T_aire_superficie) ** 4
        D = np.zeros_like(Bh)
        tr, g1 = self.trans, self.g1
        for k in range(self.n):
            D[k + 1] = D[k] * tr[k] + Bh[k + 1] * (1 - tr[k]) - (Bh[k + 1] - Bh[k]) * g1[k]
        return D, Bh

    def infrarrojo_subida(self, Bh, emision_superficie):
        """Flujo ascendente en los seminiveles. U[0] = infrarrojo que sale al espacio."""
        if HAY_NUMBA and USAR_KERNELS_IR:
            n = self.n
            U = _ir_subida(np.ascontiguousarray(Bh.reshape(n + 1, -1)),
                           np.ascontiguousarray(np.broadcast_to(emision_superficie, self.forma).reshape(-1)),
                           self._tr_2d, self._g1_2d)
            return U.reshape((n + 1,) + self.forma)
        U = np.empty_like(Bh)
        U[-1] = emision_superficie
        tr, g1 = self.trans, self.g1
        for k in range(self.n - 1, -1, -1):
            U[k] = U[k + 1] * tr[k] + Bh[k] * (1 - tr[k]) + (Bh[k + 1] - Bh[k]) * g1[k]
        return U

    @staticmethod
    def calentamiento_infrarrojo(U, D):
        """W/m2 que gana cada capa por el infrarrojo."""
        return (U[1:] - U[:-1]) + (D[:-1] - D[1:])

    def ajuste_convectivo(self, T):
        """Ajuste convectivo seco EXACTO, que conserva sum(cap*T). Devuelve una copia.

        Con theta* = T / PI_k (PI_k: perfil critico relativo a la capa mas
        baja), el perfil critico es theta* constante y "estable" es theta* no
        decreciente hacia arriba. Mezclar capas = sustituirlas por la media de
        theta* ponderada con w_k = cap_k * PI_k, que conserva la energia. El
        resultado del ajuste completo es la REGRESION ISOTONICA ponderada de
        theta* (la misma que da el algoritmo clasico de ir mezclando capas
        inestables hasta que no queda ninguna), y tiene forma cerrada:
            theta*_i = max_{j<=i} min_{l>=i} media(j..l)
        (i, j, l contados desde abajo). Se calcula a la vez en todas las
        celdas con sumas acumuladas: exacto en una sola pasada."""
        n = self.n
        thc = T / self.pi_ajuste
        inestable = (thc[:-1] < thc[1:]).any(axis=0)       # alguna capa con theta* menor que la de debajo
        T = T.copy()
        if not inestable.any():
            return T
        if HAY_NUMBA:
            th = np.ascontiguousarray(thc[::-1][:, inestable])
            w = np.ascontiguousarray((self.cap * self.pi_ajuste)[::-1][:, inestable])
            T[:, inestable] = _pav_columnas(th, w)[::-1] * self.pi_ajuste[:, inestable]
            return T
        th = thc[::-1]
        w = (self.cap * self.pi_ajuste)[::-1]
        cero = np.zeros((1,) + self.forma)
        Sw = np.concatenate([cero, np.cumsum(w, axis=0)], axis=0)        # (n+1, F, C)
        Swt = np.concatenate([cero, np.cumsum(w * th, axis=0)], axis=0)
        # media del tramo j..l (j <= l): (Swt[l+1]-Swt[j]) / (Sw[l+1]-Sw[j])
        with np.errstate(invalid='ignore', divide='ignore'):
            M = (Swt[None, 1:] - Swt[:-1, None]) / (Sw[None, 1:] - Sw[:-1, None])   # (j, l, F, C)
        jj = np.arange(n)[:, None, None, None]
        ll = np.arange(n)[None, :, None, None]
        M = np.where(ll >= jj, M, np.inf)
        # A[j, i] = min_{l >= i} M[j, l]
        A = np.minimum.accumulate(M[:, ::-1], axis=1)[:, ::-1]
        A = np.where(ll >= jj, A, -np.inf)                # solo j <= i
        th_aj = A.max(axis=0)                             # max_{j <= i}
        return np.where(inestable[None], th_aj[::-1] * self.pi_ajuste, T)

    def ajuste_convectivo_vapor(self, T, q):
        """v3.1: ajuste convectivo seco (el mismo que ajuste_convectivo) que ademas mezcla el vapor
        en los tramos mezclados, conservando su masa. Devuelve copias (T, q)."""
        thc = T / self.pi_ajuste
        inestable = (thc[:-1] < thc[1:]).any(axis=0)
        T = T.copy(); q = q.copy()
        if not inestable.any():
            return T, q
        th = np.ascontiguousarray(thc[::-1][:, inestable])
        w = np.ascontiguousarray((self.cap * self.pi_ajuste)[::-1][:, inestable])
        qv = np.ascontiguousarray(q[::-1][:, inestable])
        wq = np.ascontiguousarray(self.dp[::-1][:, inestable])
        th, qv = _pav_columnas_vapor(th, w, qv, wq)
        T[:, inestable] = th[::-1] * self.pi_ajuste[:, inestable]
        q[:, inestable] = qv[::-1]
        return T, q

    def aire_superficie(self, T):
        """Temperatura del aire junto a la superficie: extrapolacion adiabatica
        seca desde la capa mas baja (equivale a conservar s = cp*T + g*z)."""
        return T[-1] * self.factor_superficie

    def perfil_inicial(self, T_superficie):
        # estado inicial: gradiente critico desde la superficie, con un suelo de 200 K
        # en lo alto (la estratosfera se ajusta sola en pocas semanas de simulacion)
        T = T_superficie[None] * (self.pm / self.ps[None]) ** (R_AIRE * min(self.gradiente, GRADIENTE_CRITICO) / self.g)
        return np.maximum(T, 200.0)

    def media_masa(self, T, mascara_capas=None):
        if mascara_capas is None:
            return (T * self.cap).sum(axis=0) / self.cap.sum(axis=0)
        return (T[mascara_capas] * self.cap[mascara_capas]).sum(axis=0) / self.cap[mascara_capas].sum(axis=0)


# ================================================================
# DIAGNOSTICO: TRANSPORTE MERIDIONAL DE CALOR (v3.0)
# ================================================================
def transporte_meridional(convergencia, radio):
    """
    Transporte de calor hacia el norte (PW) a traves de cada borde entre
    filas, a partir de la convergencia media anual del transporte (W/m2,
    FILAS x COLUMNAS; positiva = la celda recibe). Se integra desde el
    polo norte: lo que cruza el borde sur de la fila i (hacia el norte) es lo que
    han recibido todas las filas de i hacia el norte. El area de cada
    celda es la del operador de difusion (radio^2 cos(lat) dlat dlon), con
    la que la difusion conserva exactamente la energia: asi el transporte
    en el polo sur sale cero (comprobacion). Devuelve (latitudes de los
    bordes, de 90 - 5 a -90 + 5 grados; transporte en PW).
    """
    from rejilla import LATITUDES_GRADOS, GRADOS_POR_FILA, GRADOS_POR_COL
    area = radio ** 2 * np.cos(np.radians(LATITUDES_GRADOS)) * np.radians(GRADOS_POR_FILA) * np.radians(GRADOS_POR_COL)
    por_fila = (convergencia * area[:, None]).sum(axis=1)          # W recibidos por cada fila
    hacia_norte = np.cumsum(por_fila)[:-1] / 1e15                  # el borde sur de la fila i lleva hacia el norte lo que recibieron 0..i
    bordes = LATITUDES_GRADOS[:-1] - GRADOS_POR_FILA / 2
    return bordes, hacia_norte

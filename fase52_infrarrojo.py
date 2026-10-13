# fase52_infrarrojo.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §7, T5, T9, T19 y T20): el INFRARROJO
# POR BANDAS de una columna. No se usa todavia en el modelo.
#
# Que hay:
#   - planck_banda: la emision de cuerpo negro de una banda, pi ∫ B_nu dnu (W/m^2).
#   - dos_flujos: los flujos de una columna para un conjunto de "terminos" pseudo-monocromaticos (cada punto g de cada
#     banda, o cada numero de onda en el calculo linea a linea), con el MISMO esquema de dos flujos que ya tiene M3N
#     (T5; DISENO_FASE_2.3.md §2, fase30_multicapa.infrarrojo_bajada/subida): fuente de Planck lineal en tau.
#   - combinar_rorr: el solapamiento de varios gases por "random overlap with resorting and rebinning" (Lacis y
#     Oinas 1991 ⚠️; formulas como las describen Amundsen et al. 2017, A&A 598, A97, §2 ✅ leido el 13/10/2026).
#   - columna_lbl / columna_ck: los flujos de una banda linea a linea (referencia) y con tablas k.
#
# T19 (provisional): factor de difusividad 1,66 (Elsasser 1942 ⚠️): el flujo difuso de una capa de espesor optico
#     vertical tau se transmite como exp(-1,66 tau). Es la aproximacion habitual en el infrarrojo (RRTMG la usa en
#     casi todas sus bandas ⚠️). Su error frente al calculo con varias direcciones se medira en la validacion (§4).
# T20 (provisional): la emision que cae fuera de las bandas (por debajo de 10 cm^-1 y por encima de 3250 cm^-1) se
#     suma a la primera y a la ultima banda, para que la suma de las bandas sea exactamente sigma T^4 y el balance de
#     energia de la superficie y de las capas cierre. Medido: fuera de 10-3250 cm^-1 queda menos del 0,05 % de sigma
#     T^4 entre 150 y 330 K (test_fase52_infrarrojo).

import numpy as np

C1 = 1.191042972e-8          # 2 h c^2, W m^-2 sr^-1 (cm^-1)^-4 (CODATA 2018)
C2 = 1.438776877             # h c / k, cm K (CODATA 2018)
SIGMA = 5.670374419e-8       # W m^-2 K^-4 (CODATA 2018)
DIFUSIVIDAD = 1.66           # T19

_XG, _WG = np.polynomial.legendre.leggauss(64)


def planck_integral(nu1, nu2, T):
    """pi ∫_nu1^nu2 B_nu(T) dnu en W/m^2 (nu en cm^-1), por Gauss-Legendre de 64 puntos en tramos de 50 cm^-1 como
    mucho (error < 1e-12 relativo: test). T puede ser un array."""
    T = np.asarray(T, dtype=float)
    out = np.zeros_like(T)
    n_tramos = max(int(np.ceil((nu2 - nu1) / 50.0)), 1)
    bordes = np.linspace(nu1, nu2, n_tramos + 1)
    for a, b in zip(bordes[:-1], bordes[1:]):
        nu = 0.5 * (b - a) * _XG + 0.5 * (a + b)
        x = C2 * nu[:, None] / T.reshape(1, -1)
        B = C1 * nu[:, None] ** 3 / np.expm1(x)
        out += (0.5 * (b - a) * (_WG[:, None] * B).sum(axis=0)).reshape(T.shape)
    return np.pi * out


def planck_cola_alta(nu1, T):
    """pi ∫_nu1^inf B_nu(T) dnu (W/m^2): la integral hasta donde la exponencial ya es < 1e-30 del maximo."""
    T = np.asarray(T, dtype=float)
    nu_max = nu1 + 80.0 * np.max(T) / C2
    return planck_integral(nu1, nu_max, T)


def planck_banda(bandas, T, colas=True):
    """(bandas, *T.shape): emision de cada banda; con colas (T20), la de fuera de las bandas se suma a la primera y a
    la ultima (las bandas tienen que ir seguidas y ordenadas)."""
    T = np.asarray(T, dtype=float)
    out = np.array([planck_integral(a, b, T) for a, b in bandas])
    if colas:
        for (a, b), (c, d) in zip(bandas[:-1], bandas[1:]):
            if b != c:
                raise ValueError("para sumar las colas, las bandas tienen que ir seguidas")
        out[0] += planck_integral(1e-6, bandas[0][0], T)
        out[-1] += planck_cola_alta(bandas[-1][1], T)
    return out


def _g1(x):
    """(1 - e^-x (1 + x)) / x, con su serie para x pequeño (evita la cancelacion)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    p = x < 1e-2
    xp = x[p]
    out[p] = xp / 2 - xp ** 2 / 3 + xp ** 3 / 8 - xp ** 4 / 30 + xp ** 5 / 144
    xg = x[~p]
    out[~p] = (1.0 - np.exp(-xg) * (1.0 + xg)) / xg
    return out


def dos_flujos(dtau, pesos, Bh, Bs, emisividad=1.0):
    """Flujos de infrarrojo (W/m^2) de un conjunto de terminos pseudo-monocromaticos.
    dtau: (terminos, n) espesor optico VERTICAL de cada capa (0 = la de arriba); pesos: (terminos,); Bh: (terminos, n+1)
    o (n+1,) emision de cuerpo negro (W/m^2) de cada termino en los seminiveles (la de abajo, con la temperatura del
    aire en superficie, como fase30_multicapa); Bs: emision de la superficie (escalar o (terminos,)); emisividad gris.
    Devuelve (U, D), (n+1,): suma pesada de los flujos ascendente y descendente en los seminiveles."""
    dtau = np.atleast_2d(np.asarray(dtau, dtype=float))
    m, n = dtau.shape
    x = DIFUSIVIDAD * np.maximum(dtau, 0.0)
    tr = np.exp(-x)
    g1 = _g1(x)
    Bh = np.broadcast_to(np.asarray(Bh, dtype=float), (m, n + 1))
    Bs = np.broadcast_to(np.asarray(Bs, dtype=float), (m,))
    D = np.zeros((m, n + 1))
    for k in range(n):
        D[:, k + 1] = D[:, k] * tr[:, k] + Bh[:, k + 1] * (1 - tr[:, k]) - (Bh[:, k + 1] - Bh[:, k]) * g1[:, k]
    U = np.empty((m, n + 1))
    U[:, n] = emisividad * Bs + (1.0 - emisividad) * D[:, n]
    for k in range(n - 1, -1, -1):
        U[:, k] = U[:, k + 1] * tr[:, k] + Bh[:, k] * (1 - tr[:, k]) + (Bh[:, k + 1] - Bh[:, k]) * g1[:, k]
    w = np.asarray(pesos, dtype=float)[:, None]
    return (w * U).sum(axis=0), (w * D).sum(axis=0)


def calentamiento(U, D):
    """W/m^2 que gana cada capa (como fase30_multicapa.calentamiento_infrarrojo)."""
    return (U[1:] - U[:-1]) + (D[:-1] - D[1:])


# ------------------------------------------------------------------------------------------------ solapamiento
def _reagrupar(tau, w, bordes):
    """Reparte los terminos (tau, w), ordenados por tau, en los intervalos de g 'bordes' y devuelve la media de tau
    pesada en cada uno (un termino que cae entre dos intervalos se parte, como en Amundsen et al. 2017 §2.4)."""
    orden = np.argsort(tau, kind="stable")
    t, p = tau[orden], w[orden]
    G = np.concatenate([[0.0], np.cumsum(p)])
    G /= G[-1]
    I = np.concatenate([[0.0], np.cumsum(t * p / p.sum())])     # integral de tau(g) dg
    Ib = np.interp(bordes, G, I)                                # tau(g) es constante a trozos: I es lineal a trozos
    return np.diff(Ib) / np.diff(bordes)


def combinar_rorr(taus, pesos):
    """Solapamiento aleatorio con reordenacion y reagrupamiento, capa a capa.
    taus: lista (un elemento por gas) de (n, puntos g) espesores opticos de cada capa y punto g; pesos: (puntos g,),
    los mismos para todos los gases y para el resultado. Devuelve (n, puntos g)."""
    pesos = np.asarray(pesos, dtype=float)
    bordes = np.concatenate([[0.0], np.cumsum(pesos)])
    bordes /= bordes[-1]
    acum = np.array(taus[0], dtype=float)
    ww = np.outer(pesos, pesos).ravel()
    for t in taus[1:]:
        t = np.asarray(t, dtype=float)
        nuevo = np.empty_like(acum)
        for k in range(acum.shape[0]):
            comb = (acum[k][:, None] + t[k][None, :]).ravel()
            nuevo[k] = _reagrupar(comb, ww, bordes)
        acum = nuevo
    return acum


def tabla_mapa(tau_ref_columna, dnu, pesos):
    """Agrupamiento FIJO de los numeros de onda de una banda (metodo del "mapa espectral", como ecCKD: Hogan y
    Matricardi 2022, JAMES 14, e2022MS003033 ✅t, que ordenan por la altura de maximo enfriamiento): cada numero de
    onda se ordena por la presion a la que su espesor optico desde arriba en una columna de referencia llega a 1
    (tau_ref_columna: (n+1, n_nu) espesor optico acumulado desde arriba en los seminiveles), y se reparte en los
    intervalos de g de 'pesos' por su anchura dnu. Devuelve el indice de grupo (0..puntos g-1) de cada numero de onda."""
    tau = np.asarray(tau_ref_columna, dtype=float)
    n1, m = tau.shape
    llega = tau >= 1.0
    # nivel (fraccionario) donde tau = 1: los que nunca llegan, por debajo del suelo segun su tau total
    idx = np.where(llega.any(axis=0), llega.argmax(axis=0), n1)
    clave = np.empty(m)
    dentro = idx < n1
    j = idx[dentro]
    t0 = np.where(j > 0, tau[np.maximum(j - 1, 0), np.nonzero(dentro)[0]], 0.0)
    t1 = tau[j, np.nonzero(dentro)[0]]
    clave[dentro] = (j - 1) + (1.0 - t0) / np.maximum(t1 - t0, 1e-300)
    clave[~dentro] = n1 + (1.0 - tau[-1, ~dentro])              # mas transparentes: detras, por tau total
    orden = np.argsort(-clave, kind="stable")                   # de mas transparente (g = 0) a mas opaco (g = 1)
    G = np.cumsum(np.asarray(dnu, dtype=float)[orden])
    G = (G - 0.5 * np.asarray(dnu, dtype=float)[orden]) / G[-1]   # g del centro de cada numero de onda
    bordes = np.concatenate([[0.0], np.cumsum(pesos)])
    bordes /= bordes[-1]
    grupo = np.empty(m, dtype=int)
    grupo[orden] = np.clip(np.searchsorted(bordes, G, side="right") - 1, 0, len(pesos) - 1)
    return grupo


def media_por_grupo(k_capas, dnu, grupo, n_grupos, tipo="aritmetica"):
    """(n, grupos): media de k (n, n_nu) en cada grupo, pesada por dnu; tipo "aritmetica" o "logaritmica"."""
    k = np.asarray(k_capas, dtype=float)
    w = np.asarray(dnu, dtype=float)
    out = np.zeros((k.shape[0], n_grupos))
    for j in range(n_grupos):
        s = grupo == j
        if not s.any():
            continue
        ws = w[s] / w[s].sum()
        if tipo == "aritmetica":
            out[:, j] = (k[:, s] * ws).sum(axis=1)
        else:
            out[:, j] = np.exp((np.log(np.maximum(k[:, s], 1e-300)) * ws).sum(axis=1))
    return out


def media_transmision_por_grupo(tau_capas, peso, grupo, n_grupos):
    """(n, grupos): espesor optico de cada capa y grupo que da la misma transmision difusa media, pesada con 'peso'
    ((n, n_nu) o (n_nu,), p. ej. Planck por anchura): tau = -ln(Σ p e^(-1,66 tau) / Σ p) / 1,66. Es el promedio por
    defecto de ecCKD ("transmission", average_optical_depth.cpp, leido el 13/10/2026 ✅; Apache 2.0)."""
    t = np.asarray(tau_capas, dtype=float)
    p = np.broadcast_to(np.asarray(peso, dtype=float), t.shape)
    out = np.zeros((t.shape[0], n_grupos))
    for j in range(n_grupos):
        s = grupo == j
        if not s.any():
            continue
        ps = p[:, s]
        # 1 - transmision media, con expm1/log1p: exacto tambien con espesores muy pequeños (limite: media aritmetica)
        # (si absorbe mas de la mitad no hay cancelacion: la transmision directa, con un tope de 1e-300 -> tau ~ 416)
        a = -(ps * np.expm1(-DIFUSIVIDAD * t[:, s])).sum(axis=1) / ps.sum(axis=1)
        tr = (ps * np.exp(-DIFUSIVIDAD * t[:, s])).sum(axis=1) / ps.sum(axis=1)
        out[:, j] = np.where(a < 0.5, -np.log1p(-np.minimum(a, 0.5)), -np.log(np.maximum(tr, 1e-300))) / DIFUSIVIDAD
    return out


# ------------------------------------------------------------------------------------------- con las tablas (T9)
class TablasInfrarrojo:
    """Las tablas k del infrarrojo con mapa espectral fijo (herramientas/generar_tablas_k.py, 3.23.0) y el calculo de
    los flujos de una columna con ellas: espesor optico de cada grupo = suma de los gases (u k), emision de cada grupo
    = su fraccion de Planck por la emision de la banda (T20).
    Interpolacion en la rejilla de las tablas: lineal en ln p y en T, de ln k ("log", por defecto) o de k ("lineal");
    en la fraccion de vapor, en su logaritmo (con el mismo k o ln k). Fuera de la rejilla, el extremo."""

    def __init__(self, carpeta, gases=None, interpolacion="log"):
        import os
        m = np.load(os.path.join(carpeta, "mapas_infrarrojo.npz"))
        self.bandas = [tuple(b) for b in m["bandas"]]
        self.pesos = m["pesos"]                                   # (bandas, grupos)
        self.T_fr, self.fr = m["T_fr"], m["fraccion_planck"]      # (bandas, grupos, T)
        self.tablas = {}
        for f in sorted(os.listdir(carpeta)):
            if f.startswith("tabla_k_") and f.endswith(".npz"):
                gas = f[8:-4]
                if gases is None or gas in gases:
                    d = np.load(os.path.join(carpeta, f))
                    self.tablas[gas] = {"k": d["k"], "lp": np.log(d["p_hPa"]), "T": d["T_K"], "vmr": d["vmr"]}
        if interpolacion not in ("log", "lineal"):
            raise ValueError("interpolacion: 'log' o 'lineal'")
        self.interpolacion = interpolacion

    @staticmethod
    def _indice(x, nodos):
        """(i, f) con nodos[i] <= x <= nodos[i+1] para nodos crecientes o decrecientes; fuera, el extremo."""
        nodos = np.asarray(nodos, dtype=float)
        dec = nodos[0] > nodos[-1]
        n = nodos[::-1] if dec else nodos
        x = np.clip(np.asarray(x, dtype=float), n[0], n[-1])
        i = np.clip(np.searchsorted(n, x) - 1, 0, len(n) - 2)
        f = (x - n[i]) / (n[i + 1] - n[i])
        if dec:
            i, f = len(n) - 2 - i, 1.0 - f
        return i, f

    def k(self, gas, p_hPa, T, vmr_h2o=0.0):
        """(bandas, grupos, capas): k (cm^2/molecula) del gas en cada capa."""
        t = self.tablas[gas]
        p_hPa = np.atleast_1d(np.asarray(p_hPa, dtype=float))
        T = np.broadcast_to(np.asarray(T, dtype=float), p_hPa.shape)
        vmr_h2o = np.broadcast_to(np.asarray(vmr_h2o, dtype=float), p_hPa.shape)
        ip, fp = self._indice(np.log(p_hPa), t["lp"])
        it, ft = self._indice(T, t["T"])
        if len(t["vmr"]) > 1:            # en ln de la fraccion molar (rejilla logaritmica de generar_tablas_k)
            iv, fv = self._indice(np.log(np.maximum(vmr_h2o, 1e-12)), np.log(t["vmr"]))
        else:
            iv, fv = np.zeros_like(ip), np.zeros(np.shape(ip))
        K = t["k"]                                                # (bandas, p, T, vmr, grupos)
        out = 0.0
        for a, wa in ((0, 1 - fp), (1, fp)):
            for b, wb in ((0, 1 - ft), (1, ft)):
                if self.interpolacion == "log":
                    def val(c):
                        return np.log(np.maximum(K[:, ip + a, it + b, np.minimum(iv + c, K.shape[3] - 1), :], 1e-300))
                else:
                    def val(c):
                        return K[:, ip + a, it + b, np.minimum(iv + c, K.shape[3] - 1), :]
                v = (1 - fv)[None, :, None] * val(0) + fv[None, :, None] * val(1)       # (bandas, capas, grupos)
                out = out + (wa * wb)[None, :, None] * v
        if self.interpolacion == "log":
            out = np.exp(out)
            out[out < 1e-299] = 0.0
        return np.transpose(out, (0, 2, 1))

    def espesores(self, p_hPa, T, u, vmr_h2o):
        """(bandas, grupos, capas): espesor optico vertical; u: {gas: (capas,) moleculas/cm^2}."""
        tau = 0.0
        for gas, ug in u.items():
            if gas in self.tablas:
                tau = tau + self.k(gas, p_hPa, T, vmr_h2o) * np.asarray(ug)[None, None, :]
        return tau

    def fraccion(self, ib, T):
        """(grupos, len(T)): fraccion de Planck de cada grupo de la banda ib a las temperaturas T."""
        T = np.atleast_1d(np.asarray(T, dtype=float))
        return np.array([np.interp(T, self.T_fr, self.fr[ib, j]) for j in range(self.fr.shape[1])])

    def flujos(self, p_hPa, T, Th, Ts, u, vmr_h2o, emisividad=1.0):
        """Flujos de infrarrojo (U, D) en los seminiveles (W/m^2), suma de las bandas de las tablas.
        Th: temperatura en los seminiveles (la de abajo, la del aire en superficie); Ts: la del suelo."""
        tau = self.espesores(p_hPa, T, u, vmr_h2o)
        Bh = planck_banda(self.bandas, Th) if len(self.bandas) == 16 else \
            np.array([planck_integral(a, b, Th) for a, b in self.bandas])
        Bs = planck_banda(self.bandas, Ts) if len(self.bandas) == 16 else \
            np.array([planck_integral(a, b, Ts) for a, b in self.bandas])
        U = D = 0.0
        for ib in range(len(self.bandas)):
            ok = self.pesos[ib] > 0
            w = self.pesos[ib][ok]
            frh = self.fraccion(ib, Th)[ok]
            frs = self.fraccion(ib, [Ts])[ok][:, 0]
            u_, d_ = dos_flujos(tau[ib][ok], w, frh * Bh[ib][None, :] / w[:, None],
                                frs * float(np.ravel(Bs[ib])[0]) / w, emisividad)
            U, D = U + u_, D + d_
        return U, D

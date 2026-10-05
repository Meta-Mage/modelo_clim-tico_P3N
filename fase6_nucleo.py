# fase6_nucleo.py -- Fase 6, etapa 6.2: nucleo dinamico PROPIO de M3N (ecuaciones primitivas, seco).
#
# Diseno completo y fuentes: DISENO_FASE6_2.md. Resumen:
#   - Ecuaciones primitivas hidrostaticas en coordenada sigma (Phillips 1957 ✅), forma vectorial invariante.
#   - Vertical: Simmons y Burridge (1981) en sigma pura (✅² IFS Cy25r1 + codigo GFDL), seminiveles de M3N.
#   - Horizontal: rejilla C de la 6.1 (decision de Carlos 05/10), Coriolis de Sadourny con flujos de masa.
#     DEMOSTRACION PROPIA (DISENO_FASE6_2.md §3 bis): la energia total se conserva exactamente en el sistema
#     semidiscreto; lo comprueba test_fase6_2.py.
#   - Filtro polar de la 6.1 (desde 60 grados) sobre todas las tendencias.
#   - Tiempo: leapfrog + RAW (Williams 2009 ✅); semiimplicito (DISENO_FASE6_2.md §4) en una segunda fase.
#
# Convenios: capa k = 0 arriba ... N-1 abajo. Seminiveles sh[0] = 0 (tope, p = 0) ... sh[N] = 1 (suelo).
# Estado: ps (F, C) [Pa]; T (N, F, C) [K]; u (N, F, C) [m/s] en la cara oeste de cada celda;
# v (N, F-1, C) [m/s] en la cara sur de las filas 0..F-2 (positivo hacia el norte).
# W (N+1, F, C) [Pa/s]: flujo de masa vertical HACIA ABAJO por los seminiveles (ps * dsigma/dt), 0 arriba y abajo.

import math
import numpy as np

from fase6_aguas_someras import Rejilla, _filtro_fourier, LAT_FILTRO_DEFECTO, RAW_NU, RAW_ALFA


class NucleoSeco:
    def __init__(self, radio, omega, g, R_gas, cp, sigma_semi, phis=None, filas=36,
                 lat_filtro=LAT_FILTRO_DEFECTO):
        self.Rj = Rj = Rejilla(radio, filas)
        self.omega, self.g, self.R, self.cp = omega, g, R_gas, cp
        self.kappa = R_gas / cp
        sh = np.asarray(sigma_semi, dtype=float)
        assert sh[0] == 0.0 and abs(sh[-1] - 1.0) < 1e-15 and np.all(np.diff(sh) > 0)
        self.sh = sh
        self.N = N = len(sh) - 1
        self.dsig = np.diff(sh)                                               # (N,)
        # ln(sigma_{k+1/2}/sigma_{k-1/2}); en la capa de arriba no se usa (multiplica a 0) -> 0
        self.dln = np.zeros(N)
        self.dln[1:] = np.log(sh[2:] / sh[1:-1])
        # alfa_k de Simmons y Burridge: alfa_0 = ln 2 (tope a p = 0); alfa_k = 1 - sigma_{k-1/2}/dsigma ln(...)
        self.alfa = np.empty(N)
        self.alfa[0] = math.log(2.0)
        self.alfa[1:] = 1.0 - sh[1:-1] / self.dsig[1:] * self.dln[1:]
        self.phis = np.zeros((Rj.filas, Rj.columnas)) if phis is None else phis     # geopotencial del suelo
        self.f_esq = 2 * omega * np.sin(Rj.phi_f[1:-1])[:, None] * np.ones((1, Rj.columnas))
        self.cos_ref = math.cos(math.radians(lat_filtro))
        self.filtrar = lat_filtro < 90.0
        self._d3 = self.dsig[:, None, None]

    # ------------------------------------------------------------------ utilidades
    @staticmethod
    def _oeste(x):            # valor de la celda vecina oeste (para caras oeste)
        return np.roll(x, 1, axis=-1)

    def geopotencial(self, T):
        """Phi en el centro de cada capa (SB81): Phi_k = Phis + sum_{j>k} R T_j dln_j + alfa_k R T_k."""
        R = self.R
        cap = R * T * self.dln[:, None, None]                                  # salto de Phi en cada capa
        debajo = np.cumsum(cap[::-1], axis=0)[::-1] - cap                       # sum_{j>k}
        return self.phis[None] + debajo + self.alfa[:, None, None] * R * T

    def energia_cinetica(self, u, v):
        Rj = self.Rj
        ku = Rj.A_u * u * u
        kv = np.zeros(u.shape[:-2] + (Rj.filas + 1, Rj.columnas)); kv[..., 1:-1, :] = Rj.A_v * v * v
        return (0.5 * (ku + np.roll(ku, -1, axis=-1)) + 0.5 * (kv[..., :-1, :] + kv[..., 1:, :])) / (2 * Rj.area)

    def vorticidad(self, u, v):
        Rj = self.Rj
        circ = (Rj.dy * (v - np.roll(v, 1, axis=-1)) - u[..., :-1, :] * Rj.L_n[:-1] + u[..., 1:, :] * Rj.L_n[1:])
        return circ / Rj.area_dual

    @staticmethod
    def _esquina(x):
        return 0.25 * (x[..., :-1, :] + x[..., 1:, :] + np.roll(x[..., :-1, :], 1, axis=-1)
                       + np.roll(x[..., 1:, :], 1, axis=-1))

    def _divergencia(self, F_u, F_v):
        """Divergencia de flujos por cara (salida neta por celda, mismas unidades que F)."""
        div = np.roll(F_u, -1, axis=-1) - F_u
        div[..., :-1, :] -= F_v
        div[..., 1:, :] += F_v
        return div

    # ------------------------------------------------------------------ tendencias explicitas completas
    def tendencias(self, ps, T, u, v, diag=None):
        Rj, R, cp = self.Rj, self.R, self.cp
        N = self.N
        dsig = self._d3
        # presion y espesores
        ps_u = 0.5 * (ps + self._oeste(ps))
        ps_v = 0.5 * (ps[:-1] + ps[1:])
        dp = dsig * ps[None]                                                    # (N, F, C)
        F_u = u * (dsig * ps_u[None]) * Rj.L_u                                  # Pa m2/s por la cara oeste
        F_v = v * (dsig * ps_v[None]) * Rj.L_v                                  # Pa m2/s por la cara sur
        div = self._divergencia(F_u, F_v)                                       # Pa m2/s
        D = div / Rj.area                                                       # nabla.(dp v) [Pa/s]
        dps = -D.sum(axis=0)
        # flujo de masa vertical hacia abajo por los seminiveles interiores
        W = np.zeros((N + 1,) + ps.shape)
        W[1:-1] = -np.cumsum(D, axis=0)[:-1] - self.sh[1:-1, None, None] * dps[None]
        # ---------- temperatura (forma de flujo para la entalpia) ----------
        T_u = 0.5 * (T + self._oeste(T))
        T_v = 0.5 * (T[:, :-1] + T[:, 1:])
        adv_h = self._divergencia(F_u * T_u, F_v * T_v) / Rj.area             # nabla.(dp v T)
        T_semi = np.zeros((N + 1,) + ps.shape)
        T_semi[1:-1] = 0.5 * (T[:-1] + T[1:])
        flujo_v = W * T_semi
        adv_v = flujo_v[1:] - flujo_v[:-1]
        # conversion de energia, parte 1 (SB81): -R T_k [dln_k sum_{j<k} D_j + alfa_k D_k]
        D_encima = np.cumsum(D, axis=0) - D
        conv1 = -R * T * (self.dln[:, None, None] * D_encima + self.alfa[:, None, None] * D)
        # parte 2: R T v.grad(ln ps) * dp, construida desde las caras (cada cara reparte mitad y mitad)
        g_u = (np.log(ps) - self._oeste(np.log(ps))) / Rj.dx_u                  # d ln ps / dx en la cara
        g_v = (np.log(ps[:-1]) - np.log(ps[1:])) / Rj.dy                        # d ln ps / dy (norte +)
        Y_u = Rj.A_u * (dsig * ps_u[None]) * u * R * T_u * g_u                  # trabajo por cara
        Y_v = Rj.A_v * (dsig * ps_v[None]) * v * R * T_v * g_v
        X = 0.5 * (Y_u + np.roll(Y_u, -1, axis=-1))
        X[:, :-1] += 0.5 * Y_v
        X[:, 1:] += 0.5 * Y_v
        conv2 = X / Rj.area
        d_dpT = -adv_h - adv_v + (conv1 + conv2) / cp                            # d(dp T)/dt
        d_dp = dsig * dps[None]
        dT = (d_dpT - T * d_dp) / dp
        # ---------- momento ----------
        dp_esq = self._esquina(dp)
        q = (self.vorticidad(u, v) + self.f_esq[None]) / dp_esq                  # (N, F-1, C)
        qV = q * 0.5 * (F_v + np.roll(F_v, 1, axis=-1))
        Fu = np.zeros_like(u)
        Fu[:, 1:] += 0.5 * qV
        Fu[:, :-1] += 0.5 * qV
        Fu /= Rj.dx_u
        qU = q * 0.5 * (F_u[:, :-1] + F_u[:, 1:])
        Fv = -0.5 * (qU + np.roll(qU, -1, axis=-1)) / Rj.dy
        Phi = self.geopotencial(T)
        B = self.energia_cinetica(u, v) + Phi
        du = Fu - (B - self._oeste(B)) / Rj.dx_u - R * T_u * g_u
        dv = Fv - (B[:, :-1] - B[:, 1:]) / Rj.dy - R * T_v * g_v
        # advection vertical (SB81, ec. 2.19 IFS) con W promediado a las caras
        W_u = 0.5 * (W + self._oeste(W))
        W_v = 0.5 * (W[:, :-1] + W[:, 1:])
        du -= self._adv_vertical(u, W_u, dsig * ps_u[None])
        dv -= self._adv_vertical(v, W_v, dsig * ps_v[None])
        if diag is not None:
            diag.update(W=W, D=D, Phi=Phi, F_u=F_u, F_v=F_v)
        if self.filtrar:
            dps = self._filtro3(dps, Rj.cos_c)
            dT = self._filtro3(dT, Rj.cos_c)
            du = self._filtro3(du, Rj.cos_c)
            dv = self._filtro3(dv, Rj.cos_v)
        return dps, dT, du, dv

    @staticmethod
    def _adv_vertical(X, W, dp):
        """(1/(2 dp_k)) [W_{k+1/2}(X_{k+1} - X_k) + W_{k-1/2}(X_k - X_{k-1})] (W = 0 arriba y abajo)."""
        dX = X[1:] - X[:-1]                                                     # en los seminiveles interiores
        t = np.zeros_like(X)
        t[:-1] += W[1:-1] * dX
        t[1:] += W[1:-1] * dX
        return t / (2 * dp)

    def _filtro3(self, x, cos_filas):
        """Filtro polar de _filtro_fourier aplicado a todas las capas a la vez (mismos factores S)."""
        clave = len(cos_filas)
        cache = self.__dict__.setdefault("_cache_filtro", {})
        if clave not in cache:
            S = self._S(cos_filas)
            filas = np.where((S != 1.0).any(axis=1))[0]
            cache[clave] = (filas, S[filas])
        filas, S = cache[clave]
        if len(filas) == 0:
            return x
        out = x.copy()
        out[..., filas, :] = np.fft.irfft(np.fft.rfft(x[..., filas, :], axis=-1) * S, n=x.shape[-1], axis=-1)
        return out

    # ------------------------------------------------------------------ diagnosticos
    def integrales(self, ps, T, u, v):
        """Masa total del aire y energia total (cinetica + entalpia + potencial de la superficie) [J]."""
        Rj, g = self.Rj, self.g
        dp = self._d3 * ps[None]
        K = self.energia_cinetica(u, v)
        masa = float((ps * Rj.area).sum() / g)
        energia = float((((K + self.cp * T) * dp).sum(axis=0) + self.phis * ps) .__mul__(Rj.area).sum() / g)
        return {"masa": masa, "energia": energia}

    def derivada_energia(self, ps, T, u, v, tend):
        """dE/dt exacta del sistema semidiscreto a partir de las tendencias (para la prueba de conservacion),
        separada en sus contribuciones para poder comparar la suma con su escala."""
        Rj, g, cp = self.Rj, self.g, self.cp
        dps, dT, du, dv = tend
        dsig = self._d3
        dp = dsig * ps[None]
        d_dp = dsig * dps[None]
        ps_u = 0.5 * (ps + self._oeste(ps)); ps_v = 0.5 * (ps[:-1] + ps[1:])
        dps_u = 0.5 * (dps + self._oeste(dps)); dps_v = 0.5 * (dps[:-1] + dps[1:])
        # E = sum_f A_f dp_f u^2/2 + ... + sum cp T dp + phis ps   (todo / g)
        t_u = (Rj.A_u * (dsig * ps_u[None]) * u * du).sum() + (Rj.A_u * (dsig * dps_u[None]) * u * u / 2).sum()
        t_v = (Rj.A_v * (dsig * ps_v[None]) * v * dv).sum() + (Rj.A_v * (dsig * dps_v[None]) * v * v / 2).sum()
        t_T = (cp * (dT * dp + T * d_dp) * Rj.area).sum()
        t_s = (self.phis * dps * Rj.area).sum()
        partes = np.array([t_u, t_v, t_T, t_s]) / g
        return float(partes.sum()), partes

    # ------------------------------------------------------------------ semiimplicito (DISENO_FASE6_2.md §4)
    def matrices_verticales(self, T_ref):
        """G (hidrostatica), tau (conversion de energia linealizada) y M = G tau + R T_ref 1 dsigma^T."""
        N, R = self.N, self.R
        G = np.zeros((N, N)); tau = np.zeros((N, N))
        for i in range(N):
            G[i, i + 1:] = R * self.dln[i + 1:]
            G[i, i] = R * self.alfa[i]
            tau[i, :i] = self.kappa * T_ref * self.dln[i] * self.dsig[:i] / self.dsig[i]
            tau[i, i] = self.kappa * T_ref * self.alfa[i]
        M = G @ tau + R * T_ref * np.outer(np.ones(N), self.dsig)
        return G, tau, M

    def _S(self, cos_filas):
        """Factores del filtro polar S(fila, k), identicos a _filtro_fourier."""
        Rj = self.Rj
        k = np.arange(Rj.columnas // 2 + 1)
        s = np.abs(np.sin(k * Rj.dlam / 2)); s[0] = 1.0
        S = np.ones((len(cos_filas), len(k)))
        if self.filtrar:
            f = cos_filas < self.cos_ref
            S[f] = np.minimum(1.0, cos_filas[f][:, None] / (self.cos_ref * s[None, :]))
            S[f, 0] = 1.0
        return S

    def preparar_semiimplicito(self, dt, T_ref=300.0, p_ref=1.0e5, beta=0.5):
        """Precalcula los modos verticales y, para cada modo y cada onda zonal, la inversa del operador
        de Helmholtz (I - (2 beta dt)^2 c_m^2 H_k), con H_k = S div S grad en la rejilla C (exacto)."""
        Rj = self.Rj
        F = Rj.filas
        self.si_dt, self.si_Tref, self.si_pref, self.si_beta = dt, T_ref, p_ref, beta
        G, tau, M = self.matrices_verticales(T_ref)
        lam, E = np.linalg.eig(M)
        assert np.all(np.abs(lam.imag) < 1e-9 * np.abs(lam.real).max()) and np.all(lam.real > 0)
        orden = np.argsort(-lam.real)
        self.si_c2 = lam.real[orden]; self.si_E = E.real[:, orden]; self.si_Einv = np.linalg.inv(self.si_E)
        self.si_G, self.si_tau = G, tau
        Sc = self._S(Rj.cos_c)                                                  # (F, K)
        Sv = self._S(Rj.cos_v)                                                  # (F-1, K)
        k = np.arange(Rj.columnas // 2 + 1)
        zonal = -4 * np.sin(k * Rj.dlam / 2) ** 2                               # (e^{ik}-1)(1-e^{-ik})
        dx = Rj.dx_u[:, 0]; area = Rj.area[:, 0]; Lv = Rj.L_v[:, 0]
        nK = len(k)
        H = np.zeros((nK, F, F))
        for i in range(F):
            H[:, i, i] += Sc[i] * Sc[i] * Rj.L_u / area[i] * zonal / dx[i]
            if i > 0:            # cara norte de la fila i = cara sur de la fila i-1 (indice i-1)
                c = Sc[i] * Lv[i - 1] * Sv[i - 1] / (Rj.dy * area[i])
                H[:, i, i - 1] += c; H[:, i, i] -= c
            if i < F - 1:        # cara sur de la fila i (indice i)
                c = Sc[i] * Lv[i] * Sv[i] / (Rj.dy * area[i])
                H[:, i, i + 1] += c; H[:, i, i] -= c
        a2 = (2 * beta * dt) ** 2
        I = np.eye(F)
        self.si_Hinv = np.linalg.inv(I[None, None] - a2 * self.si_c2[:, None, None, None] * H[None])  # (N,K,F,F)
        self.si_H = H

    def _grad(self, P):
        Rj = self.Rj
        return (P - self._oeste(P)) / Rj.dx_u, (P[:, :-1] - P[:, 1:]) / Rj.dy

    def _div_simple(self, u, v):
        Rj = self.Rj
        return self._divergencia(u * Rj.L_u, v * Rj.L_v) / Rj.area

    def _lineal_P(self, T, ps):
        return np.tensordot(self.si_G, T, axes=1) + (self.R * self.si_Tref / self.si_pref) * ps[None]

    def paso_semiimplicito(self, ant, act, forzamiento=None):
        """x^{n+1} = x^{n-1} + 2dt [N(x^n) - L(x^n)] + 2dt [beta L(x^{n+1}) + (1-beta) L(x^{n-1})],
        con L la parte lineal de las ondas de gravedad alrededor de (T_ref, p_ref) en reposo."""
        dt, beta = self.si_dt, self.si_beta
        Rj = self.Rj
        self.ultimo_diag = {} if getattr(self, "guardar_flujos", False) else None
        dps, dT, du, dv = self.tendencias(*act, diag=self.ultimo_diag)
        if forzamiento is not None:
            f = forzamiento(*ant)
            dps, dT, du, dv = dps + f[0], dT + f[1], du + f[2], dv + f[3]
        ps0, T0, u0, v0 = ant
        ps1, T1, u1, v1 = act
        p_ref = self.si_pref
        filt = (lambda x, c: self._filtro3(x, c)) if self.filtrar else (lambda x, c: x)
        filt2 = (lambda x: self._filtro3(x, Rj.cos_c)) if self.filtrar else (lambda x: x)
        # parte lineal en n y n-1
        def L(ps, T, u, v):
            gu, gv = self._grad(self._lineal_P(T, ps))
            delta = self._div_simple(u, v)
            return (-p_ref * np.tensordot(self.dsig, delta, axes=1),
                    -np.tensordot(self.si_tau, delta, axes=1), -gu, -gv)
        L1 = L(*act); L0 = L(*ant)
        L1 = (filt2(L1[0]), filt(L1[1], Rj.cos_c), filt(L1[2], Rj.cos_c), filt(L1[3], Rj.cos_v))
        L0 = (filt2(L0[0]), filt(L0[1], Rj.cos_c), filt(L0[2], Rj.cos_c), filt(L0[3], Rj.cos_v))
        A = [x0 + 2 * dt * (d - l1) + 2 * dt * (1 - beta) * l0
             for x0, d, l1, l0 in zip(ant, (dps, dT, du, dv), L1, L0)]
        A_ps, A_T, A_u, A_v = A
        b = 2 * dt * beta
        # P^{n+1} - b^2 M S div S grad P^{n+1} = Q A_Z - b M S div(A_u)
        M = self.si_G @ self.si_tau + self.R * self.si_Tref * np.outer(np.ones(self.N), self.dsig)
        rhs = self._lineal_P(A_T, A_ps) - b * np.tensordot(M, filt(self._div_simple(A_u, A_v), Rj.cos_c), axes=1)
        r_modo = np.tensordot(self.si_Einv, rhs, axes=1)
        rk = np.fft.rfft(r_modo, axis=-1)                                       # (N, F, K)
        rkt = np.ascontiguousarray(rk.transpose(0, 2, 1))[..., None]                      # (N, K, F, 1)
        Pk = (np.matmul(self.si_Hinv, rkt.real) + 1j * np.matmul(self.si_Hinv, rkt.imag))[..., 0].transpose(0, 2, 1)
        P_modo = np.fft.irfft(Pk, n=Rj.columnas, axis=-1)
        P = np.tensordot(self.si_E, P_modo, axes=1)
        gu, gv = self._grad(P)
        u_n = A_u - b * filt(gu, Rj.cos_c)
        v_n = A_v - b * filt(gv, Rj.cos_v)
        delta = filt(self._div_simple(u_n, v_n), Rj.cos_c)
        T_n = A_T - b * np.tensordot(self.si_tau, delta, axes=1)
        ps_n = A_ps - b * p_ref * np.tensordot(self.dsig, delta, axes=1)
        return ps_n, T_n, u_n, v_n

    def integrar_si(self, ps, T, u, v, pasos, forzamiento=None, cada=None, al_registrar=None, nu=RAW_NU,
                    reanudar=None, guardar=None, cada_guardar=None):
        """Leapfrog semiimplicito + RAW. El primer paso es un medio paso semiimplicito (ant = act).
        reanudar = (n, ant, act) continua exactamente donde se dejo; guardar(n, ant, act) se llama cada
        cada_guardar pasos (punto de control: el resultado es identico al de una corrida sin cortes)."""
        dt = self.si_dt
        if reanudar is not None:
            n0, ant, act = reanudar
        else:
            ant = (ps.copy(), T.copy(), u.copy(), v.copy())
            # arranque: paso con dt/2 desde el mismo estado (equivale a Euler semiimplicito de dt)
            self.si_dt = dt / 2
            self._rehacer_hinv()
            act = self.paso_semiimplicito(ant, ant, forzamiento)
            self.si_dt = dt
            n0 = 1
        self._rehacer_hinv()
        dif = getattr(self, "nu4", None) is not None
        for n in range(n0, pasos):
            nue = self.paso_semiimplicito(ant, act, forzamiento)
            if dif:                                  # hiperdifusion implicita sobre el salto de 2 dt
                T_d, u_d, v_d = self.aplicar_hiperdifusion(nue[1], nue[2], nue[3], 2 * dt)
                nue = (nue[0], T_d, u_d, v_d)
            corr = tuple(0.5 * nu * (a - 2 * b + c) for a, b, c in zip(ant, act, nue))
            ant = tuple(b + RAW_ALFA * c for b, c in zip(act, corr))
            act = tuple(x + (RAW_ALFA - 1) * c for x, c in zip(nue, corr))
            if al_registrar is not None and cada and (n + 1) % cada == 0:
                al_registrar(n + 1, *act)
            if guardar is not None and cada_guardar and (n + 1) % cada_guardar == 0:
                guardar(n + 1, ant, act)
        return act

    def _rehacer_hinv(self):
        F = self.Rj.filas
        a2 = (2 * self.si_beta * self.si_dt) ** 2
        if getattr(self, "_hinv_cache", None) is None:
            self._hinv_cache = {}
        clave = round(self.si_dt, 9)
        if clave not in self._hinv_cache:
            self._hinv_cache[clave] = np.linalg.inv(np.eye(F)[None, None] - a2 * self.si_c2[:, None, None, None] * self.si_H[None])
        self.si_Hinv = self._hinv_cache[clave]

    # ------------------------------------------------------------------ hiperdifusion (DISENO_FASE6_2.md §7)
    def _laplaciano_escalar(self, x):
        """div grad en la rejilla C (sin filtro). x: (..., F, C), puede ser complejo."""
        Rj = self.Rj
        gu = (x - np.roll(x, 1, axis=-1)) / Rj.dx_u
        gv = (x[..., :-1, :] - x[..., 1:, :]) / Rj.dy
        return self._divergencia(gu * Rj.L_u, gv * Rj.L_v) / Rj.area

    def _laplaciano_vector(self, u, v):
        """Laplaciano vectorial grad(div) - rot(rot) en la rejilla C. En los polos la vorticidad es la
        circulacion de la fila de u mas cercana dividida por el area del casquete."""
        Rj = self.Rj
        d = self._div_simple(u, v)
        z_int = self.vorticidad(u, v)                                          # esquinas interiores (F-1)
        cap = 2 * math.pi * Rj.a ** 2 * (1 - math.sin(Rj.phi_c[0]))
        circ_n = (u[..., 0, :] * Rj.L_n[0, 0]).sum(axis=-1, keepdims=True)
        circ_s = -(u[..., -1, :] * Rj.L_n[-1, 0]).sum(axis=-1, keepdims=True)
        z_n = np.broadcast_to(circ_n / cap, u[..., 0, :].shape)
        z_s = np.broadcast_to(circ_s / cap, u[..., 0, :].shape)
        z = np.concatenate([z_n[..., None, :], z_int, z_s[..., None, :]], axis=-2)   # (F+1) caras
        Lu = (d - np.roll(d, 1, axis=-1)) / Rj.dx_u - (z[..., :-1, :] - z[..., 1:, :]) / Rj.dy
        Lv = (d[..., :-1, :] - d[..., 1:, :]) / Rj.dy + (np.roll(z_int, -1, axis=-1) - z_int) / Rj.L_v
        return Lu, Lv

    def preparar_hiperdifusion(self, tau_dias=0.5, dt_efectivo=None):
        """Hiperdifusion nabla^4 IMPLICITA sobre T, u y v: (I + dt nu4 L^2) x_nuevo = x, resuelta exactamente
        para cada onda zonal (los operadores no dependen de la longitud). nu4 se fija para que la onda
        zonal mas corta en el ecuador (autovalor discreto -4/dx^2) se amortigue en tau_dias."""
        Rj = self.Rj
        F, C = Rj.filas, Rj.columnas
        dx_ec = Rj.a * Rj.dlam                                                   # anchura en el ecuador
        self.nu4 = 1.0 / (tau_dias * 86400.0 * (4.0 / dx_ec ** 2) ** 2)
        self.dif_dt = dt_efectivo
        nK = C // 2 + 1
        lam = Rj.lam_c                                                           # mismas fases para todas las filas
        Ls = np.zeros((nK, F, F), complex)
        Lvec = np.zeros((nK, 2 * F - 1, 2 * F - 1), complex)
        for k in range(nK):
            base_h = np.exp(1j * k * Rj.lam_c)
            base_u = np.exp(1j * k * Rj.lam_u)
            for r in range(F):
                x = np.zeros((F, C), complex); x[r] = base_h
                Ls[k, :, r] = (self._laplaciano_escalar(x) / base_h[None]).mean(axis=-1)
                uu = np.zeros((F, C), complex); uu[r] = base_u
                Lu, Lv = self._laplaciano_vector(uu, np.zeros((F - 1, C), complex))
                Lvec[k, :F, r] = (Lu / base_u[None]).mean(axis=-1)
                Lvec[k, F:, r] = (Lv / base_h[None]).mean(axis=-1)
            for r in range(F - 1):
                vv = np.zeros((F - 1, C), complex); vv[r] = base_h
                Lu, Lv = self._laplaciano_vector(np.zeros((F, C), complex), vv)
                Lvec[k, :F, F + r] = (Lu / base_u[None]).mean(axis=-1)
                Lvec[k, F:, F + r] = (Lv / base_h[None]).mean(axis=-1)
        self._Ls, self._Lvec = Ls, Lvec
        self._dif_inv = {}

    def _inversas_difusion(self, dt):
        clave = round(dt, 6)
        if clave not in self._dif_inv:
            a = dt * self.nu4
            Is = np.eye(self._Ls.shape[1]); Iv = np.eye(self._Lvec.shape[1])
            inv_s = np.linalg.inv(Is[None] + a * self._Ls @ self._Ls)
            assert np.abs(inv_s.imag).max() < 1e-12 * np.abs(inv_s.real).max()   # el laplaciano escalar es real
            self._dif_inv[clave] = (np.ascontiguousarray(inv_s.real),
                                    np.linalg.inv(Iv[None] + a * self._Lvec @ self._Lvec))
        return self._dif_inv[clave]

    def aplicar_hiperdifusion(self, T, u, v, dt):
        """Paso implicito de la hiperdifusion de duracion dt (separado del resto de la dinamica)."""
        Rj = self.Rj
        F, C = Rj.filas, Rj.columnas
        Is, Iv = self._inversas_difusion(dt)
        Tk = np.ascontiguousarray(np.fft.rfft(T, axis=-1).transpose(2, 1, 0))               # (K, F, N)
        T_n = np.fft.irfft((np.matmul(Is, Tk.real) + 1j * np.matmul(Is, Tk.imag)).transpose(2, 1, 0), n=C, axis=-1)
        # coeficientes en las bases de las matrices: e^{ik lam_u} para u y e^{ik lam_c} para v
        K = np.arange(C // 2 + 1)
        fu = np.exp(1j * K * Rj.lam_u[0]); fh = np.exp(1j * K * Rj.lam_c[0])
        w = np.concatenate([np.fft.rfft(u, axis=-1) / fu, np.fft.rfft(v, axis=-1) / fh], axis=-2)
        w_n = np.matmul(Iv, np.ascontiguousarray(w.transpose(2, 1, 0))).transpose(2, 1, 0)
        u_n = np.fft.irfft(w_n[:, :F] * fu, n=C, axis=-1)
        v_n = np.fft.irfft(w_n[:, F:] * fh, n=C, axis=-1)
        return T_n, u_n, v_n

    # ------------------------------------------------------------------ avance (explicito)
    def integrar(self, ps, T, u, v, dt, pasos, forzamiento=None, cada=None, al_registrar=None):
        """Leapfrog + RAW. forzamiento(ps, T, u, v) -> tendencias adicionales (se evaluan en n-1)."""
        ant = (ps.copy(), T.copy(), u.copy(), v.copy())
        d = self.tendencias(*ant)
        if forzamiento is not None:
            d = tuple(a + b for a, b in zip(d, forzamiento(*ant)))
        act = tuple(x + dt * dx for x, dx in zip(ant, d))
        for n in range(1, pasos):
            d = self.tendencias(*act)
            if forzamiento is not None:
                d = tuple(a + b for a, b in zip(d, forzamiento(*ant)))
            nue = tuple(x + 2 * dt * dx for x, dx in zip(ant, d))
            corr = tuple(0.5 * RAW_NU * (a - 2 * b + c) for a, b, c in zip(ant, act, nue))
            ant = tuple(b + RAW_ALFA * c for b, c in zip(act, corr))
            act = tuple(x + (RAW_ALFA - 1) * c for x, c in zip(nue, corr))
            if al_registrar is not None and cada and (n + 1) % cada == 0:
                al_registrar(n + 1, *act)
        return act

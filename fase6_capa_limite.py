# fase6_capa_limite.py -- Fase 6.3 (paso 4): capa limite y rozamiento con la superficie, en columna.
#
# Diseno: DISENO_FASE6_3.md §6.2 (decisiones 1.3, 1.4 y 1.8, aprobadas por Carlos el 06/10/2026) y §6.5.
# Modulo AUTONOMO: no esta conectado a M3N (eso es el paso 5 de la 6.3, con el interruptor I16).
#
# Que hace, en cada columna y en un paso dt (implicito, Euler hacia atras, incondicionalmente estable):
#   1. Coeficientes de intercambio con la superficie a la ALTURA REAL de la capa mas baja (z_a), no a 10 m:
#        C_N = [kappa / ln(z_a/z0m)]^2                       (neutro, momento)
#        C_HN = kappa^2 / [ln(z_a/z0m) ln(z_a/z0h)]          (neutro, calor)
#      y estabilidad de Louis, Tiedtke y Geleyn (1982), Tabla 1, b = c = d = 5 ✅ (verificada en el codigo de
#      LMDZ, cdrag_mod.F90):
#        momento: f_m = 1 - 2bRi/(1 + 3bc C_N sqrt(z_a/z0m |Ri|))   (Ri < 0);  1/(1 + 2bRi/sqrt(1 + dRi))  (Ri >= 0)
#        calor:   f_h = 1 - 3bRi/(1 + 3bc C_N sqrt(z_a/z0m |Ri|))   (Ri < 0);  1/(1 + 3bRi sqrt(1 + dRi))  (Ri >= 0)
#      (la de calor es la misma que factor_estabilidad_louis de fase2b_atmosfera.py con z_a = Z_REF).
#   2. Arrastre en superficie tau = rho C_N f_m |v_a| v_a, implicito en v_a (|v_a| del paso anterior). Sin
#      rafaga (Frierson, Held y Zurita-Gotor 2006 ✅: "zero gustiness"); solo un minimo NUMERICO del viento
#      (V_MIN_NUMERICO) para no dividir por cero en Ri.
#   3. Altura de la capa limite h: donde el numero de Richardson global, medido desde la capa mas baja,
#      supera Ri_c = 1 (interpolacion lineal), como diffusivity.F90 de GFDL/Isca con do_simple ✅.
#   4. Difusividad K(z) en las interfaces (diffusivity.F90 + monin_obukhov.F90 de GFDL/Isca ✅):
#        z < f h  (f = 0,1):     K = kappa u* z / phi(z/L)                 (Monin-Obukhov)
#        f h <= z < h:           K = K(f h) (z / f h) [1 - (z - f h)/(h - f h)]^2
#        z >= h:                 K = 0
#      phi_m = (1 - 16 zeta)^(-1/4), phi_h = (1 - 16 zeta)^(-1/2) si zeta < 0;
#      phi_m = phi_h = 1 + zeta (5 + b_stab zeta)/(1 + zeta) si zeta >= 0, b_stab = 1/rich_crit = 0,5.
#      L de Monin-Obukhov con u* = sqrt(C_N f_m) |v_a| y el flujo de calor sensible de la superficie (entrada).
#   5. Difusion vertical implicita de u, v (con el arrastre), del calor (flujo con el gradiente de theta, en
#      forma conservativa de energia: ver paso_capa_limite) y del vapor q. Los flujos de calor y de vapor con la SUPERFICIE no se aplican aqui
#      (los calcula la fisica de M3N, que ya cierra su balance): para s y q el fondo es de flujo cero.
#   6. Calor por rozamiento (como CAM, difcor.F90 ✅): la energia cinetica que pierde cada capa vuelve a esa
#      misma capa como calor, dT_k = -dKE_k / cp. Asi sum_k dp_k (cp T_k + KE_k) se conserva EXACTAMENTE.
#
# Convenios: columnas aplanadas. Niveles k = 0 (arriba) ... N-1 (abajo), como el nucleo y fase30_multicapa.
#   u, v, T, q: (N, M); ph: (N+1, M) presion en las semicapas (Pa, ph[0] = tope); pm: (N, M) presion de cada
#   capa; Ts, z0m, z0h, flujo_calor (W/m2, positivo hacia arriba): (M,).

import numpy as np

KARMAN = 0.4
LOUIS_B = LOUIS_C = LOUIS_D = 5.0
RI_CRITICO_CL = 1.0          # Frierson et al. 2006 / GFDL rich_crit_pbl
FRACCION_INTERIOR = 0.1      # f (Frierson et al. 2006 / GFDL frac_inner)
B_ESTABLE = 1.0 / 2.0        # GFDL monin_obukhov: b_stab = 1/rich_crit, rich_crit = 2
V_MIN_NUMERICO = 0.01        # m/s: solo para no dividir por cero (NO es una rafaga fisica)


def louis_momento(ri, c_n, z_sobre_z0):
    """Factor de estabilidad para el MOMENTO, LTG82 Tabla 1."""
    ri_abs = np.abs(ri)
    inestable = 1 - 2 * LOUIS_B * ri / (1 + 3 * LOUIS_B * LOUIS_C * c_n * np.sqrt(z_sobre_z0 * ri_abs))
    estable = 1 / (1 + 2 * LOUIS_B * ri / np.sqrt(1 + LOUIS_D * ri_abs))
    return np.where(ri < 0, inestable, estable)


def louis_calor(ri, c_n, z_sobre_z0):
    """Factor de estabilidad para el CALOR, LTG82 Tabla 1 (= factor_estabilidad_louis con z_a = Z_REF)."""
    ri_abs = np.abs(ri)
    inestable = 1 - 3 * LOUIS_B * ri / (1 + 3 * LOUIS_B * LOUIS_C * c_n * np.sqrt(z_sobre_z0 * ri_abs))
    estable = 1 / (1 + 3 * LOUIS_B * ri * np.sqrt(1 + LOUIS_D * ri_abs))
    return np.where(ri < 0, inestable, estable)


def phi_monin_obukhov(zeta):
    """phi_m, phi_h de GFDL monin_obukhov.F90 (stable_option = 1)."""
    z_in = np.minimum(zeta, 0.0)
    z_es = np.maximum(zeta, 0.0)
    estable = 1.0 + z_es * (5.0 + B_ESTABLE * z_es) / (1.0 + z_es)
    phi_m = np.where(zeta < 0, (1 - 16.0 * z_in) ** -0.25, estable)
    phi_h = np.where(zeta < 0, (1 - 16.0 * z_in) ** -0.5, estable)
    return phi_m, phi_h


def alturas(T, ph, pm, R_gas, g):
    """Altura sobre el suelo de cada capa (z_capa, (N, M)) y de cada interfaz interior (z_int, (N-1, M);
    z_int[i] separa las capas i y i+1). Hidrostatica con la T de cada capa."""
    N = T.shape[0]
    z_semi = np.zeros_like(ph)
    for k in range(N - 1, 0, -1):                                   # ph[0] puede ser 0: no se usa
        z_semi[k] = z_semi[k + 1] + R_gas * T[k] / g * np.log(ph[k + 1] / ph[k])
    z_capa = z_semi[1:] + R_gas * T / g * np.log(ph[1:] / pm)
    return z_capa, z_semi[1:N]


def altura_capa_limite(theta, u, v, z_capa, g):
    """h donde el Ri global (desde la capa mas baja, como GFDL: z g (theta_k - theta_b)/theta_b / |v_k|^2;
    GFDL usa s/cp, equivalente en el continuo) supera RI_CRITICO_CL, con interpolacion lineal; si no se
    supera en ninguna capa, la altura de la capa mas alta."""
    N, M = theta.shape
    base = theta[-1]
    h = np.full(M, np.nan)
    ri_ant = np.zeros(M)
    z_ant = z_capa[-1].copy()
    for k in range(N - 2, -1, -1):
        v2 = u[k] ** 2 + v[k] ** 2 + V_MIN_NUMERICO ** 2
        ri_k = z_capa[k] * g * (theta[k] - base) / base / v2
        cruza = np.isnan(h) & (ri_k > RI_CRITICO_CL)
        den = np.where(cruza, ri_k - ri_ant, 1.0)          # si cruza, ri_ant <= Ri_c < ri_k: den > 0
        h = np.where(cruza, z_capa[k] + (z_ant - z_capa[k]) * (ri_k - RI_CRITICO_CL) / den, h)
        ri_ant, z_ant = ri_k, z_capa[k]
    return np.where(np.isnan(h), z_capa[0], h)


def difusividades(z_int, h, u_estrella, flujo_cinematico, T_ref, g):
    """K_m, K_h en las interfaces interiores (N-1, M). flujo_cinematico = w'theta' en superficie (K m/s)."""
    # longitud de Monin-Obukhov: zeta = z / L = -kappa g z w'theta' / (T u*^3)
    u3 = np.maximum(u_estrella, 1e-6) ** 3
    def k_mo(z):
        zeta = -KARMAN * g * z * flujo_cinematico / (T_ref * u3)
        pm_, ph_ = phi_monin_obukhov(zeta)
        return KARMAN * u_estrella * z / pm_, KARMAN * u_estrella * z / ph_
    h_in = FRACCION_INTERIOR * h
    km_ref, kh_ref = k_mo(h_in)
    km_sup, kh_sup = k_mo(z_int)
    factor = (z_int / h_in) * (1.0 - (z_int - h_in) / np.maximum(h - h_in, 1e-9)) ** 2
    exterior = (z_int >= h_in) & (z_int < h)
    K_m = np.where(z_int < h_in, km_sup, np.where(exterior, km_ref * factor, 0.0))
    K_h = np.where(z_int < h_in, kh_sup, np.where(exterior, kh_ref * factor, 0.0))
    return K_m, K_h


def _tridiagonal(a, b, c, d):
    """Resuelve a_k x_{k-1} + b_k x_k + c_k x_{k+1} = d_k por columnas (algoritmo de Thomas). (N, M)."""
    N = b.shape[0]
    cp_ = np.empty_like(b); dp_ = np.empty_like(d)
    cp_[0] = c[0] / b[0]; dp_[0] = d[0] / b[0]
    for k in range(1, N):
        den = b[k] - a[k] * cp_[k - 1]
        cp_[k] = c[k] / den if k < N - 1 else 0.0
        dp_[k] = (d[k] - a[k] * dp_[k - 1]) / den
    x = np.empty_like(d)
    x[-1] = dp_[-1]
    for k in range(N - 2, -1, -1):
        x[k] = dp_[k] - cp_[k] * x[k + 1]
    return x


def _difundir(X, masa, D, dt, arrastre=None):
    """masa_k (X_k^new - X_k)/dt = D_{k+1/2}(X_{k+1} - X_k) - D_{k-1/2}(X_k - X_{k-1}) [- arrastre X_{N-1}^new].
    D: (N-1, M) conductancias de las interfaces. Implicito (Euler hacia atras) y conservativo: sum masa X
    solo cambia por el arrastre. Con D >= 0 la matriz es una M-matriz: no crea maximos ni minimos (q >= 0)."""
    N = X.shape[0]
    a = np.zeros_like(X); c = np.zeros_like(X)
    a[1:] = -dt * D / masa[1:]
    c[:-1] = -dt * D / masa[:-1]
    b = 1.0 - a - c
    if arrastre is not None:
        b[-1] = b[-1] + dt * arrastre / masa[-1]
    return _tridiagonal(a, b, c, X)


def paso_capa_limite(u, v, T, ph, pm, Ts, z0m, z0h, flujo_calor, dt, g, R_gas, cp, q=None):
    """Un paso implicito de capa limite. Devuelve (u, v, T, q, diag)."""
    N, M = T.shape
    masa = (ph[1:] - ph[:-1]) / g                                     # kg/m2
    z_capa, z_int = alturas(T, ph, pm, R_gas, g)
    z_a = z_capa[-1]
    # coeficientes de superficie a la altura real de la capa mas baja
    lm = np.log(z_a / z0m); lh = np.log(z_a / z0h)
    c_n = (KARMAN / lm) ** 2
    c_hn = KARMAN ** 2 / (lm * lh)
    theta_s = Ts                                                      # superficie: p = ps
    theta_a = T[-1] * (ph[-1] / pm[-1]) ** (R_gas / cp)              # T de la capa mas baja llevada a ps
    v_a = np.sqrt(u[-1] ** 2 + v[-1] ** 2 + V_MIN_NUMERICO ** 2)
    ri = g * z_a * (theta_a - theta_s) / (0.5 * (theta_a + theta_s) * v_a ** 2)
    c_m = c_n * louis_momento(ri, c_n, z_a / z0m)
    c_h = c_hn * louis_calor(ri, c_n, z_a / z0m)
    rho_s = ph[-1] / (R_gas * theta_a)
    u_est = np.sqrt(c_m) * v_a
    arrastre = rho_s * c_m * v_a                                      # kg/m2/s: tau = arrastre * v_a(nuevo)
    # capa limite y difusividades (h con la temperatura potencial: neutro = theta constante)
    kappa_r = R_gas / cp
    exner = (pm / ph[-1]) ** kappa_r                                  # Pi_k, referido a ps (T = Pi theta)
    exner_int = (ph[1:N] / ph[-1]) ** kappa_r                         # Pi en las interfaces interiores
    theta = T / exner
    h = altura_capa_limite(theta, u, v, z_capa, g)
    K_m, K_h = difusividades(z_int, h, u_est, flujo_calor / (rho_s * cp), theta_a, g)
    T_int = 0.5 * (T[:-1] + T[1:])
    rho_int = ph[1:N] / (R_gas * T_int)
    dz = z_capa[:-1] - z_capa[1:]
    D_m = rho_int * K_m / dz
    D_h = rho_int * K_h / dz
    # momento (con arrastre)
    u1 = _difundir(u, masa, D_m, dt, arrastre)
    v1 = _difundir(v, masa, D_m, dt, arrastre)
    # calor por rozamiento, POSITIVO y EXACTO (demostracion propia, DISENO_FASE6_3.md §6.5). Para Euler
    # implicito, sumando por partes:
    #   sum_k m_k (KE1_k - KE0_k) = -dt sum_i D_i |v1_{i} - v1_{i+1}|^2 - dt arrastre |v1_{N-1}|^2
    #                               - 1/2 sum_k m_k |v1_k - v0_k|^2
    # Cada termino es >= 0 y se deposita donde ocurre: la disipacion de cada interfaz, a medias entre sus dos
    # capas; la del suelo, en la capa mas baja; la del paso implicito, en su capa. Su suma es exactamente la
    # energia cinetica perdida, asi que sum_k m_k (cp T_k + KE_k) se conserva a redondeo.
    dis_int = dt * D_m * ((u1[:-1] - u1[1:]) ** 2 + (v1[:-1] - v1[1:]) ** 2)       # J/m2, (N-1, M)
    calor = 0.5 * masa * ((u1 - u) ** 2 + (v1 - v) ** 2)                            # J/m2, (N, M)
    calor[:-1] += 0.5 * dis_int
    calor[1:] += 0.5 * dis_int
    calor[-1] += dt * arrastre * (u1[-1] ** 2 + v1[-1] ** 2)
    # calor: flujo de la teoria K con el gradiente de theta, H = -rho cp K Pi dtheta/dz, en forma conservativa
    # de energia: masa_k cp Pi_k dtheta_k/dt = D_{k+1/2} cp Pi_{k+1/2} (theta_{k+1} - theta_k) - (...)_{k-1/2}.
    # Asi sum_k masa_k cp T_k se conserva exactamente y una columna neutra (theta constante) no tiene flujo
    # (difundir s = cp T + g z no lo cumple: con la hidrostatica discreta s no es constante en una columna
    # isentropica; hasta 77 J/kg entre las 5 capas bajas de M3N, medido el 06/10).
    theta1 = _difundir(theta, masa * exner, D_h * exner_int, dt)
    T1 = exner * theta1 + calor / (masa * cp)
    q1 = None if q is None else _difundir(q, masa, D_h, dt)
    diag = {"h": h, "u_estrella": u_est, "c_m": c_m, "c_h": c_h, "z_a": z_a, "ri": ri,
            "tau_x": arrastre * u1[-1], "tau_y": arrastre * v1[-1], "K_m": K_m, "K_h": K_h,
            "calor_rozamiento": calor.sum(0) / dt}                    # W/m2 de la columna
    return u1, v1, T1, q1, diag

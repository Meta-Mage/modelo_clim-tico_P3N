# test_fase6_capa_limite.py -- pruebas del modulo de capa limite de la Fase 6.3 (paso 4). Ejecutar con pytest.
import math
import numpy as np

import fase30_multicapa as Mm
import fase2b_atmosfera as F
import fase6_capa_limite as CL

G, R_GAS, CP = 9.80665, 287.05, 1004.0


def _columnas_m3n(M=40, semilla=7, estable=None):
    """Columnas con los 20 niveles de M3N, viento y temperatura aleatorios razonables."""
    rng = np.random.default_rng(semilla)
    sh = np.asarray(Mm.sigma_seminiveles())
    ps = 1e5 * (1 + 0.02 * rng.standard_normal(M))
    ph = sh[:, None] * ps[None]
    pm = np.empty((20, M))
    pm[0] = 0.5 * ph[1]
    pm[1:] = np.exp(0.5 * (np.log(ph[1:-1]) + np.log(ph[2:])))
    T = 290.0 * (pm / ps) ** 0.19 + 2 * rng.standard_normal((20, M))
    T = np.maximum(T, 200.0)
    u = 10 * rng.standard_normal((20, M)); v = 10 * rng.standard_normal((20, M))
    if estable is None:
        Ts = T[-1] + 4 * rng.standard_normal(M)
    else:
        Ts = T[-1] + (-5.0 if estable else 5.0)
    z0m = np.where(rng.random(M) < 0.5, 0.01, 2e-4); z0h = z0m / 10
    flujo = 50 * rng.standard_normal(M)
    q = np.abs(0.01 * rng.standard_normal((20, M)))
    return u, v, T, ph, pm, Ts, z0m, z0h, flujo, q


def test_louis_calor_igual_que_el_de_m3n():
    ri = np.linspace(-5, 5, 2001)
    for z0m in (0.01, 2e-4):
        c_n = (F.KARMAN / np.log(F.Z_REF / z0m)) ** 2
        a = CL.louis_calor(ri, c_n, F.Z_REF / z0m)
        b = F.factor_estabilidad_louis(ri, c_n, z0m)
        assert np.array_equal(a, b)


def test_louis_neutro_monotono_y_momento_menos_frenado():
    ri = np.linspace(-3, 3, 601)
    c_n = (0.4 / math.log(184 / 0.01)) ** 2
    fm = CL.louis_momento(ri, c_n, 184 / 0.01); fh = CL.louis_calor(ri, c_n, 184 / 0.01)
    i0 = 300
    assert fm[i0] == 1.0 and fh[i0] == 1.0
    assert np.all(np.diff(fm) < 0) and np.all(np.diff(fh) < 0)
    assert np.all(fm[i0 + 1:] > fh[i0 + 1:])          # estable: el momento se frena menos que el calor


def test_phi_monin_obukhov():
    pm_, ph_ = CL.phi_monin_obukhov(np.array([-1.0, 0.0, 1e-6, 1.0]))
    assert pm_[0] == 17 ** -0.25 and ph_[0] == 17 ** -0.5
    assert pm_[1] == 1.0 and ph_[1] == 1.0
    assert abs((pm_[2] - 1) / 1e-6 - 5.0) < 1e-4                      # 1 + 5 zeta para zeta pequeno
    assert abs(pm_[3] - (1 + 1 * (5 + 0.5) / 2)) < 1e-15


def test_energia_exacta_calor_positivo_momento_y_agua():
    for estable in (None, True, False):
        u, v, T, ph, pm, Ts, z0m, z0h, flujo, q = _columnas_m3n(estable=estable)
        dt = 992.0
        masa = (ph[1:] - ph[:-1]) / G
        u1, v1, T1, q1, d = CL.paso_capa_limite(u, v, T, ph, pm, Ts, z0m, z0h, flujo, dt, G, R_GAS, CP, q=q)
        e0 = (masa * (CP * T + 0.5 * (u ** 2 + v ** 2))).sum(0)
        e1 = (masa * (CP * T1 + 0.5 * (u1 ** 2 + v1 ** 2))).sum(0)
        assert np.max(np.abs(e1 - e0) / e0) < 1e-14
        assert np.all(d["calor_rozamiento"] >= 0)
        # momento: solo lo quita el rozamiento con el suelo
        assert np.allclose((masa * (u1 - u)).sum(0), -dt * d["tau_x"], rtol=1e-10, atol=1e-9)
        assert np.allclose((masa * (v1 - v)).sum(0), -dt * d["tau_y"], rtol=1e-10, atol=1e-9)
        # vapor: conservado y positivo
        assert np.max(np.abs((masa * (q1 - q)).sum(0)) / (masa * q).sum(0)) < 1e-13
        assert q1.min() >= 0


def test_calor_por_rozamiento_positivo_capa_a_capa():
    u, v, T, ph, pm, Ts, z0m, z0h, flujo, q = _columnas_m3n(semilla=11)
    masa = (ph[1:] - ph[:-1]) / G
    u1, v1, T1, _, d = CL.paso_capa_limite(u, v, T, ph, pm, Ts, z0m, z0h, flujo, 992.0, G, R_GAS, CP)
    # la misma difusion del calor sin el rozamiento: la diferencia (el calor por rozamiento) es >= 0 en cada capa
    z = CL.alturas(T, ph, pm, R_GAS, G)[0]
    ex = (pm / ph[-1]) ** (R_GAS / CP); ex_i = (ph[1:20] / ph[-1]) ** (R_GAS / CP)
    D_h = (ph[1:20] / (R_GAS * 0.5 * (T[:-1] + T[1:]))) * d["K_h"] / (z[:-1] - z[1:])
    T_solo_dif = ex * CL._difundir(T / ex, masa * ex, D_h * ex_i, 992.0)
    assert np.all(T1 - T_solo_dif >= -1e-12)
    assert np.allclose(((T1 - T_solo_dif) * masa * CP).sum(0) / 992.0, d["calor_rozamiento"], rtol=1e-9)


def test_reposo_sin_cambios():
    u, v, T, ph, pm, Ts, z0m, z0h, flujo, q = _columnas_m3n()
    # aire en reposo, columna isentropica y sin flujo: no cambia nada
    th = 300.0
    T = th * (pm / ph[-1]) ** (R_GAS / CP)
    cero = np.zeros_like(T)
    u1, v1, T1, _, d = CL.paso_capa_limite(cero, cero, T, ph, pm, np.full(T.shape[1], th), z0m, z0h,
                                           np.zeros(T.shape[1]), 992.0, G, R_GAS, CP)
    assert np.all(u1 == 0) and np.all(v1 == 0)
    assert np.max(np.abs(T1 - T)) < 1e-10          # neutra: sin flujo de calor (exacto salvo redondeo)


def _columna_fina(n=80, z_tope=2000.0, z_a=2.0, T0=290.0):
    """Columna neutra fina (alturas en progresion geometrica) para la prueba del perfil logaritmico."""
    zs = np.geomspace(z_a, z_tope, n)[::-1]                              # centros, de arriba abajo
    zi = np.sqrt(zs[:-1] * zs[1:])                                       # interfaces interiores
    zsemi = np.concatenate(([2 * zs[0] - zi[0]], zi, [0.0]))            # semicapas (tope arbitrario)
    th = T0
    Tz = lambda z: th - G / CP * z                                       # neutra: theta constante
    p = lambda z: 1e5 * (Tz(z) / th) ** (CP / R_GAS)
    ph = p(zsemi)[:, None]; pm = p(zs)[:, None]; T = Tz(zs)[:, None]
    return zs, ph, pm, T


def test_perfil_logaritmico_en_la_capa_superficial_neutra():
    zs, ph, pm, T = _columna_fina()
    n = len(zs)
    u = np.zeros((n, 1)); v = np.zeros((n, 1))
    z0 = np.array([0.01]); Ts = np.array([290.0])
    Fx = 1e-3                                                            # m/s2, fuerza uniforme (como un gradiente de presion)
    dt = 600.0
    for _ in range(6000):
        u = u + Fx * dt
        u, v, _, _, d = CL.paso_capa_limite(u, v, T, ph, pm, Ts, z0, z0 / 10, np.zeros(1), dt, G, R_GAS, CP)
    u_est = float(d["u_estrella"][0]); h = float(d["h"][0])
    assert h > 1900.0                                                    # neutra: el Ri no llega a 1
    # Solucion estacionaria: du/dz = tau(z) / (rho(z) K(z)), con K = kappa u* z (neutra, capa superficial) y
    # tau(z)/tau_s = (p(z) - p_tope)/(p_s - p_tope) (fuerza uniforme por unidad de masa, hidrostatica). Se usa
    # la tension que el modulo aplica de verdad (tau_x): el arrastre usa |v| de ANTES del paso y la prueba
    # suma la fuerza antes de cada paso, asi que tau_s/rho_s no es exactamente u*^2 (un 3 %, medido).
    th = 290.0
    p = lambda z: 1e5 * ((th - G / CP * z) / th) ** (CP / R_GAS)
    rho = lambda z: p(z) / (R_GAS * (th - G / CP * z))
    rho_s = ph[-1, 0] / (R_GAS * th)
    p_top = ph[0, 0]
    tau_s = float(d["tau_x"][0])
    zz = np.geomspace(4.0, 80.0, 8001)
    tasa = (tau_s * (p(zz) - p_top) / (ph[-1, 0] - p_top)) / (rho(zz) * CL.KARMAN * u_est * zz)
    z_c = CL.alturas(T, ph, pm, R_GAS, G)[0][:, 0]
    sel = np.where((z_c > 5.0) & (z_c < 60.0))[0]
    z_alto, z_bajo = z_c[sel[0]], z_c[sel[-1]]
    du_num = float(u[sel[0], 0] - u[sel[-1], 0])
    m = (zz >= z_bajo) & (zz <= z_alto)
    du_ex = float(np.trapezoid(tasa[m], zz[m]))
    assert abs(du_num / du_ex - 1) < 0.01, (du_num, du_ex)

# test_fase6_3.py -- pruebas del transporte de trazadores (Fase 6.3, paso 3). Ejecutar con pytest.
import math
import numpy as np

from fase6_aguas_someras import Rejilla, campana, caso1_exacta, A_TIERRA, DIA, errores_normalizados
from fase6_trazadores import transportar


def _flujos_caso1(R, alfa):
    """Flujos de volumen por cara exactamente no divergentes, a partir de la funcion de corriente del caso 1
    de Williamson en las esquinas."""
    a = R.a; u0 = 2 * math.pi * a / (12 * DIA)
    lam = np.concatenate([R.lam_u, [R.lam_u[0] + 2 * math.pi]])
    L, P = np.meshgrid(lam, R.phi_f)
    psi = -a * u0 * (np.sin(P) * math.cos(alfa) - np.cos(L) * np.cos(P) * math.sin(alfa))
    return psi[1:, :-1] - psi[:-1, :-1], psi[1:-1, 1:] - psi[1:-1, :-1]


def _corre(alfa, filas=36, dias=12, dt=3600.0):
    R = Rejilla(A_TIERRA, filas)
    Fx, Fy = _flujos_caso1(R, alfa)
    m = R.area.copy(); q = campana(R, -math.pi / 2, 0.0)
    M0 = (m * q).sum(); qmin = 0.0
    for k in range(int(round(dias * DIA / dt))):
        q, m = transportar(q, m, Fx * dt, Fy * dt, orden_xy=(k % 2 == 0))
        qmin = min(qmin, q.min())
    return R, q, m, M0, qmin


def test_campana_por_encima_de_los_polos_positiva_y_conservativa():
    R, q, m, M0, qmin = _corre(math.pi / 2)
    assert abs((m * q).sum() / M0 - 1) < 1e-13
    assert qmin > -1e-100                                    # solo redondeo de denormales
    e = errores_normalizados(q, caso1_exacta(12 * DIA, math.pi / 2), R.area)
    assert e[1] < 0.6                                        # 0,53 medido (el esquema centrado da ~1)


def test_campo_constante_con_viento_divergente_sigue_constante():
    R = Rejilla(A_TIERRA, 36)
    rng = np.random.default_rng(0)
    m = R.area * (1 + 0.1 * rng.random(R.area.shape))
    Mx = 0.3 * m * rng.standard_normal(m.shape); My = 0.2 * np.minimum(m[:-1], m[1:]) * rng.standard_normal((35, 72))
    q = np.full(m.shape, 3.7)
    q1, m1 = transportar(q, m, Mx, My)
    assert np.abs(q1 - 3.7).max() < 1e-12


def test_courant_zonal_mayor_que_uno():
    R = Rejilla(A_TIERRA, 36)
    m = R.area.copy()
    q = np.zeros_like(m); q[0, 10] = 1.0
    Mx = np.zeros_like(m); Mx[0] = 3.0 * m[0, 0]               # tres celdas por paso en la fila polar
    q1, m1 = transportar(q, m, Mx, np.zeros((35, 72)))
    assert abs(q1[0, 13] - 1.0) < 1e-12 and abs((m1 * q1).sum() - (m * q).sum()) < 1e-6

# fase6_trazadores.py -- Fase 6.3 (paso 3): transporte de trazadores (vapor) POSITIVO y CONSERVATIVO.
#
# Diseno (DISENO_FASE_6.3.md §1.5):
#   - Forma de flujo con las masas de aire de cada celda m y los flujos de masa por cara M (los del nucleo):
#     la masa de trazador se conserva exactamente y un campo constante sigue constante.
#   - Separacion por direcciones "consistente con la masa": en cada barrido 1D se actualizan a la vez la masa
#     del trazador (m q) y la del aire (m), y q = (m q) / m. Asi q = cte se conserva exactamente aunque cada
#     barrido por separado sea divergente (idea de Easter 1993 ⚠️ cita por verificar). Orden alterno
#     (XY, YX) en pasos sucesivos.
#   - Reconstruccion lineal a trozos con pendiente limitada (monotonizada, van Leer 1977 ⚠️): sin nuevos
#     maximos ni minimos en 1D -> q >= 0.
#   - Direccion zonal con numero de Courant arbitrario (celdas estrechas de cerca del polo): el flujo por una
#     cara suma las celdas enteras que la cruzan y la fraccion de la siguiente (idea "flux-form semi-
#     Lagrangian" de Lin y Rood 1996 ✅ cita, MWR 124, 2046). En la meridiana el Courant es siempre < 1.
#
# Convenios como el nucleo: m (..., F, C); M_x (..., F, C) masa que cruza la cara OESTE de cada celda hacia
# el este en el paso; M_y (..., F-1, C) masa que cruza la cara SUR de las filas 0..F-2 hacia el norte.

import numpy as np


def _pendiente_mc(q, eje):
    """Pendiente por celda (diferencia por celda) monotonizada (MC: minmod(2 dl, dc, 2 dr))."""
    if eje == -1:
        dl = q - np.roll(q, 1, axis=-1)
        dr = np.roll(q, -1, axis=-1) - q
    else:                                    # meridiano, sentido norte: "derecha" = fila i-1, "izquierda" = fila i+1
        dl = np.zeros_like(q); dr = np.zeros_like(q)          # en los polos (sin vecino) la pendiente es 0
        dr[..., 1:, :] = q[..., :-1, :] - q[..., 1:, :]      # q(fila i-1) - q(fila i)
        dl[..., :-1, :] = q[..., :-1, :] - q[..., 1:, :]     # q(fila i) - q(fila i+1)
    dc = 0.5 * (dl + dr)
    s = np.sign(dc)
    lim = np.minimum(np.minimum(2 * np.abs(dl), 2 * np.abs(dr)), np.abs(dc))
    return np.where(dl * dr > 0, s * lim, 0.0)


def _flujo_zonal(q, m, M):
    """Masa de trazador que cruza la cara oeste de cada celda (hacia el este) cuando cruza la masa de aire M,
    con M de cualquier tamano (celdas enteras + fraccion). Filas independientes y periodicas."""
    C = q.shape[-1]
    dq = _pendiente_mc(q, -1)
    F = np.zeros_like(M)
    # Recorre las celdas aguas arriba acumulando masa; como mucho C celdas (en la practica 1-3).
    resto = np.abs(M)
    pos = M > 0
    k = 0
    while np.any(resto > 0) and k < C:
        # celda donante k-esima aguas arriba: j-1-k si M > 0 (viene del oeste), j+k si M < 0 (del este)
        q_d = np.where(pos, np.roll(q, 1 + k, axis=-1), np.roll(q, -k, axis=-1))
        m_d = np.where(pos, np.roll(m, 1 + k, axis=-1), np.roll(m, -k, axis=-1))
        dq_d = np.where(pos, np.roll(dq, 1 + k, axis=-1), np.roll(dq, -k, axis=-1))
        toma = np.minimum(resto, m_d)
        c = np.where(m_d > 0, toma / np.where(m_d > 0, m_d, 1.0), 0.0)
        entera = toma >= m_d
        # fraccion c de la celda pegada a la cara: lado este de la donante si M > 0, lado oeste si M < 0
        q_frac = np.where(pos, q_d + 0.5 * (1 - c) * dq_d, q_d - 0.5 * (1 - c) * dq_d)
        q_med = np.where(entera, q_d, q_frac)
        F = F + np.sign(M) * toma * q_med
        resto = resto - toma
        k += 1
    return F


def barrido_x(mq, m, M_x):
    q = mq / m
    F = _flujo_zonal(q, m, M_x)
    mq_n = mq - (np.roll(F, -1, axis=-1) - F)
    m_n = m - (np.roll(M_x, -1, axis=-1) - M_x)
    return mq_n, m_n


def barrido_y(mq, m, M_y):
    q = mq / m
    dq = _pendiente_mc(q, -2)               # diferencia por celda en sentido norte (fila i-1 menos fila i)
    # cara sur de la fila i: separa la fila i (norte) de la i+1 (sur)
    m_s = m[..., 1:, :]; q_s = q[..., 1:, :]; dq_s = dq[..., 1:, :]          # celda sur (donante si M > 0)
    m_n = m[..., :-1, :]; q_n = q[..., :-1, :]; dq_n = dq[..., :-1, :]       # celda norte (donante si M < 0)
    c_s = np.clip(M_y / m_s, 0, 1)
    c_n = np.clip(-M_y / m_n, 0, 1)
    q_cara = np.where(M_y > 0, q_s + 0.5 * (1 - c_s) * dq_s, q_n - 0.5 * (1 - c_n) * dq_n)
    F = M_y * q_cara
    div_t = np.zeros_like(mq); div_m = np.zeros_like(m)
    div_t[..., :-1, :] -= F; div_t[..., 1:, :] += F
    div_m[..., :-1, :] -= M_y; div_m[..., 1:, :] += M_y
    return mq - div_t, m - div_m


def transportar(q, m, M_x, M_y, orden_xy=True):
    """Un paso de transporte horizontal. Devuelve (q_nuevo, m_nuevo). m_nuevo = m - div(M) exactamente."""
    mq = m * q
    if orden_xy:
        mq, m1 = barrido_x(mq, m, M_x)
        mq, m2 = barrido_y(mq, m1, M_y)
    else:
        mq, m1 = barrido_y(mq, m, M_y)
        mq, m2 = barrido_x(mq, m1, M_x)
    return mq / m2, m2


def _pendiente_vertical(q):
    """Pendiente por capa en sentido descendente (capa k+1 menos capa k), monotonizada; 0 en tope y fondo."""
    dl = np.zeros_like(q); dr = np.zeros_like(q)
    dl[1:] = q[1:] - q[:-1]          # q_k - q_{k-1}
    dr[:-1] = q[1:] - q[:-1]         # q_{k+1} - q_k
    dc = 0.5 * (dl + dr)
    lim = np.minimum(np.minimum(2 * np.abs(dl), 2 * np.abs(dr)), np.abs(dc))
    return np.where(dl * dr > 0, np.sign(dc) * lim, 0.0)


def barrido_z(mq, m, Mz):
    """Mz (N-1, ...): masa que cruza hacia ABAJO el seminivel entre la capa k y la k+1 en el paso.
    |Courant| < 1 (se comprueba)."""
    q = mq / m
    dq = _pendiente_vertical(q)
    arriba, abajo = slice(None, -1), slice(1, None)
    c_a = Mz / m[arriba]                         # baja: donante la capa de arriba (k)
    c_b = -Mz / m[abajo]                         # sube: donante la de abajo (k+1)
    assert np.all(np.where(Mz > 0, c_a, c_b) <= 1.0 + 1e-12), "Courant vertical > 1"
    q_cara = np.where(Mz > 0, q[arriba] + 0.5 * (1 - c_a) * dq[arriba], q[abajo] - 0.5 * (1 - c_b) * dq[abajo])
    F = Mz * q_cara
    mq = mq.copy(); m = m.copy()
    mq[arriba] -= F; mq[abajo] += F
    m[arriba] -= Mz; m[abajo] += Mz
    return mq, m


def transportar_3d(q, m, M_x, M_y, M_z, paso_par=True):
    """Transporte 3D (q, m: (N, F, C)). Orden simetrico alterno: X Y Z / Z Y X. Devuelve (q, m_nuevo)."""
    mq = m * q
    if paso_par:
        mq, m = barrido_x(mq, m, M_x); mq, m = barrido_y(mq, m, M_y); mq, m = barrido_z(mq, m, M_z)
    else:
        mq, m = barrido_z(mq, m, M_z); mq, m = barrido_y(mq, m, M_y); mq, m = barrido_x(mq, m, M_x)
    return mq / m, m

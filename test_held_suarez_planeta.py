# test_held_suarez_planeta.py -- M3N 3.20.0 (Fase 6.4, DISENO_FASE_6.4.md §8.4): Held y Suarez con la geometria
# de P3N (held_suarez.py --planeta p3n). Solo cambian el radio, la rotacion y la gravedad; el forzamiento, no.
#
# Uso:  python -m pytest test_held_suarez_planeta.py      (segundos)

import os

import numpy as np

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import parametros as P
from fase6_aguas_someras import A_TIERRA, OMEGA_TIERRA, G_TIERRA
from fase6_forzamientos import HeldSuarez
from fase6_nucleo import NucleoSeco
from fase30_multicapa import sigma_seminiveles
from held_suarez import constantes_planeta


def test_constantes_de_cada_planeta():
    """'tierra' son las de Williamson et al. (1992), como antes de la 3.20.0; 'p3n', las de parametros.py."""
    assert constantes_planeta("tierra") == (A_TIERRA, OMEGA_TIERRA, G_TIERRA)
    assert not P.MODO_TIERRA
    a, om, g = constantes_planeta("p3n")
    assert (a, om, g) == (P.P3N_RADIO, P.OMEGA_SIDERAL, P.P3N_GRAVEDAD)
    # los valores de DISENO_FASE_6.4.md §2
    assert abs(a / 1e3 - 4676.7) < 0.1 and abs(om - 8.824e-5) < 1e-8 and abs(g - 9.053) < 1e-3


def test_forzamiento_igual_en_los_dos_planetas():
    """HS94 no depende del planeta: misma T de equilibrio y mismos tiempos de relajacion (en segundos)."""
    sh = np.asarray(sigma_seminiveles(), dtype=float)
    hs = [HeldSuarez(NucleoSeco(*constantes_planeta(n), 287.0, 1004.0, sh)) for n in ("tierra", "p3n")]
    ps = np.full((36, 72), 1e5)
    assert np.array_equal(hs[0].T_equilibrio(ps), hs[1].T_equilibrio(ps))
    assert np.array_equal(hs[0].kT, hs[1].kT) and np.array_equal(hs[0].kv, hs[1].kv)


def test_p3n_estable_y_con_su_coriolis():
    """El nucleo con la geometria de P3N usa su rotacion y su radio, y 20 pasos de HS94 quedan finitos."""
    sh = np.asarray(sigma_seminiveles(), dtype=float)
    a, om, g = constantes_planeta("p3n")
    m = NucleoSeco(a, om, g, 287.0, 1004.0, sh)
    assert np.allclose(m.f_esq, 2 * om * np.sin(m.Rj.phi_f[1:-1])[:, None])
    assert m.Rj.a == a
    hs = HeldSuarez(m)
    m.preparar_semiimplicito(450.0)
    m.preparar_hiperdifusion(0.5)
    rng = np.random.default_rng(1)
    ps = np.full((36, 72), 1e5)
    T = hs.T_equilibrio(ps) + 0.1 * rng.standard_normal((m.N, 36, 72))
    fin = m.integrar_si(ps, T, np.zeros((m.N, 36, 72)), np.zeros((m.N, 35, 72)), 20, forzamiento=hs)
    assert all(np.isfinite(x).all() for x in fin)

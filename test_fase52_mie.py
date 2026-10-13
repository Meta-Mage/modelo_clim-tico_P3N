# test_fase52_mie.py -- la teoria de Mie propia (Fase 5.2, borrador).
import math

import numpy as np
import pytest

import fase52_mie as M

# Referencia: miepython 3.3 (efficiencies_mx, convenio m = n - i k), calculada el 12/10/2026.
# ((n, k), x) -> (Q_ext, Q_sca, g)
REFERENCIA = {
    ((1.33, 1e-08), 0.1): (1.10928809619e-05, 1.10906253621e-05, 0.00183195882088),
    ((1.33, 1e-08), 1.0): (0.0939240291967, 0.0939240005155, 0.184516674655),
    ((1.33, 1e-08), 10.0): (2.20654875442, 2.20654829916, 0.712459314548),
    ((1.33, 1e-08), 100.0): (2.10108983456, 2.10108502725, 0.868315509183),
    ((1.33, 1e-08), 1000.0): (2.01657862804, 2.01654442178, 0.883095885764),
    ((1.33, 0.01), 10.0): (2.24924090806, 1.87211206054, 0.754141067603),
    ((1.33, 0.1), 50.0): (2.13769696713, 1.1104088641, 0.96748401332),
    ((1.5, 0.5), 5.0): (2.53094027737, 1.18858883898, 0.862028001782),
    ((1.29, 0.0472), 30.0): (2.22036439176, 1.13833148323, 0.969053024349),
    ((1.55, 0.0), 3.0): (3.70220134746, 3.70220134746, 0.707863653071),
}


@pytest.mark.parametrize("caso,ref", list(REFERENCIA.items()))
def test_igual_que_miepython(caso, ref):
    (n, k), x = caso
    q = M.mie(complex(n, k), x)
    for a, b in zip(q, ref):
        assert a == pytest.approx(b, rel=1e-9, abs=1e-14)


def test_limite_de_rayleigh():
    """x << 1: Q_sca = 8/3 x^4 |(m^2-1)/(m^2+2)|^2 y Q_abs = 4 x Im[(m^2-1)/(m^2+2)] (Bohren y Huffman 1983, §5.2)."""
    m, x = complex(1.5, 0.01), 1e-3
    a = (m * m - 1) / (m * m + 2)
    qext, qsca, g = M.mie(m, x)
    assert qsca == pytest.approx(8.0 / 3.0 * x ** 4 * abs(a) ** 2, rel=1e-4)
    assert qext - qsca == pytest.approx(4.0 * x * a.imag, rel=1e-4)
    assert abs(g) < 1e-5


def test_paradoja_de_la_extincion():
    """Esferas muy grandes: Q_ext -> 2 (van de Hulst 1957)."""
    qext, _, _ = M.mie(complex(1.33, 1e-8), 5000.0)
    assert abs(qext - 2.0) < 0.01


def test_poblacion_de_gotas():
    """La gamma de Hansen y Travis tiene el radio efectivo pedido; en el visible una nube de r_e = 10 micras tiene
    Q_ext ~ 2, omega ~ 1 y g ~ 0,86-0,87 (lo tipico de las nubes de agua); a 3 micras el agua absorbe (omega ~ 0,5)."""
    r = np.linspace(1e-4, 60.0, 20000)
    n = M.distribucion_gamma(r, 10.0, 0.1)
    assert np.trapezoid(r ** 3 * n, r) / np.trapezoid(r ** 2 * n, r) == pytest.approx(10.0, rel=1e-6)
    q, w, g = M.optica_poblacion(complex(1.33, 1e-8), 0.55, 10.0, 0.1)
    assert 2.0 < q < 2.15 and w > 0.9999 and 0.85 < g < 0.88
    q, w, g = M.optica_poblacion(complex(1.29, 0.0472), 3.0, 10.0, 0.1)
    assert 0.45 < w < 0.6


def test_referencia_reproducible():
    mp = pytest.importorskip("miepython")
    for ((n, k), x), ref in REFERENCIA.items():
        r = mp.efficiencies_mx(complex(n, -k), x)
        assert (r[0], r[1], r[3]) == pytest.approx(ref, rel=1e-9)

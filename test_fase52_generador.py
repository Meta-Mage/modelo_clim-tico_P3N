# test_fase52_generador.py -- el generador de tablas k (Fase 5.2, borrador): el calculo rapido de k(nu) frente al de
# referencia (fase52_lineas, Voigt completo en una rejilla uniforme fina), el descarte de lineas y una casilla entera.
import numpy as np
import pytest

import fase52_generador as G
import fase52_k as K
import fase52_lineas as L


def _lista(rng, n, a, b):
    return {"nu0": np.sort(rng.uniform(a, b, n)), "S": 10 ** rng.uniform(-24, -19, n), "E": rng.uniform(0, 2000, n),
            "g_air": rng.uniform(0.05, 0.1, n), "g_self": rng.uniform(0.2, 0.5, n), "n_air": rng.uniform(0.5, 0.8, n),
            "d_air": rng.uniform(-0.005, 0.005, n)}


CASOS = [((300, 350), 1.0, 290.0), ((300, 350), 0.1, 250.0), ((300, 350), 1e-3, 210.0), ((300, 350), 1e-5, 200.0),
         ((1500, 1550), 1e-5, 200.0), ((1500, 1550), 0.03, 230.0)]


@pytest.mark.parametrize("pedestal", [False, True])
@pytest.mark.parametrize("banda,p,T", CASOS)
def test_rapido_igual_que_la_referencia(banda, p, T, pedestal):
    """150 lineas al azar (tambien fuera de la banda, que entran con sus alas), de la superficie a la alta atmosfera:
    la transmision de la banda con la rejilla adaptada y el calculo rapido difiere de la de referencia < 2e-4 (medido:
    hasta 1e-4), y la de la tabla k < 5e-4 (medido: hasta 3,3e-4, incluida la cuadratura)."""
    rng = np.random.default_rng(int(1000 * p) + int(T))
    a, b = banda
    lin = _lista(rng, 150, a - 30, b + 30)
    S, c, gL, aD = G.parametros(lin, T, p, 0.0, 1.0, 1.0, 18.0)
    nu = G.rejilla_fina(a, b, gL, aD)
    ke = L.coeficiente_absorcion(nu, lin, T, p, 0.0, 1.0, 1.0, 18.0, quitar_pedestal=pedestal)
    nua, w = G.rejilla_adaptada(a, b, c, gL, aD)
    ka = G.k_rapido(nua, S, c, gL, aD, quitar_pedestal=pedestal)
    assert abs(w.sum() - (b - a)) < 1e-9 and np.all(w > 0)
    _, wq, kgr = K.tabla_k(ke)
    _, _, kga = K.tabla_k(ka, dnu=w)
    for u in (1e18, 1e20, 1e22, 1e24):
        assert abs((np.exp(-ka * u) * w).sum() / w.sum() - np.exp(-ke * u).mean()) < 2e-4
        assert abs(K.transmision_ck(wq, kga, u) - K.transmision_ck(wq, kgr, u)) < 5e-4


def test_descarte_y_su_cota():
    """Las lineas descartadas absorben menos que la cota devuelta (suma de S u_max), medido en la banda."""
    rng = np.random.default_rng(7)
    lin = _lista(rng, 200, 600, 700)
    lin["S"] = 10 ** rng.uniform(-30, -20, 200)
    S, c, gL, aD = G.parametros(lin, 260.0, 0.3, 0.0, 1.0, 1.0, 44.0)
    u = 1e22
    quedan, cota = G.descartar(S, u, 1e-7)
    assert 0 < len(quedan) < 200 and cota > 0
    nu, w = G.rejilla_adaptada(620, 680, c, gL, aD)
    todo = G.k_rapido(nu, S, c, gL, aD)
    parte = G.k_rapido(nu, S[quedan], c[quedan], gL[quedan], aD[quedan])
    perdida = ((np.exp(-parte * u) - np.exp(-todo * u)) * w).sum()        # cm^-1 de absorcion perdida en la banda
    assert 0 <= perdida <= cota


def test_casilla_completa_con_dos_isotopologos():
    """k_gas_banda: sumar dos isotopologos (cada uno con su Q y su masa) da la misma tabla que hacerlo a mano."""
    rng = np.random.default_rng(9)
    lin = _lista(rng, 120, 640, 760)
    lin["iso"] = rng.choice([1, 2], 120)
    Ts = np.arange(60.0, 401.0)
    isotopos = {1: {"masa_u": 44.0, "T": Ts, "Q": 1.0 + Ts}, 2: {"masa_u": 45.0, "T": Ts, "Q": 2.0 + 2 * Ts}}
    kg, perdida, n, npts = G.k_gas_banda(lin, isotopos, (650, 750), 0.2, 240.0, 0.0, 1e22, 1e-12)
    partes = []
    for iso, inf in isotopos.items():
        m = lin["iso"] == iso
        sub = {k: v[m] for k, v in lin.items()}
        partes.append(G.parametros(sub, 240.0, 0.2, 0.0, 1.0 + 240.0 if iso == 1 else 2.0 + 480.0,
                                   1.0 + 296.0 if iso == 1 else 2.0 + 592.0, inf["masa_u"]))
    S, c, gL, aD = (np.concatenate([p[i] for p in partes]) for i in range(4))
    nu, w = G.rejilla_adaptada(650, 750, c, gL, aD)
    _, _, ref = K.tabla_k(G.k_rapido(nu, S, c, gL, aD), dnu=w)
    assert n == 120 and perdida == 0.0 and np.allclose(kg, ref, rtol=1e-12, atol=0)
    assert np.all(np.diff(kg) >= 0)                           # k(g) creciente


def test_banda_sin_lineas():
    lin = _lista(np.random.default_rng(1), 10, 5000, 5100)
    lin["iso"] = np.ones(10, dtype=int)
    kg, perdida, n, _ = G.k_gas_banda(lin, {1: {"masa_u": 18.0, "T": np.arange(60.0, 401.0),
                                                "Q": np.ones(341)}}, (500, 600), 1.0, 280.0, 0.0, 1e22, 1e-7)
    assert n == 0 and np.all(kg == 0)


def test_k_extra_se_suma_antes_de_ordenar():
    """Un continuo constante: con lineas, la tabla es la de las lineas mas la constante (k(g) se desplaza); sin
    lineas en la banda, la tabla es la constante."""
    rng = np.random.default_rng(4)
    lin = _lista(rng, 60, 640, 760)
    lin["iso"] = np.ones(60, dtype=int)
    iso = {1: {"masa_u": 18.0, "T": np.arange(60.0, 401.0), "Q": np.ones(341)}}
    base, _, _, _ = G.k_gas_banda(lin, iso, (650, 750), 0.5, 260.0, 0.0, 1e22, 1e-12)
    mas, _, _, _ = G.k_gas_banda(lin, iso, (650, 750), 0.5, 260.0, 0.0, 1e22, 1e-12,
                                 k_extra=lambda nu: np.full(len(nu), 1e-24))
    assert np.allclose(mas, base + 1e-24, rtol=1e-12, atol=1e-36)
    solo, _, n, _ = G.k_gas_banda(lin, iso, (5000, 5100), 0.5, 260.0, 0.0, 1e22, 1e-12,
                                  k_extra=lambda nu: np.full(len(nu), 2e-25))
    assert n == 0 and np.allclose(solo, 2e-25)

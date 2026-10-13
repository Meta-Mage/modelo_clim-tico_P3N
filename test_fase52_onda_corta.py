# test_fase52_onda_corta.py -- pruebas del dos flujos delta-Eddington (Fase 5.2, borrador).
import math

import numpy as np
import pytest

import fase52_onda_corta as O

# Referencia: ordenadas discretas con 32 direcciones (PythonicDISORT 1.8, funcion de fase de Henyey-Greenstein con
# escalado delta-M, suelo negro; calculado el 12/10/2026 en el entorno de la IA; con 16 direcciones cambia < 1e-4, asi
# que esta convergido). (tau, omega, g, mu0) -> (R, T total)
REFERENCIA = {
    (10, 0.999, 0.85, 1.0): (0.4133, 0.5665),     # nube, sol en lo alto
    (10, 0.999, 0.85, 0.5): (0.5940, 0.3864),     # nube, sol a 60 grados
    (10, 0.999, 0.85, 0.2): (0.7239, 0.2606),     # nube, sol bajo (el caso dificil del delta-Eddington)
    (1, 0.999, 0.85, 0.5): (0.1642, 0.8334),      # nube fina
    (30, 0.9999, 0.85, 0.5): (0.8014, 0.1932),    # nube gruesa
    (0.2, 0.9, 0.7, 0.5): (0.0576, 0.8991),       # aerosol
    (0.1, 0.999999, 0.0, 0.5): (0.0910, 0.9090),  # dispersion casi isotropa, sin absorcion (tipo Rayleigh)
    (5, 0.5, 0.5, 0.8): (0.0617, 0.0135),         # muy absorbente
}


def test_sin_dispersion_solo_beer():
    d = O.capa(2.0, 0.0, 0.5, 0.5)
    assert d["R"] == pytest.approx(0, abs=1e-12) and d["T_dif"] == pytest.approx(0, abs=1e-12)
    assert d["T_dir"] == pytest.approx(math.exp(-4.0), rel=1e-12)
    assert d["A_capa"] == pytest.approx(1 - math.exp(-4.0), rel=1e-12)


@pytest.mark.parametrize("albedo", [0.0, 0.3, 0.9])
@pytest.mark.parametrize("omega", [1.0, 0.999, 0.8])
def test_energia(albedo, omega):
    """Lo que entra sale al espacio o se absorbe (en la capa o en el suelo); sin absorcion, nada en la capa."""
    for tau, g, mu0 in ((0.1, 0.0, 0.9), (5.0, 0.85, 0.6), (50.0, 0.85, 0.3)):
        d = O.capa(tau, omega, g, mu0, albedo)
        assert d["R"] + d["A_capa"] + d["A_suelo"] == pytest.approx(1.0, abs=1e-9)
        assert d["A_capa"] >= -1e-6 and d["R"] >= 0 and d["T_dif"] >= -1e-12
        if omega == 1.0:
            assert abs(d["A_capa"]) < 1e-6


def test_escalado_delta():
    tau, omega, g = O.escalado_delta(10.0, 0.999, 0.85)
    f = 0.85 ** 2
    assert tau == pytest.approx((1 - 0.999 * f) * 10.0) and g == pytest.approx((0.85 - f) / (1 - f))
    assert omega == pytest.approx((1 - f) * 0.999 / (1 - 0.999 * f))


def test_mas_suelo_blanco_mas_reflexion():
    R = [O.capa(5.0, 0.999, 0.85, 0.5, a)["R"] for a in (0.0, 0.3, 0.6, 0.9)]
    assert all(b > a for a, b in zip(R, R[1:]))


@pytest.mark.parametrize("caso,ref", list(REFERENCIA.items()))
def test_frente_a_ordenadas_discretas(caso, ref):
    """El error conocido del delta-Eddington (King y Harshvardhan 1986 ⚠️): pequeño con el sol alto o capas finas y
    mayor con el sol bajo. Medido aqui: hasta 0,038 en la reflexion y la transmision (nube con mu0 = 0,2)."""
    d = O.capa(*caso)
    R, T = ref
    facil = caso[3] >= 0.5 and caso[0] <= 1.0
    tol = 0.02 if facil else 0.04
    assert abs(d["R"] - R) < tol and abs(d["T_dir"] + d["T_dif"] - T) < tol


def test_referencia_reproducible():
    """Si PythonicDISORT esta instalado (pip install PythonicDISORT), se recalcula la tabla de referencia."""
    pydisort = pytest.importorskip("PythonicDISORT").pydisort
    import warnings
    for (tau, omega, g, mu0), (R, T) in REFERENCIA.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, Fp, Fm = pydisort(np.array([tau]), np.array([omega]), 32, np.array([[g ** l for l in range(128)]]),
                                 mu0, 1.0, 0, NLeg=32, f_arr=g ** 32, only_flux=True)[:3]
        dif, dirc = Fm(tau)
        assert Fp(0.0) / mu0 == pytest.approx(R, abs=6e-5) and (dif + dirc) / mu0 == pytest.approx(T, abs=6e-5)


# ---------------------------------------------------------------- varias capas
# Referencia de varias capas (mismo calculo: PythonicDISORT 1.8, 32 direcciones; suelo de Lambert con el albedo dado).
# (capas: tau, omega, g; mu0; albedo) -> (R, T total, absorcion de cada capa)
REFERENCIA_COLUMNA = [
    (([0.05, 10, 0.3], [0.999999, 0.999, 0.85], [0.0, 0.85, 0.7], 0.6, 0.0), (0.5727, 0.3738, [0.0, 0.0203, 0.0331])),
    (([0.05, 10, 0.3], [0.999999, 0.999, 0.85], [0.0, 0.85, 0.7], 0.6, 0.2), (0.6052, 0.4120, [0.0, 0.0218, 0.0435])),
    (([0.1, 0.5, 2.0, 0.2], [0.999999, 0.95, 0.999, 0.9], [0.0, 0.7, 0.85, 0.7], 0.3, 0.15),
     (0.5211, 0.4367, [0.0, 0.0832, 0.0035, 0.0211])),
]


@pytest.mark.parametrize("caso", [(10, 0.999, 0.85, 0.5, 0.0), (2, 0.8, 0.6, 0.3, 0.4), (0.3, 0.0, 0.5, 0.7, 0.2)])
def test_columna_de_una_capa_es_la_capa(caso):
    a = O.capa(*caso)
    b = O.columna([caso[0]], [caso[1]], [caso[2]], caso[3], caso[4])
    for k in ("R", "T_dir", "T_dif", "A_suelo"):
        assert b[k] == pytest.approx(a[k], abs=1e-12)
    assert b["A_capas"][0] == pytest.approx(a["A_capa"], abs=1e-12)


def test_dividir_una_capa_no_cambia_nada():
    """La solucion es exacta para las ecuaciones: partir una capa homogenea en cinco da lo mismo (hasta el redondeo)."""
    una = O.columna([10.0], [0.99], [0.85], 0.5, 0.3)
    cinco = O.columna([2.0] * 5, [0.99] * 5, [0.85] * 5, 0.5, 0.3)
    assert cinco["R"] == pytest.approx(una["R"], abs=1e-12)
    assert cinco["T_dif"] == pytest.approx(una["T_dif"], abs=1e-12)
    assert cinco["A_capas"].sum() == pytest.approx(una["A_capas"][0], abs=1e-12)
    assert np.all(np.diff(cinco["abajo_directo"]) < 0)


def test_columna_energia_y_estabilidad():
    """Muchas capas y una nube muy gruesa (tau 100) sin desbordamientos; la energia se conserva."""
    d = O.columna([2.5] * 40, [0.9999] * 40, [0.85] * 40, 0.4, 0.6)
    assert np.isfinite(d["R"]) and 0.9 < d["R"] < 1
    assert d["R"] + d["A_capas"].sum() + d["A_suelo"] == pytest.approx(1.0, abs=1e-9)
    d = O.columna([0.05, 10, 0.3], [1.0, 0.999, 0.85], [0.0, 0.85, 0.7], 0.6, 0.2)
    assert d["R"] + d["A_capas"].sum() + d["A_suelo"] == pytest.approx(1.0, abs=1e-9)
    assert abs(d["A_capas"][0]) < 1e-6 and np.all(d["A_capas"] > -1e-6)


@pytest.mark.parametrize("caso,ref", REFERENCIA_COLUMNA)
def test_columna_frente_a_ordenadas_discretas(caso, ref):
    """Mismo margen que una capa (el error del delta-Eddington): 0,04 en la reflexion y la transmision; 0,02 en la
    absorcion de cada capa. Medido: hasta 0,034 (transmision con mu0 = 0,3) y 0,013 (absorcion)."""
    d = O.columna(*caso)
    R, T, A = ref
    assert abs(d["R"] - R) < 0.04 and abs(d["T_dir"] + d["T_dif"] - T) < 0.04
    assert np.all(np.abs(d["A_capas"] - np.array(A)) < 0.02)


def test_referencia_columna_reproducible():
    pydisort = pytest.importorskip("PythonicDISORT").pydisort
    import warnings
    for (taus, om, gg, mu0, alb), (R, T, A) in REFERENCIA_COLUMNA:
        leg = np.array([[g ** l for l in range(128)] for g in gg])
        kw = {"BDRF_Fourier_modes": [alb]} if alb else {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, Fp, Fm = pydisort(np.cumsum(taus), np.array(om), 32, leg, mu0, 1.0, 0, NLeg=32,
                                 f_arr=np.array(gg) ** 32, only_flux=True, **kw)[:3]
        b = np.concatenate([[0], np.cumsum(taus)])
        up = Fp(b) / mu0
        dif, dirc = Fm(b)
        net = (dif + dirc) / mu0 - up
        assert up[0] == pytest.approx(R, abs=6e-5) and (dif[-1] + dirc[-1]) / mu0 == pytest.approx(T, abs=6e-5)
        assert np.allclose(net[:-1] - net[1:], A, atol=6e-5)

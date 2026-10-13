# test_fase52_rayleigh.py -- la dispersion de Rayleigh gas a gas (Fase 5.2, T22).
import numpy as np
import pytest

import fase52_rayleigh as R


def _peck_reeder_300ppm(lam):
    """Indice del aire seco con 300 ppm de CO2 a 15 °C y 1013,25 hPa (Peck y Reeder 1972; Bodhaine et al. 1999, ec. 18)."""
    return (8060.51 + 2480990.0 / (132.274 - lam ** -2) + 17455.7 / (39.32957 - lam ** -2)) * 1e-8


@pytest.mark.parametrize("lam", [0.25, 0.3, 0.4, 0.55, 0.8, 1.0, 1.5])
def test_la_mezcla_de_la_tierra_reproduce_el_indice_del_aire(lam):
    """Sumando la polarizabilidad de cada gas (Lorentz-Lorenz, aqui en el limite n ~ 1: refractividad proporcional a
    la densidad) con la composicion de la Tierra sale el indice del aire de Peck y Reeder dentro del 0,1 %."""
    N15 = 101325.0 / (R.K_B * 288.15)
    x = R.composicion(0.20946, 300e-6)
    mezcla = 0.0
    for g, f in x.items():
        nm1, N = R.refractividad(g, lam)
        mezcla += f * nm1 * N15 / N
    assert mezcla == pytest.approx(_peck_reeder_300ppm(lam), rel=1e-3)


def test_seccion_y_espesor_a_550_nm_como_bodhaine():
    """Bodhaine et al. 1999, tabla 3 (aire seco, 360 ppm de CO2): sigma(0,55 µm) = 4,5105e-27 cm^2 y espesor optico
    a nivel del mar a 45° ~0,0971 (medido aqui: 4,5112e-27 y 0,0969; la diferencia del espesor es la gravedad y la
    masa molecular que ellos ajustan a la altura y la latitud)."""
    x = R.composicion(0.20946, 360e-6)
    assert R.seccion_aire(0.55, x) == pytest.approx(4.5105e-27, rel=1e-3)
    assert R.espesor_optico(0.55, 101325.0, 9.80616, x) == pytest.approx(0.0971, rel=5e-3)


def test_o2_continuo_entre_la_tabla_y_la_formula():
    a, _ = R.refractividad("O2", R.LAMBDA_SMITH_MAX * (1 - 1e-9))
    b, _ = R.refractividad("O2", R.LAMBDA_SMITH_MAX)
    # la tabla esta a 0 °C y la formula a 20 °C: se comparan por molecula (refractividad / densidad)
    na = a / (101325.0 / (R.K_B * 273.15))
    nb = b / (101325.0 / (R.K_B * 293.15))
    assert na == pytest.approx(nb, rel=2e-3)


def test_casi_lambda_menos_4_y_king():
    """Entre 0,4 y 0,8 µm la seccion cae algo mas que lambda^-4 (la dispersion del indice: ~lambda^-4,08);
    Ar no se despolariza (F = 1) y N2 y O2 si (F > 1)."""
    x = R.composicion(0.20946, 400e-6)
    r = R.seccion_aire(0.4, x) / R.seccion_aire(0.8, x)
    assert 16.0 < r < 2 ** 4.15
    assert R.factor_king("Ar", 0.5) == 1.0 and R.factor_king("N2", 0.5) > 1.03 and R.factor_king("O2", 0.5) > 1.09


def test_aire_de_p3n():
    """Con la composicion sorteada de P3N (18,88 % de O2, 174,7 ppm de CO2), la seccion por molecula apenas cambia:
    +0,1 a +0,3 % (el N2 dispersa algo mas que el O2 al que sustituye)."""
    t = R.seccion_aire(np.array([0.3, 0.55, 1.0]), R.composicion(0.20946, 280e-6))
    p = R.seccion_aire(np.array([0.3, 0.55, 1.0]), R.composicion(0.1888, 174.7e-6))
    assert np.all(p / t > 1.0005) and np.all(p / t < 1.004)


def test_gas_desconocido():
    with pytest.raises(ValueError):
        R.factor_king("CH4", 0.5)

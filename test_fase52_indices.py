# test_fase52_indices.py -- indices de refraccion del agua y del hielo (T17) y optica de las gotas por banda (T18),
# Fase 5.2, borrador.
import hashlib
import os

import numpy as np
import pytest

import fase52_indices as I
import fase52_optica_nubes as O


def test_archivo_de_datos_es_el_preparado():
    """El archivo del repositorio es exactamente el que escribio herramientas/preparar_indices_refraccion.py el
    13/10/2026 con la base de refractiveindex.info fijada (sus huellas, en ese programa)."""
    with open(I.RUTA, "rb") as f:
        assert hashlib.sha256(f.read()).hexdigest() == \
            "8841dda15097878cfc730abbc7337eadaa952f0a3c51eeb3a3ddc84a3dc1aa93"


# Tabla impresa de Warren y Brandt 2008 (atmos.washington.edu/ice_optical_constants/IOPtable_2008.pdf, leida el
# 13/10/2026): lambda (µm), n, k.
TABLA_WB = [(0.5, 1.3130, 5.889e-10), (1.0, 1.3015, 1.620e-6), (2.0, 1.2744, 1.640e-3), (3.003, 1.0390, 4.380e-1),
            (3.096, 1.4225, 5.930e-1), (10.0, 1.1926, 5.008e-2), (11.9, 1.2546, 4.090e-1), (50.03, 1.7505, 6.191e-1),
            (100.0, 1.8654, 1.706e-1)]


@pytest.mark.parametrize("lam,n,k", TABLA_WB)
def test_hielo_igual_que_la_tabla_de_warren_y_brandt(lam, n, k):
    m = I.indice_hielo(lam)
    assert m.real == pytest.approx(n, rel=1e-12)
    assert m.imag == pytest.approx(k, rel=1e-9)


def test_agua_en_los_nodos_es_el_dato():
    """En un punto de la rejilla de Rowe et al. y a una de sus temperaturas, el valor es el de la tabla; a 25 °C y en
    un punto de Hale y Querry, el de Hale y Querry."""
    d = I.datos()
    for i, T in enumerate(d["rowe_T"]):
        for j in (0, 5000, 12000, len(d["rowe_lam"]) - 1):
            m = I.indice_agua(d["rowe_lam"][j], T)
            assert m.real == pytest.approx(d["rowe_n"][i, j], rel=1e-12)
            assert m.imag == pytest.approx(d["rowe_k"][i, j], rel=1e-9)
    for j in range(len(d["hq_lam"])):
        m = I.indice_agua(d["hq_lam"][j], I.T_HQ)
        assert m.real == pytest.approx(d["hq_n"][j], rel=1e-12)
        assert m.imag == pytest.approx(d["hq_k"][j], rel=1e-9)


def test_agua_lineal_en_la_temperatura_y_extremos():
    lam = np.array([1.0, 3.0, 10.0, 12.0, 20.0, 150.0])
    a, b = I.indice_agua(lam, 253.0), I.indice_agua(lam, 263.0)
    c = I.indice_agua(lam, 258.0)
    assert np.allclose(c.real, 0.5 * (a.real + b.real), rtol=1e-12)
    assert np.allclose(np.log(c.imag), 0.5 * (np.log(a.imag) + np.log(b.imag)), rtol=1e-12)
    assert np.array_equal(I.indice_agua(lam, 200.0), I.indice_agua(lam, 240.0))    # por debajo: el de 240 K
    assert np.array_equal(I.indice_agua(lam, 330.0), I.indice_agua(lam, I.T_HQ))   # por encima: el de 298 K
    assert np.array_equal(I.indice_agua(500.0, 290.0), I.indice_agua(500.0, 273.0))  # > 200 µm: solo Rowe


def test_fuera_de_los_datos_da_error():
    with pytest.raises(ValueError):
        I.indice_agua(0.1, 280.0)
    with pytest.raises(ValueError):
        I.indice_hielo(0.01)


def test_cambios_en_las_costuras_medidos():
    """Donde se pasa de un conjunto a otro, el salto es pequeño (medido el 13/10/2026):
    - a 0,667 µm (Hale y Querry -> Rowe et al.): n, < 0,005; k (~2·10^-8, el agua casi no absorbe), < 25 %;
    - a 200 µm y por encima de 273 K (fin de Hale y Querry): hasta 0,13 en n y 0,14 en k, pero por debajo de 50 cm^-1
      esta menos del 0,1 % de la emision termica (Planck a 260-300 K)."""
    d = I.datos()
    l0 = d["rowe_lam"][0]
    for T in (240.0, 273.0, 298.0):
        a, b = I.indice_agua(l0 * (1 - 1e-9), T), I.indice_agua(l0, T)
        assert abs(a.real - b.real) < 0.005
        assert 0.8 < b.imag / a.imag < 1.25
    for T in (240.0, 263.0, 273.0):                        # sin nodo de 298 K: continuo
        a, b = I.indice_agua(200.0, T), I.indice_agua(200.0 * (1 + 1e-9), T)
        assert abs(a - b) < 1e-6
    a, b = I.indice_agua(200.0, I.T_HQ), I.indice_agua(200.0 * (1 + 1e-9), I.T_HQ)
    assert abs(a.real - b.real) < 0.13 and abs(a.imag - b.imag) < 0.14
    nu = np.linspace(0.01, 5000, 500001)
    for T in (260.0, 300.0):
        B = O.planck_nu(nu, T)
        assert np.trapezoid(B[nu < 50], nu[nu < 50]) / np.trapezoid(B, nu) < 1e-3


def test_k_del_agua_crece_al_enfriarse_en_la_ventana():
    """Lo que motiva usar Rowe et al. 2020: a 12 µm, k es ~50 % mayor a 240 K que a 298 K (0,31 frente a 0,20)."""
    assert I.indice_agua(12.0, 240.0).imag == pytest.approx(0.3075, abs=0.002)
    assert I.indice_agua(12.0, I.T_HQ).imag == pytest.approx(0.199, abs=1e-6)


# --------------------------------------------------------------------------------------- optica de las gotas (T18)
def test_visible_extincion_de_stephens():
    """En el visible, Q_ext ~ 2, omega ~ 1 y g ~ 0,86 (nube de r_e = 10 µm): beta = 3 Q_ext/(4 rho r_e) ~ 157 m^2/kg."""
    beta, omega, g = O.optica_monocromatica_agua(1e4 / 0.55, 10.0, 280.0)
    assert beta == pytest.approx(3 * 2.09 / (4 * 1000 * 10e-6), rel=0.01)
    assert omega > 0.99999
    assert 0.85 < g < 0.88


def test_promedio_de_banda_con_los_pesos_de_T18(monkeypatch):
    """Con propiedades conocidas en cada punto, el promedio es el de las formulas de T18."""
    vals = {}

    def falso(nu, r_e, T, v_e, n_radios):
        b, w, g = 100.0 + nu, 0.5 + 1e-4 * nu, 0.8 + 1e-5 * nu
        vals[nu] = (b, w, g)
        return b, w, g
    monkeypatch.setattr(O, "optica_monocromatica_agua", falso)
    beta, omega, g = O.optica_banda_agua(800.0, 1000.0, 10.0, 260.0, n_nu=4)
    nu = O.puntos_banda(800.0, 1000.0, 4)
    S = O.planck_nu(nu, 260.0)
    b = np.array([vals[x][0] for x in nu]); w = np.array([vals[x][1] for x in nu]); gg = np.array([vals[x][2] for x in nu])
    assert beta == pytest.approx(np.sum(S * b) / np.sum(S), rel=1e-14)
    assert omega == pytest.approx(np.sum(S * b * w) / np.sum(S * b), rel=1e-14)
    assert g == pytest.approx(np.sum(S * b * w * gg) / np.sum(S * b * w), rel=1e-14)
    with pytest.raises(ValueError):
        O.optica_banda_agua(800.0, 1000.0, 10.0, 260.0, peso=lambda x: -np.ones_like(x), n_nu=4)


def test_efecto_de_la_temperatura_en_las_bandas_de_rowe():
    """Rowe et al. 2020 (resumen): ignorar la temperatura sesga los flujos de las nubes subfundidas hacia arriba entre
    250 y 580 cm^-1 y hacia abajo entre 710 y 940 cm^-1; es decir, el agua fria absorbe menos en la primera region
    y mas en la segunda. Con la optica por banda: omega a 240 K es mayor que a 298 K en 400-550 cm^-1 y menor en
    820-980 cm^-1 (medido: 0,517 frente a 0,445 y 0,422 frente a 0,476, con r_e = 10 µm)."""
    a = O.optica_banda_agua(400.0, 550.0, 10.0, 240.0, n_nu=6, n_radios=200)
    b = O.optica_banda_agua(400.0, 550.0, 10.0, 298.0, n_nu=6, n_radios=200)
    assert a[1] > b[1] + 0.05
    a = O.optica_banda_agua(820.0, 980.0, 10.0, 240.0, n_nu=6, n_radios=200)
    b = O.optica_banda_agua(820.0, 980.0, 10.0, 298.0, n_nu=6, n_radios=200)
    assert a[1] < b[1] - 0.03


# ------------------------------------------------------------------------------------- tabla de las gotas y radio
def test_tabla_de_las_gotas_en_el_infrarrojo():
    """La tabla del repositorio (herramientas/generar_optica_gotas.py): en sus nodos, TablaGotas devuelve lo
    guardado; una casilla recalculada da lo mismo; y los valores son fisicos."""
    T = O.TablaGotas()
    assert len(T.bandas) == 16 and T.beta.shape == (16, len(T.r_e), len(T.T))
    b, w, g = T.propiedades(T.r_e[5], T.T[3])
    assert np.allclose(b, T.beta[:, 5, 3], rtol=1e-12) and np.allclose(w, T.omega[:, 5, 3]) and np.allclose(g, T.g[:, 5, 3])
    ref = O.optica_banda_agua(*T.bandas[6], float(T.r_e[5]), float(T.T[3]), n_nu=20, n_radios=400)
    assert np.allclose([T.beta[6, 5, 3], T.omega[6, 5, 3], T.g[6, 5, 3]], ref, rtol=1e-12)
    assert np.all((T.omega > 0) & (T.omega < 1)) and np.all((T.g > 0) & (T.g < 1)) and np.all(T.beta > 0)
    # entre nodos, entre los valores vecinos; fuera de la rejilla, el extremo
    b2 = T.propiedades(np.sqrt(T.r_e[5] * T.r_e[6]), T.T[3])[0]
    lo, hi = np.minimum(T.beta[:, 5, 3], T.beta[:, 6, 3]), np.maximum(T.beta[:, 5, 3], T.beta[:, 6, 3])
    assert np.all(b2 >= lo * 0.98) and np.all(b2 <= hi * 1.02)
    assert np.allclose(T.propiedades(100.0, 400.0)[1], T.omega[:, -1, -1])


def test_radio_efectivo_de_martin_con_las_constantes_del_ifs():
    assert O.gotas_por_cm3(False) == pytest.approx(50.575) and O.gotas_por_cm3(True) == pytest.approx(313.2)
    # 0,3 g/m^3: (3 * 3e-7 / (4 pi 50,575 * 0,77))^(1/3) cm = 12,25 µm
    assert O.radio_efectivo_um(0.3, False) == pytest.approx(12.25, abs=0.01)
    assert O.radio_efectivo_um(0.3, True) == pytest.approx(6.92, abs=0.01)
    assert O.radio_efectivo_um(0.3, False, rwc_gm3=0.1) > O.radio_efectivo_um(0.3, False)
    assert O.radio_efectivo_um(1e-6, False) == O.RE_MIN_UM and O.radio_efectivo_um(50.0, False) == O.RE_MAX_UM

# test_fase52_continuo.py -- el continuo del vapor MT_CKD 4.3 (Fase 5.2, borrador).
import math
import os

import numpy as np
import pytest

import fase52_continuo as C

RUTA = os.environ.get("M3N_MTCKD", os.path.expanduser("~/Documentos/B3N/datos/hitran/absco-ref_wv-mt-ckd.nc"))

# Salida de referencia del propio MT_CKD (run_example/WATER.COEF_mt_ckd_AER del repositorio de AER, commit 679f0a4):
# coeficientes a 296 K y la densidad de referencia, con el termino de radiacion (propio, ajeno), cm^2/molecula.
REFERENCIA = {100.0: (2.3555e-20, 2.5328e-21), 1000.0: (1.3107e-22, 2.4133e-25), 2500.0: (8.5275e-24, 1.8283e-27),
              5000.0: (1.3225e-23, 8.3500e-26), 10000.0: (3.8380e-26, 1.0160e-27)}


def test_termino_de_radiacion():
    T = 250.0
    xkt = T / C.RADCN2
    nu = np.array([0.001, 100.0, 5000.0])
    r = C.termino_radiacion(nu, T)
    assert r[0] == pytest.approx(0.5 * nu[0] / xkt * nu[0])                    # x <= 0,01
    e = math.exp(-100.0 / xkt)
    assert r[1] == pytest.approx(100.0 * (1 - e) / (1 + e))
    assert r[2] == 5000.0                                                      # x > 10


def test_interpolacion_exacta_en_los_nodos_y_para_rectas():
    rej = np.arange(-20.0, 200.0, 10.0)
    a = 3.0 + 0.5 * rej
    assert np.allclose(C._interpolar(rej[2:-3], rej, a), a[2:-3])
    nu = np.linspace(5.0, 150.0, 77)
    assert np.allclose(C._interpolar(nu, rej, a), 3.0 + 0.5 * nu)


@pytest.mark.skipif(not os.path.exists(RUTA), reason="sin el archivo de MT_CKD (lo descarga descargar_hitran.py)")
def test_igual_que_mt_ckd():
    d = C.leer_mtckd(RUTA)
    for nu, (s_ref, f_ref) in REFERENCIA.items():
        s, _ = C.coeficientes(np.array([nu]), 296.0, 1013.0, 1.0, d)     # vmr 1: el propio sin escalar
        _, f = C.coeficientes(np.array([nu]), 296.0, 1013.0, 0.0, d)     # vmr 0: el ajeno sin escalar
        assert s[0] == pytest.approx(s_ref, rel=1e-4) and f[0] == pytest.approx(f_ref, rel=1e-4)
    # escalado con la densidad y la fraccion molar
    s1, f1 = C.coeficientes(np.array([1000.0]), 250.0, 500.0, 0.01, d)
    s0, f0 = C.coeficientes(np.array([1000.0]), 250.0, 1013.0, 1.0, d)
    _, f00 = C.coeficientes(np.array([1000.0]), 250.0, 1013.0, 0.0, d)
    rho = (500.0 / 1013.0) * (296.0 / 250.0) / ((1013.0 / 1013.0) * (296.0 / 250.0))
    assert s1[0] == pytest.approx(s0[0] * 0.01 * rho) and f1[0] == pytest.approx(f00[0] * 0.99 * rho)

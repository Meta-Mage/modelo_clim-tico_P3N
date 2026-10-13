# test_fase52_hitran.py -- el lector propio de HITRAN y las bandas (Fase 5.2, borrador).
import numpy as np
import pytest

import fase52_hitran as H

# Cinco lineas con el formato .par (los valores imitan lineas reales; la ultima, del isotopologo 11 de CO2, 'A', y con
# la intensidad en el formato raro sin "E" que tambien acepta hitran-api).
FILAS = [
    (2, "1", 667.379960, "1.050E-19", "1.567E+00", 0.0748, 0.095, 0.0000, 0.75, -0.002100),
    (2, "2", 648.478600, "2.660E-21", "1.538E+00", 0.0753, 0.096, 1.5626, 0.74, -0.001500),
    (1, "1", 1594.746300, "1.130E-19", "2.070E+01", 0.0950, 0.433, 23.7944, 0.68, -0.007760),
    (1, "0", 1600.000000, "1.000E-25", "1.000E-02", 0.0800, 0.300, 300.0000, 0.70, 0.000000),
    (2, "A", 2300.000000, "2.700-164", "1.000E-03", 0.0700, 0.090, 3000.0000, 0.70, -0.001000),
]


def _fortran(x, ancho, dec):
    """Como Fortran Fw.d: si no cabe, sin el cero de delante (".0748", "-.002100"), como en los archivos de HITRAN."""
    s = f"{x:{ancho}.{dec}f}"
    if len(s) > ancho:
        s = s.replace("0.", ".", 1)
    assert len(s) == ancho
    return s


def _linea(f):
    m, i, nu, S, A, ga, gs, E, n, d = f
    l = (f"{m:2d}{i}{nu:12.6f}{S:>10s}{A:>10s}{_fortran(ga, 5, 4)}{_fortran(gs, 5, 3)}{E:10.4f}{_fortran(n, 4, 2)}"
         f"{_fortran(d, 8, 6)}")
    return l + " " * (146 - len(l)) + f"{3.0:7.1f}{1.0:7.1f}"      # pesos estadisticos g', g'' (columnas 146-159)


def _escribir(tmp_path):
    ruta = tmp_path / "prueba.data"
    ruta.write_text("\n".join(_linea(f) for f in FILAS) + "\n", encoding="ascii")
    return ruta


def test_lector_lee_las_columnas_y_los_codigos_de_isotopologo(tmp_path):
    d = H.leer_par(_escribir(tmp_path))
    assert list(d["molec"]) == [2, 2, 1, 1, 2] and list(d["iso"]) == [1, 2, 1, 10, 11]
    assert d["nu0"][0] == 667.37996 and d["S"][0] == 1.050e-19 and d["S"][4] == 2.7e-164
    assert d["g_air"][2] == 0.095 and d["g_self"][2] == 0.433 and d["E"][2] == 23.7944
    assert d["n_air"][1] == 0.74 and d["d_air"][2] == -0.00776
    sub = H.leer_par(_escribir(tmp_path), 600, 1000)
    assert len(sub["nu0"]) == 2 and set(sub["iso"]) == {1, 2}
    assert len(H.lineas_de_isotopologo(d, 1)["nu0"]) == 2


def test_lector_igual_que_hitran_api(tmp_path):
    """Mismos valores, linea a linea, que el lector oficial (si hitran-api esta instalado)."""
    hapi = pytest.importorskip("hapi")
    import contextlib
    import io
    import shutil
    carpeta = tmp_path / "db"
    carpeta.mkdir()
    shutil.copy(_escribir(tmp_path), carpeta / "T.data")
    import copy
    import json
    cab = copy.deepcopy(hapi.HITRAN_DEFAULT_HEADER)   # la cabecera del formato .par, como la que guarda hitran-api
    cab["table_name"], cab["number_of_rows"] = "T", len(FILAS)
    (carpeta / "T.header").write_text(json.dumps(cab), encoding="ascii")
    with contextlib.redirect_stdout(io.StringIO()):
        hapi.db_begin(str(carpeta))
    ref = hapi.LOCAL_TABLE_CACHE["T"]["data"]
    d = H.leer_par(carpeta / "T.data")
    pares = {"molec": "molec_id", "iso": "local_iso_id", "nu0": "nu", "S": "sw", "A": "a", "g_air": "gamma_air",
             "g_self": "gamma_self", "E": "elower", "n_air": "n_air", "d_air": "delta_air"}
    for mio, suyo in pares.items():
        assert np.array_equal(np.asarray(d[mio]), np.asarray(ref[suyo])), mio


def test_bandas_contiguas_y_ordenadas():
    for bandas in (H.BANDAS_IR, H.BANDAS_LUZ[:-1]):
        assert all(a < b for a, b in bandas) and all(b == c for (_, b), (c, _) in zip(bandas, bandas[1:]))
    assert H.BANDAS_IR[0][0] == 10 and H.BANDAS_IR[-1][1] == 3250 and H.BANDAS_LUZ[-2][1] == 50000
    assert len(H.BANDAS_IR) == 16 and len(H.BANDAS_LUZ) == 14

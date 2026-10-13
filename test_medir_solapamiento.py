# test_medir_solapamiento.py -- la herramienta que mide el error del solapamiento de gases (Fase 5.2, T9), con lineas
# sinteticas (las de verdad estan en el PC de Carlos).
import json
import os
import subprocess
import sys

import numpy as np
import pytest

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "herramientas"))
import medir_solapamiento as M     # noqa: E402


def _lineas(rng, n, a, b, s0):
    return {"nu0": np.sort(rng.uniform(a, b, n)), "S": s0 * rng.lognormal(0, 2, n), "E": rng.uniform(0, 2000, n),
            "g_air": rng.uniform(0.05, 0.1, n), "g_self": rng.uniform(0.2, 0.5, n), "n_air": rng.uniform(0.6, 0.8, n),
            "d_air": rng.uniform(-0.005, 0.0, n), "iso": np.ones(n, dtype=int), "molec": np.ones(n, dtype=int)}


def _iso(masa):
    T = np.arange(60.0, 401.0)
    return {1: {"masa_u": masa, "T": T, "Q": 100.0 * (T / 296.0) ** 1.5}}


@pytest.fixture(scope="module")
def datos():
    rng = np.random.default_rng(11)
    return {"H2O": (_lineas(rng, 150, 700, 820, 3e-21), _iso(18.0)),
            "CO2": (_lineas(rng, 150, 690, 830, 3e-20), _iso(44.0))}


def test_un_solo_gas_rorr_y_suma_son_la_tabla_del_gas(datos):
    col = M.columna(6)
    res = M.medir_banda({"CO2": datos["CO2"]}, (700, 820), col, None)
    for m in ("RORR", "SUMA"):
        assert np.allclose(res[m][0], res["PM"][0], rtol=1e-12) and np.allclose(res[m][1], res["PM"][1], rtol=1e-12)


def test_errores_pequenos_y_finitos_con_dos_gases(datos):
    """Con lineas sinteticas, las tablas k quedan cerca del linea a linea (< 0,5 W/m2 en una banda de ~20 W/m2) y
    sumar los gases como si estuvieran correlacionados es lo que mas se equivoca."""
    col = M.columna(6)
    res = M.medir_banda(datos, (700, 820), col, None)
    r = M.resumen(res, col)
    assert 5.0 < r["LBL"]["OLR"] < 40.0
    for m, v in r.items():
        assert all(np.isfinite(x) for x in v.values())
    for m in ("PM", "PM+fr", "RORR", "MAPA", "MAPA-log", "MAPA-tr"):
        assert abs(r[m]["dOLR"]) < 0.5 and abs(r[m]["dDLR"]) < 0.5, m
    assert abs(r["SUMA"]["dDLR"]) > abs(r["RORR"]["dDLR"])


def test_mapa_de_otra_columna(datos):
    """El mapa hecho con la columna de referencia se puede usar en otra (la polar): la herramienta lo hace sin errores
    y el linea a linea de la polar no depende del mapa."""
    ref, pol = M.columna(6, Ts=288.0), M.columna(6, Ts=255.0)
    rej = M.rejilla_comun(datos, (700, 820), ref)
    a = M.medir_banda(datos, (700, 820), pol, None, col_ref=ref, rejilla=rej)
    b = M.medir_banda(datos, (700, 820), pol, None, rejilla=rej)
    assert np.allclose(a["LBL"][0], b["LBL"][0]) and np.allclose(a["PM"][0], b["PM"][0])


def test_programa_completo(tmp_path, datos):
    """El programa lee una carpeta como la de descargar_hitran.py y escribe informe y resultados."""
    from test_fase52_hitran import _fortran
    carpeta = tmp_path / "hitran"
    carpeta.mkdir()
    tips = {"T_K": list(range(60, 401)), "Q": {}}
    iso = {}
    for gas, m, masa in (("H2O", 1, 18.010565), ("CO2", 2, 43.98983)):
        L = datos[gas][0]
        filas = []
        for i in range(len(L["nu0"])):
            l = (f"{m:2d}1{L['nu0'][i]:12.6f}{L['S'][i]:10.3E}{1.0:10.3E}{_fortran(L['g_air'][i], 5, 4)}"
                 f"{_fortran(L['g_self'][i], 5, 3)}{L['E'][i]:10.4f}{_fortran(L['n_air'][i], 4, 2)}"
                 f"{_fortran(L['d_air'][i], 8, 6)}")
            filas.append(l + " " * (146 - len(l)) + f"{3.0:7.1f}{1.0:7.1f}")
        (carpeta / f"{gas}.data").write_text("\n".join(filas) + "\n", encoding="ascii")
        tips["Q"][f"{m}_1"] = [100.0 * (t / 296.0) ** 1.5 for t in range(60, 401)]
        iso[f"{m}_1"] = {"masa_u": masa}
    (carpeta / "tips.json").write_text(json.dumps(tips))
    (carpeta / "isotopologos.json").write_text(json.dumps(iso))
    salida = tmp_path / "salida"
    r = subprocess.run([sys.executable, os.path.join(AQUI, "herramientas", "medir_solapamiento.py"), "--hitran",
                        str(carpeta), "--bandas", "4", "3", "--hilos", "2", "--gases", "H2O", "CO2", "--capas", "5", "--columnas", "ref",
                        "polar", "--salida", str(salida)], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    texto = (salida / "informe.txt").read_text(encoding="utf-8")
    assert "columna polar" in texto and "MAPA-tr" in texto
    res = json.loads((salida / "resultados.json").read_text())
    assert set(res) == {"ref", "polar"} and "700-820" in res["ref"] and "630-700" in res["polar"]
    assert texto.index("Banda 4") < texto.index("Banda 3")          # en el orden pedido, aunque acaben en otro


def test_interpolar_en_la_rejilla_de_las_tablas(datos):
    """Con --interpolar salen TABLA y TABLA-log; si la capa cae justo en un nodo (p, T) de la rejilla, coinciden con
    MAPA-trR."""
    sys.path.insert(0, os.path.join(AQUI, "herramientas"))
    import generar_tablas_k as GT
    ref = M.columna(4)
    col = M.columna(4, Ts=270.0)
    col["pm"] = np.array([GT.P_HPA[30], GT.P_HPA[25], GT.P_HPA[12], GT.P_HPA[5]])
    col["Tm"] = np.array([210.0, 240.0, 240.0, 270.0])
    res = M.medir_banda(datos, (700, 820), col, None, col_ref=ref, interpolar=True)
    for m in ("TABLA", "TABLA-log", "TABLA-1/T", "TABLA-15K"):
        assert np.allclose(res[m][0], res["MAPA-trR"][0], rtol=1e-9)
        assert np.allclose(res[m][1], res["MAPA-trR"][1], rtol=1e-9)

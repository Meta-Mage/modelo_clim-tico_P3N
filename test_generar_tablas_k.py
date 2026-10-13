# test_generar_tablas_k.py -- el generador de tablas k del infrarrojo con mapa espectral fijo (Fase 5.2, T9), con
# lineas sinteticas y un continuo MT_CKD sintetico.
import json
import os
import subprocess
import sys

import numpy as np
import pytest

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "herramientas"))

import fase52_infrarrojo as F       # noqa: E402
import fase52_mapa as MP            # noqa: E402
import fase52_referencia as R       # noqa: E402


def _lineas(rng, n, a, b, s0):
    return {"nu0": np.sort(rng.uniform(a, b, n)), "S": s0 * rng.lognormal(0, 2, n), "E": rng.uniform(0, 2000, n),
            "g_air": rng.uniform(0.05, 0.1, n), "g_self": rng.uniform(0.2, 0.5, n), "n_air": rng.uniform(0.6, 0.8, n),
            "d_air": rng.uniform(-0.005, 0.0, n), "iso": np.ones(n, dtype=int), "molec": np.ones(n, dtype=int)}


def _iso(masa):
    T = np.arange(60.0, 401.0)
    return {1: {"masa_u": masa, "T": T, "Q": 100.0 * (T / 296.0) ** 1.5}}


def _escribir_carpeta(carpeta, datos):
    from scipy.io import netcdf_file
    from test_fase52_hitran import _fortran
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
    # continuo sintetico, con las variables del archivo de MT_CKD (constante: 1e-24 propio y 1e-25 ajeno)
    with netcdf_file(str(carpeta / "absco-ref_wv-mt-ckd.nc"), "w") as f:
        w = np.arange(-20.0, 20020.0, 10.0)
        f.createDimension("wavenumbers", len(w))
        for nombre, val in (("wavenumbers", w), ("self_absco_ref", np.full(len(w), 1e-24)),
                            ("for_absco_ref", np.full(len(w), 1e-25)), ("self_texp", np.full(len(w), 4.0))):
            v = f.createVariable(nombre, "d", ("wavenumbers",))
            v[:] = val
        for nombre, val in (("ref_press", 1013.0), ("ref_temp", 296.0)):
            v = f.createVariable(nombre, "d", ())
            v[()] = val


@pytest.fixture(scope="module")
def datos():
    rng = np.random.default_rng(21)
    return {"H2O": (_lineas(rng, 120, 610, 720, 3e-21), _iso(18.0)),
            "CO2": (_lineas(rng, 120, 620, 710, 3e-20), _iso(44.0))}


def test_mapa_y_k_por_grupo_dan_el_linea_a_linea_en_la_columna_de_referencia(datos):
    """Con el k de cada grupo de la propia tabla (en las capas de la columna de referencia) se recupera la transmision
    media de cada grupo de esa columna exactamente (es la definicion del promedio), y los flujos quedan cerca del linea
    a linea (medido con estas lineas sinteticas y 8 capas: 0,04 W/m2 al espacio y 0,17 al suelo, en una banda de
    ~10 y ~29 W/m2)."""
    banda = (630, 700)
    mapa = MP.mapa_banda(datos, banda, None, n_capas=8)
    assert mapa["pesos"].sum() == pytest.approx(1.0, abs=1e-12)
    col = R.columna(8)
    taus = R.espesores(datos, mapa["nu"], col, None)
    tau_g = 0.0
    for gas, (lin, iso) in datos.items():
        kg = np.array([MP.k_grupos(lin, iso, mapa, gas, col["pm"][c], col["Tm"][c], col["vmr"][gas][c], 1e30, 0.0)[0]
                       for c in range(8)])
        u = col["vmr"][gas] * col["n_aire"]
        peso = mapa["nu"][None, :] ** 3 / np.expm1(F.C2 * mapa["nu"][None, :] / col["Tm"][:, None]) * mapa["dnu"]
        directo = F.media_transmision_por_grupo(taus[gas], peso, mapa["grupo"], len(mapa["pesos"]))
        assert np.allclose(kg * u[:, None], directo, rtol=1e-9)
        tau_g = tau_g + kg * u[:, None]
    fr = MP.fracciones_planck(mapa, col["Th"])
    assert np.allclose(fr.sum(axis=0), 1.0, rtol=1e-12)
    Bb = F.planck_integral(*banda, col["Th"])
    ok = mapa["pesos"] > 0
    U, D = F.dos_flujos(tau_g.T[ok], mapa["pesos"][ok], fr[ok] * Bb[None, :] / mapa["pesos"][ok][:, None],
                        (fr[ok] * Bb[None, :] / mapa["pesos"][ok][:, None])[:, -1])
    import medir_solapamiento as M
    U0, D0 = M.flujos_lbl(sum(taus.values()), mapa["nu"], mapa["dnu"], col)
    assert abs(U[0] - U0[0]) < 0.1 and abs(D[-1] - D0[-1]) < 0.3


def test_u_referencia_interpola_en_log(datos):
    mapa = MP.mapa_banda(datos, (630, 700), None, n_capas=8)
    pm, u = mapa["pm_ref"], mapa["u_ref"]["CO2"]
    assert MP.u_referencia(mapa, "CO2", pm[3]) == pytest.approx(u[3], rel=1e-12)
    medio = np.sqrt(pm[3] * pm[4])
    assert MP.u_referencia(mapa, "CO2", medio) == pytest.approx(np.sqrt(u[3] * u[4]), rel=1e-9)
    assert MP.u_referencia(mapa, "CO2", 1e-3) == pytest.approx(u[0], rel=1e-12)


def test_programa_prueba_y_completo(tmp_path, datos):
    carpeta = tmp_path / "hitran"
    _escribir_carpeta(carpeta, datos)
    salida = tmp_path / "tablas"
    prog = os.path.join(AQUI, "herramientas", "generar_tablas_k.py")
    r = subprocess.run([sys.executable, prog, "--hitran", str(carpeta), "--salida", str(salida), "--prueba",
                        "--hilos", "2"], capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr
    assert "PRUEBA:" in r.stdout and (salida / "mapas" / "mapa_b03.npz").exists()
    # completo, solo la banda 3 y el CO2 (y continuar: lo hecho en la prueba no se repite)
    r = subprocess.run([sys.executable, prog, "--hitran", str(carpeta), "--salida", str(salida), "--bandas", "3",
                        "--gases", "CO2", "--hilos", "2"], capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stderr
    t = np.load(salida / "tabla_k_CO2.npz")
    assert t["k"].shape == (1, 36, 13, 1, 16) and np.all(t["k"] >= 0) and t["k"].max() > 0
    m = np.load(salida / "mapas_infrarrojo.npz")
    assert np.allclose(m["fraccion_planck"].sum(axis=1), 1.0)
    assert np.allclose(m["pesos"].sum(axis=1), 1.0)


def test_flujos_con_las_tablas_frente_al_linea_a_linea(tmp_path_factory, datos):
    """De punta a punta: tablas del generador (CO2 y vapor, banda 3) -> TablasInfrarrojo -> flujos, frente al linea a linea en
    la columna de referencia y en la polar (otra temperatura: interpolacion en la rejilla). Y en un nodo de la rejilla,
    k es el de la tabla."""
    import medir_solapamiento as M
    carpeta = tmp_path_factory.mktemp("h") / "hitran"
    _escribir_carpeta(carpeta, datos)
    salida = carpeta.parent / "tablas"
    prog = os.path.join(AQUI, "herramientas", "generar_tablas_k.py")
    r = subprocess.run([sys.executable, prog, "--hitran", str(carpeta), "--salida", str(salida), "--bandas", "3",
                        "--gases", "CO2", "H2O", "--hilos", "2"], capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stderr
    t = np.load(salida / "tabla_k_CO2.npz")
    th = np.load(salida / "tabla_k_H2O.npz")
    assert th["k"].shape == (1, 36, 13, 7, 16)
    assert np.array_equal(th["k"][0, 35, :, 5], th["k"][0, 35, :, 2])      # arriba: copiadas de la de 1e-4
    for modo in ("log", "lineal"):
        TI = F.TablasInfrarrojo(str(salida), interpolacion=modo)
        k = TI.k("CO2", np.array([t["p_hPa"][20]]), np.array([t["T_K"][3]]))
        assert np.allclose(k[0, :, 0], t["k"][0, 20, 3, 0], rtol=1e-9, atol=1e-300)
    TI = F.TablasInfrarrojo(str(salida))
    for Ts in (288.0, 255.0):
        col = R.columna(12, Ts=Ts)
        u = {g: col["vmr"][g] * col["n_aire"] for g in ("CO2", "H2O")}
        U, D = TI.flujos(col["pm"], col["Tm"], col["Th"], col["Th"][-1], u, col["vmr"]["H2O"])
        nu, dnu = R.rejilla_comun(datos, (630, 700), col)
        tau = sum(R.espesores(datos, nu, col, None).values())
        U0, D0 = M.flujos_lbl(tau, nu, dnu, col)
        # mismo resultado salvo el mapa, el promedio y la interpolacion (medido: -0,036 W/m2 al espacio y -0,041 al
        # suelo con 288 K; -0,018 y -0,033 con 255 K; en ~10 y ~18-29 W/m2)
        assert abs(U[0] - U0[0]) < 0.3 and abs(D[-1] - D0[-1]) < 0.3, (Ts, U[0] - U0[0], D[-1] - D0[-1])

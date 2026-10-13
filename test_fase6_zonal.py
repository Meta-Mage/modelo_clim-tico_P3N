# test_fase6_zonal.py -- M3N 3.21.0: medias zonales por nivel en la climatologia (fase6_zonal.py), sin cambiar el
# modelo; climatologia continuada desde otra (clima_dinamico.py --continuar-desde); su exportacion.
import json
import os
import pickle
import subprocess
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
ORDEN = [sys.executable, os.path.join(AQUI, "clima_dinamico.py"), "--tierra", "--dias", "1"]


def correr(carpeta, *extra, zonal=True):
    env = dict(os.environ, M3N_MODO="tierra", M3N_ZONAL="1" if zonal else "0")
    return subprocess.run(ORDEN + ["--carpeta", str(carpeta)] + list(extra), capture_output=True, text=True,
                          cwd=AQUI, env=env)


def leer(carpeta):
    with open(os.path.join(carpeta, "resultado.pkl"), "rb") as f:
        return pickle.load(f)


def _bits(x):
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64)).view(np.uint64)


def test_funcion_corriente_igual_que_validar_i16():
    import validar_i16 as V
    from fase6_zonal import funcion_corriente
    rng = np.random.default_rng(3)
    vps = rng.standard_normal((20, 35))
    sig = np.linspace(0, 1, 21)
    cosc = np.cos(np.radians(np.linspace(85, -85, 35)))
    assert np.array_equal(funcion_corriente(vps, sig, cosc, 6.4e6, 9.8), V.funcion_corriente(vps, sig, cosc, 6.4e6, 9.8))


def test_zonal_no_cambia_el_modelo_y_se_guarda(tmp_path):
    con, sin = tmp_path / "con", tmp_path / "sin"
    r1, r2 = correr(con), correr(sin, zonal=False)
    assert r1.returncode == 0, r1.stdout + r1.stderr
    assert r2.returncode == 0, r2.stdout + r2.stderr
    a, b = leer(con), leer(sin)
    for k in ("reg_media", "reg_min", "reg_max", "suelo_media", "hielo_espesor", "viento_u_baja", "humedad_relativa_baja"):
        assert np.array_equal(_bits(a[k]), _bits(b[k])), k                       # identico bit a bit
    for k in ("precipitacion", "evaporacion"):
        assert np.array_equal(_bits(a["agua"][k]), _bits(b["agua"][k])), k
    assert b["zonal"] is None
    z = a["zonal"]
    assert z["anos"] == 2 and z["u"].shape[1:] == (20, 36) and z["vps"].shape[1:] == (20, 35) and z["T"].shape[1:] == (20, 36)
    pasos = len(a["reg_media"]) * 96                         # modo Tierra: 96 pasos de 900 s por dia; 2 subpasos
    assert z["subpasos_por_dia"].sum() == 2 * 2 * pasos      # dos años registrados
    assert np.isfinite(z["u"]).all() and 150 < z["T"].min() and z["T"].max() < 340
    assert "medias zonales por nivel (corte vertical de H3N): si, 2 años" in (con / "resumen.txt").read_text(encoding="utf-8")


def test_continuar_desde_otra_climatologia(tmp_path):
    base, nueva = tmp_path / "base", tmp_path / "nueva"
    assert correr(base).returncode == 0
    r = correr(nueva, "--continuar-desde", str(base))
    assert r.returncode == 0, r.stdout + r.stderr
    og = json.loads((nueva / "origen.json").read_text(encoding="utf-8"))
    assert og["carpeta"] == str(base) and og["anos_climatologia_origen"] == 2
    rn = leer(nueva)
    assert rn["origen"]["carpeta"] == str(base) and rn["zonal"]["anos"] == 2
    texto = (nueva / "resumen.txt").read_text(encoding="utf-8")
    assert "continua desde su estado final" in texto and "AVISO" not in texto.split("PROCEDENCIA")[1].split("BANDAS")[0]
    # 3.22.0: tampoco en la terminal (antes salia un falso "AVISO: esta simulacion empezo antes de la 3.14.0")
    assert "empezo antes de la 3.14.0" not in r.stdout and "son de " + str(base) in r.stdout
    # no se puede reutilizar una carpeta con simulacion propia, ni la misma carpeta
    r = correr(base, "--continuar-desde", str(nueva))
    assert r.returncode != 0 and "simulacion propia" in (r.stdout + r.stderr)


def test_exportacion_del_corte():
    from fase6_zonal import funcion_corriente
    from exportar_clima import exportar_zonal, leer_serie
    dias, N, F = 5, 20, 36
    sig = np.linspace(0, 1, N + 1)
    lat = 90 - (np.arange(F) + 0.5) * 5
    latc = 90 - np.arange(1, F) * 5.0
    rng = np.random.default_rng(0)
    z = {"u": rng.standard_normal((dias + 1, N, F)), "vps": rng.standard_normal((dias + 1, N, F - 1)) * 1e3,
         "T": 250 + rng.standard_normal((dias + 1, N, F)), "ps": np.full((dias + 1, F), 1e5),
         "sigma_semi": sig, "lat": lat, "lat_caras": latc, "cos_caras": np.cos(np.radians(latc)),
         "radio": 6.4e6, "gravedad": 9.8, "anos": 30}
    e = exportar_zonal(z, num_dias=dias)
    psi = leer_serie(e["funcion_corriente_1e10_kg_s"])
    assert psi.shape == (dias, N + 1, F - 1)
    esperado = funcion_corriente(z["vps"][2], sig, z["cos_caras"], 6.4e6, 9.8) / 1e10
    assert np.allclose(psi[2], esperado, atol=0.006)
    assert leer_serie(e["temperatura_C"]).shape == (dias, N, F)
    assert np.allclose(leer_serie(e["presion_superficie_kPa"]), 100.0)
    assert e["anos"] == 30 and len(e["sigma_capas"]) == N

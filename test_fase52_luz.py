# test_fase52_luz.py -- lectores y primeros calculos de los datos de la luz (Fase 5.2), con archivos sinteticos.
import numpy as np
import pytest

import fase52_hitran as H
import fase52_infrarrojo as F
import fase52_luz as L


def _fits(ruta, datos, bitpix=-32):
    """FITS minimo: cabecera de tarjetas de 80 caracteres en bloques de 2880 bytes, datos big-endian."""
    tarj = [f"{'SIMPLE':8s}= {'T':>20s}", f"{'BITPIX':8s}= {bitpix:>20d}", f"{'NAXIS':8s}= {1:>20d}",
            f"{'NAXIS1':8s}= {len(datos):>20d}", "END"]
    cab = "".join(t.ljust(80) for t in tarj)
    cab = cab.ljust(((len(cab) + 2879) // 2880) * 2880)
    d = np.asarray(datos, dtype=">f4" if bitpix == -32 else ">f8").tobytes()
    d = d + b"\0" * ((-len(d)) % 2880)
    with open(ruta, "wb") as f:
        f.write(cab.encode("ascii") + d)


def _planck_lam(lam_um, T):
    return lam_um ** -5 / np.expm1(L.C2_UM_K / (lam_um * T))


def test_fits_y_espectro_interpolado(tmp_path):
    lam_A = np.geomspace(500.0, 55000.0, 40000)
    _fits(tmp_path / "w.fits", lam_A, -64)
    for T in (5300, 5400):
        _fits(tmp_path / f"m{T}.fits", _planck_lam(lam_A * 1e-4, T) * 1e20)
    lam, fl = L.espectro_phoenix(tmp_path / "w.fits", {5300: tmp_path / "m5300.fits", 5400: tmp_path / "m5400.fits"},
                                 5340.0)
    assert np.allclose(lam, lam_A * 1e-4, rtol=1e-12)
    a, b = _planck_lam(lam, 5300) * 1e20, _planck_lam(lam, 5400) * 1e20
    assert np.allclose(fl, np.exp(0.6 * np.log(a) + 0.4 * np.log(b)), rtol=1e-6)
    with pytest.raises(ValueError):
        L.espectro_phoenix(tmp_path / "w.fits", {5300: tmp_path / "m5300.fits", 5400: tmp_path / "m5400.fits"}, 6000.0)


def test_fracciones_de_bandas_de_un_cuerpo_negro():
    """Con un espectro de cuerpo negro de 0,05 a 5,5 µm (como PHOENIX) y la cola completada, las fracciones de cada
    banda de la luz son las de Planck (< 1e-4) y todas las bandas (de 820 a 50 000 cm^-1) recogen casi toda la luz."""
    T = 5420.0
    lam = np.geomspace(0.05, 5.5, 200000)
    bandas = sorted(set(H.BANDAS_LUZ) - {(2600, 3250)}) + [(2600, 3250)]
    bandas = sorted(bandas)
    fr, tot = L.fracciones_bandas(lam, _planck_lam(lam, T), bandas, T)
    sT4 = F.SIGMA * T ** 4
    exactas = np.array([F.planck_integral(a, b, T) for a, b in bandas]) / sT4
    assert np.allclose(fr, exactas, atol=1e-4)
    assert 0.98 < tot < 1.0


def test_ozono_serdyuchenko(tmp_path):
    lam = np.arange(213.0, 1100.0, 0.5)
    filas = []
    for x in lam:
        filas.append(f"{x:.2f} " + " ".join(f"{(1e-18 + 1e-21 * k) * np.exp(-((x - 255) / 20) ** 2):.5e}"
                                             for k in range(11)))
    (tmp_path / "o3.dat").write_text("cabecera\nWavelength 293K ...\n" + "\n".join(filas))
    d = L.leer_serdyuchenko(tmp_path / "o3.dat")
    assert d[2].shape == (11, len(lam)) and d[1][0] == 293.0 and d[1][-1] == 193.0
    s288 = L.seccion_ozono(np.array([255.0]), 288.0, d)[0]
    assert s288 == pytest.approx(1e-18 + 0.5e-21, rel=1e-9)            # a medio camino de 293 y 283 K
    assert L.seccion_ozono(np.array([255.0]), 400.0, d)[0] == pytest.approx(1e-18, rel=1e-12)
    assert np.isnan(L.seccion_ozono(np.array([200.0]), 250.0, d)[0])


def test_cia(tmp_path):
    def bloque(par, a, b, n, T, f):
        nu = np.linspace(a, b, n)
        cab = f"{par:>20s}{a:10.3f}{b:10.3f}{n:7d}{T:7.1f}{f(nu).max():10.3E}{0.5:6.3f}{'sintetico':>27s}{1:3d}"
        return [cab] + [f"{x:10.3f} {f(x):10.3E}" for x in nu]
    lin = (bloque("O2-O2", 7500, 8500, 101, 200.0, lambda x: 1e-46 + 0 * x)
           + bloque("O2-O2", 7500, 8500, 101, 300.0, lambda x: 3e-46 + 0 * x)
           + bloque("O2-O2", 9000, 9500, 51, 296.0, lambda x: 5e-47 + 0 * x))
    (tmp_path / "o2.cia").write_text("\n".join(lin) + "\n")
    b = L.leer_cia(tmp_path / "o2.cia")
    assert len(b) == 3 and b[0]["par"] == "O2-O2" and b[1]["T"] == 300.0
    k = L.k_cia(b, np.array([8000.0, 9200.0, 8800.0]), 250.0)
    assert k[0] == pytest.approx(2e-46) and k[1] == pytest.approx(5e-47) and k[2] == 0.0
    solapado = b + [{"par": "O2-O2", "T": 296.0, "nu": np.array([8400.0, 8600.0]), "k": np.array([1e-46, 1e-46])}]
    with pytest.raises(ValueError):
        L.k_cia(solapado, np.array([8000.0]), 250.0)


def test_descargar_luz_con_archivos_locales(tmp_path, monkeypatch):
    """El programa de descarga (con direcciones file:// en lugar de las webs): baja, apunta la huella, no repite lo
    ya bajado y dice claramente lo que falla."""
    import json
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "herramientas"))
    import descargar_luz as D
    origen = tmp_path / "web"
    origen.mkdir()
    (origen / "a.dat").write_bytes(b"uno")
    (origen / "b.dat").write_bytes(b"dos")
    monkeypatch.setattr(D, "ARCHIVOS", [("a.dat", (origen / "a.dat").as_uri(), "ozono"),
                                        ("b.dat", (origen / "b.dat").as_uri(), "estrella")])
    monkeypatch.setattr(D.time, "sleep", lambda s: None)
    destino = tmp_path / "luz"
    monkeypatch.setattr(sys, "argv", ["x", "--carpeta", str(destino)])
    D.main()
    p = json.loads((destino / "procedencia.json").read_text())
    assert p["archivos"]["a.dat"]["sha256"] == D.sha256(str(origen / "a.dat")) and (destino / "b.dat").exists()
    monkeypatch.setattr(D, "ARCHIVOS", D.ARCHIVOS + [("c.dat", (origen / "no_existe.dat").as_uri(), "colisiones")])
    with pytest.raises(SystemExit) as e:
        D.main()
    assert "no_existe" in str(e.value)

# test_sensibilidad_oceano.py -- 4.0.0 (Fase 6.4, DISENO_FASE_6.4.md §4.2): las pruebas de sensibilidad del oceano
# de clima_dinamico.py (--oceano-difusion, --oceano-capa-mezcla).
import json
import os
import pickle
import subprocess
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
ORDEN = [sys.executable, os.path.join(AQUI, "clima_dinamico.py"), "--tierra"]


def correr(*extra, carpeta):
    return subprocess.run(ORDEN + list(extra) + ["--carpeta", str(carpeta)], capture_output=True, text=True,
                          cwd=AQUI, env=dict(os.environ, M3N_MODO="tierra"))


def test_sensibilidad_exige_climatologia_exploratoria(tmp_path):
    r = correr("--oceano-difusion", "2", carpeta=tmp_path / "a")
    assert r.returncode != 0 and "--anos 12" in (r.stdout + r.stderr)
    r = correr("--oceano-capa-mezcla", "0", "--anos", "12", carpeta=tmp_path / "b")
    assert r.returncode != 0 and "mayores que 0" in (r.stdout + r.stderr)


def test_sensibilidad_se_aplica_se_guarda_y_no_se_mezcla(tmp_path):
    base, sens = tmp_path / "base", tmp_path / "sens"
    assert correr("--dias", "1", carpeta=base).returncode == 0
    assert correr("--dias", "1", "--oceano-difusion", "2", carpeta=sens).returncode == 0
    with open(sens / "sensibilidad.json", encoding="utf-8") as f:
        assert json.load(f) == {"oceano_difusion": 2.0, "oceano_capa_mezcla": 1.0}
    with open(sens / "resultado.pkl", "rb") as f:
        rs = pickle.load(f)
    with open(base / "resultado.pkl", "rb") as f:
        rb = pickle.load(f)
    assert rs["sensibilidad"]["oceano_difusion"] == 2.0 and rb["sensibilidad"] == {"oceano_difusion": 1.0,
                                                                                  "oceano_capa_mezcla": 1.0}
    assert not np.array_equal(np.asarray(rs["reg_media"]), np.asarray(rb["reg_media"]))   # D_oc distinto: otro clima
    assert "SENSIBILIDAD del oceano: D_oc x2" in (sens / "resumen.txt").read_text(encoding="utf-8")
    assert not (base / "sensibilidad.json").exists()
    # no se puede seguir una simulacion con otro oceano, en ningun sentido
    r = correr("--dias", "1", carpeta=sens)
    assert r.returncode != 0 and "otra configuracion del oceano" in (r.stdout + r.stderr)
    r = correr("--dias", "1", "--oceano-capa-mezcla", "0.5", carpeta=base)
    assert r.returncode != 0 and "otra configuracion del oceano" in (r.stdout + r.stderr)

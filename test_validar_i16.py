# test_validar_i16.py -- v3.12.0: pruebas de la validacion del modo Tierra con I16 (validar_i16.py,
# referencias_tierra.py; DISENO_FASE_6.3.md §6.14).
#
# Uso:  python -m pytest test_validar_i16.py      (alrededor de dos minutos: la ultima prueba hace una
#                                                  simulacion de prueba completa en modo Tierra)

import hashlib
import math
import os
import pickle
import subprocess
import sys

import numpy as np

import referencias_tierra as R
import validar_i16 as V

AQUI = os.path.dirname(os.path.abspath(__file__))


def test_referencias_tablas_y_medias():
    """Las tablas tienen las 36 filas de M3N y dan las medias globales de sus fuentes."""
    assert len(R.CRU_T_ANUAL) == 36 and np.array(R.CRU_T_MES).shape == (12, 36)
    assert len(R.GPCP_P_ANUAL) == 36 and np.array(R.GPCP_P_MES).shape == (12, 36)
    # CRU 1961-1990: 14,0 C global, 14,6 norte, 13,4 sur (Jones et al. 1999)
    assert abs(V.media_filas(R.CRU_T_ANUAL) - 13.97) < 0.01
    norte = (V.LAT > 0).astype(float)
    assert abs(V.media_filas(R.CRU_T_ANUAL, norte) - 14.59) < 0.01
    assert abs(V.media_filas(R.CRU_T_ANUAL, 1 - norte) - 13.36) < 0.01
    # el año es la media de los meses (CRU, sin pesar) y la media pesada por dias (GPCP)
    assert np.abs(np.array(R.CRU_T_MES).mean(axis=0) - R.CRU_T_ANUAL).max() < 0.006
    dias = np.array(R.DIAS_MES)
    p = (np.array(R.GPCP_P_MES) * dias[:, None]).sum(axis=0) / dias.sum()
    assert np.abs(p - R.GPCP_P_ANUAL).max() < 0.0006
    assert abs(V.media_filas(R.GPCP_P_ANUAL) - 2.674) < 0.001
    # la ITCZ de GPCP esta al norte del ecuador (7,5 N en la media anual)
    assert V.LAT[int(np.argmax(R.GPCP_P_ANUAL))] == 7.5


def test_calendario_del_modo_tierra():
    """Dia 0 = solsticio de diciembre (~21 dic): el 1 de enero es el dia 11."""
    assert V.mes_del_dia(0) == 11 and V.mes_del_dia(10) == 11
    assert V.mes_del_dia(11) == 0 and V.mes_del_dia(11 + 30) == 0 and V.mes_del_dia(11 + 31) == 1
    meses = V.mes_de_cada_dia(365)
    assert all(V.mes_del_dia(d) == meses[d] for d in range(365))
    assert np.bincount(meses, minlength=12)[0] == 31 and np.bincount(meses, minlength=12)[6] == 31


def test_transporte_de_una_convergencia_simetrica():
    """Convergencia (lo que recibe cada celda) positiva en el ecuador y negativa en los polos, simetrica: el
    transporte va hacia el ecuador (negativo en el norte, positivo en el sur), es antisimetrico y lo que se quita antes de integrar es la media global que se le añadio."""
    from rejilla import COLUMNAS
    phi = np.radians(V.LAT)
    w = np.cos(phi)
    z = np.cos(2 * phi)
    z = z - float((z * w).sum() / w.sum())
    c = (z + 2.5)[:, None] * np.ones((1, COLUMNAS))
    bordes, t, resto = V.transporte(c)
    assert abs(resto - 2.5) < 1e-12
    escala = np.abs(t).max()
    assert np.allclose(t, -t[::-1], atol=1e-12 * escala)                     # antisimetrico
    assert (t[bordes > 0] < 0).all() and (t[bordes < 0] > 0).all()


def test_funcion_corriente_de_una_celula_cerrada():
    """v p_s hacia el norte arriba y hacia el sur abajo (masa neta 0): psi = 0 en el suelo y en el tope, y en
    medio vale 2 pi a cos(phi) / g * (v p_s) * delta sigma."""
    sig = np.linspace(0.0, 1.0, 11)
    vps = np.zeros((10, 3))
    vps[:5, 1], vps[5:, 1] = 1e5, -1e5
    cosc = np.array([0.9, 0.95, 0.9])
    a, g = 6.371e6, 9.80665
    psi = V.funcion_corriente(vps, sig, cosc, a, g)
    assert psi.shape == (11, 3)
    assert abs(psi[0, 1]) == 0 and abs(psi[-1, 1]) < 1e-6 * abs(psi[5, 1])
    assert math.isclose(psi[5, 1], 2 * math.pi * a * 0.95 / g * 1e5 * 0.5, rel_tol=1e-12)
    assert (psi[:, 0] == 0).all()
    assert V.cruce_hacia_el_polo(np.array([-1.0, 1.0, 3.0]), np.array([10.0, 5.0, 0.0]), 2, -1) == 7.5


def test_validacion_de_punta_a_punta(tmp_path):
    """Simulacion de prueba completa (modo Tierra, "años" de 1 dia) y validacion con el año extra: el informe
    sale entero, el punto de control de la simulacion larga NO cambia, y leer los vientos no cambia el año
    extra (mismo resultado, bit a bit, que el año simulado sin leerlos)."""
    carpeta, salida = str(tmp_path / "clima"), str(tmp_path / "val")
    entorno = {k: v for k, v in os.environ.items() if k != "M3N_MODO"}
    subprocess.run([sys.executable, os.path.join(AQUI, "clima_dinamico.py"), "--tierra", "--dias", "1",
                    "--carpeta", carpeta], cwd=AQUI, env=entorno, check=True, capture_output=True)
    f_sim = os.path.join(carpeta, "estado_i16.pkl")
    antes = hashlib.sha256(open(f_sim, "rb").read()).hexdigest()
    subprocess.run([sys.executable, os.path.join(AQUI, "validar_i16.py"), "--dias", "1", "--carpeta", carpeta,
                    "--salida", salida], cwd=AQUI, env=entorno, check=True, capture_output=True)
    assert hashlib.sha256(open(f_sim, "rb").read()).hexdigest() == antes
    texto = open(os.path.join(salida, "informe.txt"), encoding="utf-8").read()
    for parte in ("1. AIRE A 2 m", "2. PRECIPITACION", "3. TRANSPORTE", "4. BALANCE DE RADIACION",
                  "5. HIELO MARINO", "6. PERFIL VERTICAL", "7. VIENTOS"):
        assert parte in texto, parte
    assert not os.path.exists(os.path.join(salida, "estado_i16_copia.pkl"))
    # el mismo año extra, sin leer los vientos
    guion = f"""
import os, sys, pickle, shutil
os.environ["M3N_MODO"] = "tierra"
sys.path.insert(0, {AQUI!r})
import validar_i16 as V, fase2b_atmosfera as F
V.args.dias = 1
argumentos, I = V.argumentos_simulacion()
copia = {str(tmp_path / "copia.pkl")!r}
shutil.copyfile({f_sim!r}, copia)
ck = pickle.load(open(copia, "rb"))
ini = pickle.load(open(os.path.join({carpeta!r}, "arranque_dos_capas.pkl"), "rb"))["estado"]
r = F.simular_fase2b(*argumentos, max_anos=ck["ano"], interruptores=I, estado_inicial=ini, archivo_estado=copia,
                     guardar_al_terminar=True)
pickle.dump(r["flujos"], open({str(tmp_path / "sin_vientos.pkl")!r}, "wb"))
"""
    entorno2 = dict(entorno, M3N_MODO="tierra")
    subprocess.run([sys.executable, "-c", guion], cwd=AQUI, env=entorno2, check=True, capture_output=True)
    sin = pickle.load(open(tmp_path / "sin_vientos.pkl", "rb"))
    con = pickle.load(open(os.path.join(salida, "ano_extra.pkl"), "rb"))["r"]["flujos"]
    for k in ("toa_neto", "olr_celda", "T_atm", "sw_suelo"):
        assert np.array_equal(sin[k], con[k]), k

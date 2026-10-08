# test_validar_i16.py -- 3.12.0 (3.13.0: celulas de Hadley, tropopausa y ERA-Interim; 3.14.0: chorros y trazabilidad): pruebas de la validacion del modo Tierra con I16 (validar_i16.py,
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


def test_celulas_de_hadley_cada_una_en_su_lado():
    """3.13.0: la celula del sur no puede ser la de Ferrel del norte (mismo signo), ni la del norte la de Ferrel
    del sur, aunque sean mas intensas (el fallo de diciembre-febrero en la 3.12.0)."""
    lat_c = np.arange(85.0, -86.0, -5.0)                       # caras de norte a sur
    nivel = np.linspace(0.0, 1.0, 21)
    perfil = np.sin(np.pi * nivel)[:, None]
    def celula(a, centro, ancho):
        return a * np.exp(-((lat_c - centro) / ancho) ** 2)[None, :] * perfil
    psi = (celula(14.0, 5.0, 12.0) + celula(-2.0, -12.0, 8.0)         # celulas de Hadley de diciembre-febrero
           + celula(-6.0, 40.0, 6.0) + celula(5.0, -45.0, 6.0))       # Ferrel: norte (<0) y sur (>0), mas fuertes
    k5 = int(np.argmin(np.abs(nivel - 0.5)))
    c = V.celulas_hadley(psi, lat_c, nivel, k5)
    assert c["norte"][0] > 13 and abs(c["norte"][1] - 5.0) < 0.1
    assert c["sur"][0] < 0 and -40 < c["sur"][1] < 0               # la del sur, no la Ferrel del norte (+35)
    # con la busqueda de la 3.12.0 (las dos en 40 S-40 N), la "del sur" habria salido en +35
    trop = (np.abs(lat_c) < 40)[None, :] & ((nivel > 0.15) & (nivel < 0.95))[:, None]
    assert lat_c[np.unravel_index(np.nanargmin(np.where(trop, psi, np.nan)), psi.shape)[1]] == 35.0
    assert 5 < c["norte"][2] < 40 and -45 < c["sur"][2] < -12      # bordes hacia el polo de cada maximo


def test_tropopausa_interpolada_y_punto_mas_frio():
    """Perfil de 6,5 K/km hasta 11 km e isotermo encima, con niveles cada 2 km (10 y 12 km alrededor): la de la
    OMM en un nivel da 12 km; la interpolada queda entre 11 y 12 km; el punto mas frio es el isotermo."""
    zz = np.arange(0.0, 30001.0, 2000.0)
    Tz = 288.15 - 6.5 * np.minimum(zz, 11000.0) / 1000
    pz = 101325.0 * np.exp(-zz / 7000.0)
    d = V.perfil_banda(Tz, zz, pz)
    assert abs(d["tropopausa_km"] - 12.0) < 1e-9
    assert 11.0 < d["tropopausa_interp_km"] < 12.0
    assert abs(d["T_frio"] - (216.65 - 273.15)) < 1e-9 and abs(d["T_25km"] - (216.65 - 273.15)) < 1e-9


def test_referencias_era_interim_coherentes_con_otras_fuentes():
    """REF_ERAI frente a fuentes independientes: tropopausa del ecuador (Seidel et al. 2001: ~16,5 km, punto
    mas frio a unos -81 C), celula de invierno del norte (Dima y Wallace 2003: 14-16 en NCEP; Oort y Yienger
    1996: 20; unidades 10^10 kg/s), bordes hacia 30 grados y chorros de la troposfera hacia 30-35 grados."""
    E = R.REF_ERAI
    t = E["tropopausa"]["10..-10"]
    assert abs(t["tropopausa_km"] - 16.5) < 0.5 and abs(t["T_frio"] + 81) < 2
    assert 14 <= E["hadley"]["dic-feb"]["norte"][0] <= 22 and -26 <= E["hadley"]["jun-ago"]["sur"][0] <= -14
    a = E["hadley"]["anual"]
    assert a["norte"][0] > 0 > a["sur"][0] and 28 < a["norte"][2] < 36 and -36 < a["sur"][2] < -28
    for h, signo in (("norte", 1), ("sur", -1)):
        v, la, niv = E["chorros"]["troposfera_" + h]
        assert 20 < v < 40 and 25 < signo * la < 45 and 0.1 <= niv <= 0.5
    assert len(E["u_bandas"][250]) == 36 and R.REF_TRANSPORTE_ECUADOR["atm_PW"] < 0
    # 3.14.0: el chorro de invierno es mas fuerte que el de la media anual, y la estratosfera de verano tiene
    # viento del este (negativo)
    T = E["chorros_temporada"]
    assert T["dic-feb"]["troposfera_norte"][0] > E["chorros"]["troposfera_norte"][0] + 5
    assert T["jun-ago"]["troposfera_sur"][0] > E["chorros"]["troposfera_sur"][0] + 5
    assert T["dic-feb"]["estratosfera_sur"][0] < 0 and T["jun-ago"]["estratosfera_norte"][0] < 0


def test_chorro_de_la_troposfera_es_un_nucleo():
    """3.14.0: el chorro de la troposfera es un maximo local en latitud y altura. Un viento que crece con la altura
    hasta el tope (la parte de abajo de un chorro de la estratosfera) no es un chorro de la troposfera, aunque sea
    el mayor de la ventana: con la definicion de la 3.13.0 salia ahi, en el borde de arriba."""
    nivel = np.array([0.004, 0.01, 0.03, 0.08, 0.11, 0.15, 0.2, 0.26, 0.33, 0.41, 0.5, 0.6, 0.8, 0.98])
    lat = V.LAT
    phi = np.radians(lat)
    z = -np.log(nivel)                                                  # ~altura en escalas
    # chorro de la estratosfera: crece con la altura hasta el tope, centrado en 60 grados de cada hemisferio
    estrat = 20.0 * z[:, None] * np.exp(-((np.abs(lat)[None] - 60) / 10) ** 2)
    # chorro subtropical: nucleo a 0,2 (200 hPa) y 30 grados
    sub = 30.0 * np.exp(-((np.log(nivel / 0.2)) / 0.5) ** 2)[:, None] * np.exp(-((np.abs(lat)[None] - 30) / 8) ** 2)
    u = estrat + sub
    k011 = int(np.where(nivel == 0.11)[0][0])
    assert u[k011, np.argmin(np.abs(lat - 57.5))] > 30.0                # la cola de la estratosfera gana en la ventana
    c = V.chorros(u, lat, nivel)
    for h, signo in (("norte", 1), ("sur", -1)):
        v, la, niv = c[("troposfera", h)]
        assert abs(la - signo * 32.5) <= 2.5 and niv == 0.2 and 25 < v < 35, (h, c[("troposfera", h)])
        assert c[("estratosfera", h)][2] <= 0.0101
    # sin chorro subtropical: no hay nucleo en la troposfera
    c = V.chorros(estrat, lat, nivel)
    assert all(np.isnan(c[("troposfera", h)][0]) for h in ("norte", "sur"))
    assert V.texto_chorro(c[("troposfera", "norte")]).startswith("sin nucleo")


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
    for parte in ("A traves del ecuador", "Chorros, media anual", "Troposfera (norte)", "Estratosfera (sur)",
                  "Troposfera (norte, dic-feb)", "obs. 975 hPa", "interpolada", "dic-feb  norte  M3N"):   # 3.13-3.14
        assert parte in texto, parte
    # 3.14.0: trazabilidad: la simulacion apunta su version y su codigo; el año extra guarda el codigo
    from cache_simulacion import huella_codigo, MODULOS_I16
    import json
    tramos = json.load(open(os.path.join(carpeta, "procedencia.json"), encoding="utf-8"))
    import re
    version = re.search(r'^VERSION = "([^"]+)"', open(os.path.join(AQUI, "clima_dinamico.py"), encoding="utf-8").read(),
                        re.M).group(1)
    assert len(tramos) == 1 and tramos[0]["version"] == version and tramos[0]["ano_sim"] == 0
    assert tramos[0]["codigo"] == huella_codigo(MODULOS_I16)
    resumen = open(os.path.join(carpeta, "resumen.txt"), encoding="utf-8").read()
    assert "PROCEDENCIA" in resumen and huella_codigo(MODULOS_I16) in resumen and "AVISO" not in resumen
    assert huella_codigo(MODULOS_I16) in texto
    assert pickle.load(open(os.path.join(salida, "ano_extra.pkl"), "rb"))["codigo"] == huella_codigo(MODULOS_I16)
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

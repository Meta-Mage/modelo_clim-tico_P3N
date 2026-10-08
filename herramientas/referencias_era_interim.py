# herramientas/referencias_era_interim.py -- M3N 3.13.0 (3.14.0: chorros como nucleos y de cada temporada): calcula las referencias de vientos, celulas de Hadley
# y tropopausa (referencias_tierra.REF_ERAI) a partir del reanalisis ERA-Interim (Dee et al. 2011, QJRMS 137,
# 553-597), medias zonales mensuales de enero de 1979 a diciembre de 2016, 37 niveles de presion, rejilla de
# 1,5 grados, tal como las distribuye el paquete TropD (Adam et al. 2018, GMD 11, 4339) en su repositorio:
#     https://github.com/tropd/pytropd  (carpeta ValidationData; commit c25ce1118c748f678869b9e69d19456a73fab1f4)
# Comprueba la huella sha256 de cada archivo antes de usarlo.
#
# Mide la Tierra con LAS MISMAS funciones con las que validar_i16.py mide M3N, despues de pasar ERA-Interim a
# la rejilla de 5 grados de M3N (media por area en cada banda de 5 grados; la funcion de corriente,
# interpolada a las caras de M3N): asi la referencia tiene la misma resolucion horizontal que el modelo.
#
# Uso (no hace falta para validar; solo para rehacer o revisar REF_ERAI):
#     git clone --depth 1 https://github.com/tropd/pytropd.git /tmp/pytropd
#     python herramientas/referencias_era_interim.py /tmp/pytropd/ValidationData
# Escribe en pantalla el diccionario REF_ERAI que esta copiado en referencias_tierra.py.

import hashlib
import math
import os
import pprint
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import validar_i16 as V                       # noqa: E402  (importado: solo las funciones)
import validar_v30 as V30                     # noqa: E402
from rejilla import LATITUDES_GRADOS as LAT   # noqa: E402

HUELLAS = {"ua.nc": "308526acf9ff80deafc866073188d2a4f9eabcf3e32e06d8050788dbdb85a34a",
           "va.nc": "0af946e33b633f24ed4c53cdb2a6658430641c185152e10b430d89f7dd2f20b6",
           "ta.nc": "05385bccd007afc7b7650b9cb47a4e48f27f4b3311f958f130930562b82c6b7d",
           "zg.nc": "9889879b14c3b1b7d0c92a194659b9675b5f29ada7adc372661775d9e3cafda3",
           "uas.nc": "72cb459a7703e80ab5aafd6417bb590b0a9bd05182433c4393703d340e855b8c"}
RADIO = 6.371e6          # m (radio medio de la Tierra)
GRAVEDAD = 9.80665       # m/s2
TEMPORADAS = {"anual": tuple(range(1, 13)), "dic-feb": (12, 1, 2), "jun-ago": (6, 7, 8)}


def leer(carpeta, nombre):
    import netCDF4
    ruta = os.path.join(carpeta, nombre)
    with open(ruta, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    if h != HUELLAS[nombre]:
        sys.exit(f"{ruta}: la huella sha256 no coincide ({h}). No sigo.")
    with netCDF4.Dataset(ruta) as d:
        var = nombre[:-3]
        return (np.asarray(d[var][:], dtype=float), np.asarray(d["lat"][:], dtype=float),
                np.asarray(d["lev"][:], dtype=float) if "lev" in d.variables else None)


def media_temporada(x, clave):
    """x (..., meses) desde enero de 1979: media de los meses de la temporada en todos los años."""
    mes = np.arange(x.shape[-1]) % 12 + 1
    return x[..., np.isin(mes, TEMPORADAS[clave])].mean(axis=-1)


FINO = np.arange(-89.95, 90.0, 0.1)


def a_filas(z, lat):
    """z (..., lat) en la rejilla de ERA-Interim -> (..., 36) media por area en cada banda de 5 grados de M3N."""
    orden = np.argsort(lat)
    zf = np.apply_along_axis(lambda y: np.interp(FINO, lat[orden], y[orden]), -1, z)
    w = np.cos(np.radians(FINO))
    out = np.empty(z.shape[:-1] + (len(LAT),))
    for i, c in enumerate(LAT):
        m = np.abs(FINO - c) < 2.5
        out[..., i] = (zf[..., m] * w[m]).sum(-1) / w[m].sum()
    return out


def funcion_corriente_p(v, lat, lev):
    """psi (niveles, lat) en kg/s: 2 pi a cos(phi) / g * integral de v dp desde el tope (p = 0; de 0 a 1 hPa
    con el v de 1 hPa), por trapecios entre niveles."""
    p = lev * 100.0
    acum = np.empty_like(v)
    acum[0] = v[0] * p[0]
    for k in range(1, len(p)):
        acum[k] = acum[k - 1] + 0.5 * (v[k] + v[k - 1]) * (p[k] - p[k - 1])
    return 2 * math.pi * RADIO * np.cos(np.radians(lat))[None, :] / GRAVEDAD * acum


def main(carpeta):
    ua, lat, lev = leer(carpeta, "ua.nc")
    va, _, _ = leer(carpeta, "va.nc")
    ta, _, _ = leer(carpeta, "ta.nc")
    zg, _, _ = leer(carpeta, "zg.nc")
    uas, _, _ = leer(carpeta, "uas.nc")
    assert ua.shape == (37, 121, 456) and lev[0] == 1 and lev[-1] == 1000 and lat[0] == 90
    nivel = lev / 1000.0
    caras = 0.5 * (LAT[1:] + LAT[:-1])                     # caras de M3N, de norte a sur
    k5 = int(np.argmin(np.abs(nivel - 0.5)))
    R = {"fuente": "ERA-Interim (Dee et al. 2011), medias zonales mensuales 1979-2016 de TropD "
                   "(github.com/tropd/pytropd, ValidationData, commit c25ce11), en la rejilla de 5 grados de M3N",
         "hadley": {}, "chorros": {}, "chorros_temporada": {}, "vientos_bajos_975hPa": {}, "vientos_bajos_10m": {}, "u_bandas": {},
         "tropopausa": {}}
    for clave in TEMPORADAS:
        psi = funcion_corriente_p(media_temporada(va, clave), lat, lev) / 1e10
        orden = np.argsort(lat)
        psi_c = np.array([np.interp(caras, lat[orden], fila[orden]) for fila in psi])
        R["hadley"][clave] = {k: tuple(round(x, 2) for x in t) for k, t in
                              V.celulas_hadley(psi_c, caras, nivel, k5).items()}
    u = a_filas(media_temporada(ua, "anual"), lat)        # (37, 36)
    redondeo = lambda t: tuple(round(x, 4 if i == 2 else 2) for i, x in enumerate(t))
    R["chorros"] = {f"{c}_{h}": redondeo(t) for (c, h), t in V.chorros(u, LAT, nivel).items()}
    for clave in ("dic-feb", "jun-ago"):                   # 3.14.0: chorros de cada temporada
        ut = a_filas(media_temporada(ua, clave), lat)
        R["chorros_temporada"][clave] = {f"{c}_{h}": redondeo(t) for (c, h), t in V.chorros(ut, LAT, nivel).items()}
    k975 = int(np.where(lev == 975)[0][0])
    R["vientos_bajos_975hPa"] = {k: tuple(round(x, 2) for x in t) for k, t in V.vientos_bajos(u[k975], LAT).items()}
    R["vientos_bajos_10m"] = {k: tuple(round(x, 2) for x in t)
                              for k, t in V.vientos_bajos(a_filas(media_temporada(uas, "anual"), lat), LAT).items()}
    for p_ in (250, 500, 975):
        k = int(np.where(lev == p_)[0][0])
        R["u_bandas"][p_] = [round(float(x), 2) for x in u[k]]
    T = media_temporada(ta, "anual")
    Z = media_temporada(zg, "anual")
    for (n, s) in V30.BANDAS:
        m = (FINO < n) & (FINO > s)
        w = np.cos(np.radians(FINO[m]))
        orden = np.argsort(lat)
        Tb = np.array([(np.interp(FINO, lat[orden], f[orden])[m] * w).sum() / w.sum() for f in T])
        Zb = np.array([(np.interp(FINO, lat[orden], f[orden])[m] * w).sum() / w.sum() for f in Z])
        pb = lev * 100.0
        d = V.perfil_banda(Tb[::-1], Zb[::-1], pb[::-1])     # de abajo arriba
        R["tropopausa"][f"{n}..{s}"] = {k: round(float(x), 2) for k, x in d.items()}
    pprint.pprint(R, width=118, sort_dicts=False)
    return R


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Uso: python herramientas/referencias_era_interim.py RUTA/pytropd/ValidationData")
    main(sys.argv[1])

# descargar_luz.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7): descarga, UNA vez y en el PC de Carlos (desde el entorno
# de la IA no se llega a estas webs), los datos de la LUZ que no son lineas de HITRAN. No forma parte del modelo.
#
# Uso (con el venv de M3N activado):
#     python herramientas/descargar_luz.py                       # en ~/Documentos/B3N/datos/luz
#     python herramientas/descargar_luz.py --carpeta OTRA
# Si se corta, se vuelve a lanzar: lo que ya esta descargado no se repite.
#
# Que descarga (~50 MB):
#   1. EL ESPECTRO DE LA ESTRELLA: modelos de atmosfera estelar PHOENIX (Husser et al. 2013, A&A 553, A6; la web pide
#      citarlo), version ACES-AGSS-COND-2011, alta resolucion, de 0,05 a 5,5 µm: la rejilla de longitudes de onda y los
#      modelos de 5300, 5400 y 5500 K con log g = 4,5 y [Fe/H] = -0,5 (la estrella de referencia de S3N es tau Ceti,
#      G8 V, pobre en metales: FICHA_COMPOSICION_P3N.md §6.3; Teff de tau Ceti entre 5264 y 5525 K segun la fuente,
#      Tang y Gai 2011, A&A 526, A35, tabla 1; S3N, ~5420 K). De ftp.gwdg.de (el servidor de los datos de PHOENIX).
#   2. EL OZONO EN EL ULTRAVIOLETA Y EL VISIBLE: secciones eficaces de Serdyuchenko et al. 2014 (Atmos. Meas. Tech. 7,
#      625; IUP Bremen), 11 temperaturas de 193 a 293 K, de 213 a 1100 nm (longitud de onda en el vacio), version del
#      25/07/2012.
#   3. LA ABSORCION INDUCIDA POR COLISIONES de HITRAN2024 (Terragni et al. 2025, J. Quant. Spectrosc. Radiat. Transf.
#      347, 109631): N2-N2, O2-O2 y O2-N2 (carpeta "main" de hitran.org/cia).
# Escribe procedencia.json con la fecha, la direccion y la huella sha256 de cada archivo (para la trazabilidad).

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.request

PHOENIX = "https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/"
MODELOS = [f"lte0{t}-4.50-0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits" for t in (5300, 5400, 5500)]
ARCHIVOS = (
    [("WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", PHOENIX + "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits", "estrella")]
    + [(m, PHOENIX + "PHOENIX-ACES-AGSS-COND-2011/Z-0.5/" + m, "estrella") for m in MODELOS]
    + [("serdyuchenkogorshelev5digits.dat",
        "https://www.iup.uni-bremen.de/gruppen/molspec/downloads/serdyuchenkogorshelev5digits.dat", "ozono")]
    + [(f, "https://hitran.org/data/CIA/main/" + f, "colisiones")
       for f in ("N2-N2_2021.cia", "O2-O2_2024.cia", "O2-N2_2024.cia")]
)


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def bajar(url, destino, intentos=3):
    tmp = destino + ".tmp"
    for i in range(intentos):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "M3N-descargar_luz/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(tmp) == 0:
                raise IOError("archivo vacio")
            os.replace(tmp, destino)
            return
        except Exception as e:                       # se reintenta; al final, se dice cual fallo y por que
            if os.path.exists(tmp):
                os.remove(tmp)
            if i == intentos - 1:
                raise RuntimeError(f"no se pudo descargar {url}: {e}")
            time.sleep(5)


def main():
    ap = argparse.ArgumentParser(description="Descarga los datos de la luz para la Fase 5.2 de M3N")
    ap.add_argument("--carpeta", default=os.path.expanduser(os.path.join("~", "Documentos", "B3N", "datos", "luz")))
    a = ap.parse_args()
    os.makedirs(a.carpeta, exist_ok=True)
    ruta_proc = os.path.join(a.carpeta, "procedencia.json")
    proc = {"archivos": {}}
    if os.path.exists(ruta_proc):
        with open(ruta_proc, encoding="utf-8") as f:
            proc = json.load(f)
    fallos = []
    for nombre, url, que in ARCHIVOS:
        destino = os.path.join(a.carpeta, nombre)
        if os.path.exists(destino) and nombre in proc["archivos"]:
            print(f"{nombre}: ya descargado, no se repite", flush=True)
            continue
        print(f"{nombre} ({que}): descargando...", flush=True)
        try:
            bajar(url, destino)
        except RuntimeError as e:
            fallos.append(str(e))
            print(f"  FALLO: {e}", flush=True)
            continue
        proc["archivos"][nombre] = {"que": que, "url": url, "bytes": os.path.getsize(destino),
                                    "sha256": sha256(destino), "fecha": time.strftime("%Y-%m-%d %H:%M:%S")}
        print(f"  {os.path.getsize(destino) / 1e6:.1f} MB", flush=True)
        with open(ruta_proc, "w", encoding="utf-8") as f:
            json.dump(proc, f, indent=1, ensure_ascii=False)
    if fallos:
        sys.exit("No se pudo descargar todo (lo demas queda guardado; vuelve a lanzarlo):\n  " + "\n  ".join(fallos))
    print(f"\nListo en {a.carpeta}. Pasale a Claude el archivo procedencia.json.")


if __name__ == "__main__":
    main()

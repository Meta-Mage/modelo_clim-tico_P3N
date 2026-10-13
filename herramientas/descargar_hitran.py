# descargar_hitran.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7): descarga las LINEAS de HITRAN de los gases que
# importan para la radiacion de P3N y guarda las FUNCIONES DE PARTICION (TIPS) como tabla de datos (decision T3).
# Se ejecuta UNA vez en el PC de Carlos (desde el entorno de la IA no se llega a hitran.org). No forma parte del modelo.
#
# Uso (con el venv de M3N activado y "pip install hitran-api"):
#     python descargar_hitran.py                       # en ~/Documentos/B3N/datos/hitran
#     python descargar_hitran.py --carpeta OTRA        # en otra carpeta
# Si se corta, se vuelve a lanzar: lo que ya esta descargado no se repite.
#
# Que descarga: todos los isotopologos que conoce hitran-api 1.3.0.0 de
#     H2O (1), CO2 (2), O3 (3), N2O (4), CH4 (6), O2 (7)
# de 0 a 50 000 cm^-1 (infrarrojo y luz de la estrella), con los parametros del formato .par de 160 caracteres
# (numero de onda, intensidad, anchuras con el aire y consigo mismo, energia del nivel bajo, exponente de la
# temperatura, desplazamiento por la presion), que son los que usa fase52_lineas.py.
# Que escribe en la carpeta: <GAS>.data y <GAS>.header (formato de hitran-api), tips.json (Q(T) de cada isotopologo de
# 60 a 400 K, cada 1 K), isotopologos.json (masa y abundancia terrestre de cada uno), absco-ref_wv-mt-ckd.nc (el
# continuo del vapor MT_CKD 4.3 de AER, comprobado por su huella) y procedencia.json (fecha, version de hitran-api, rangos, numero de lineas y huella sha256 de
# cada archivo, para la trazabilidad).
#
# Datos: HITRAN2020 (Gordon et al. 2022, J. Quant. Spectrosc. Radiat. Transf. 277, 107949), de uso libre para la
# investigacion; hay que citarlo. hitran-api: Kochanov et al. 2016, JQSRT 177, 15 (licencia MIT).

import argparse
import hashlib
import json
import os
import sys
import time

GASES = {"H2O": 1, "CO2": 2, "O3": 3, "N2O": 4, "CH4": 6, "O2": 7}
NU_MIN, NU_MAX = 0.0, 50000.0           # cm^-1
T_TIPS = range(60, 401)                 # K
MTCKD_URL = ("https://raw.githubusercontent.com/AER-RC/MT_CKD/679f0a44129547d8c2baeea9cfd45ffb574356f1/"
             "data/absco-ref_wv-mt-ckd.nc")
MTCKD_SHA256 = "69944eb8b045c268e2daeb2cddf536b99f9e9efe7e91c425067b46a5b287ddb3"


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description="Descarga las lineas de HITRAN para la Fase 5.2 de M3N")
    ap.add_argument("--carpeta", default=os.path.expanduser(os.path.join("~", "Documentos", "B3N", "datos", "hitran")))
    args = ap.parse_args()
    try:
        import hapi
    except ImportError:
        sys.exit("Falta hitran-api: con el venv de M3N activado, 'pip install hitran-api' y vuelve a lanzarlo.")
    os.makedirs(args.carpeta, exist_ok=True)
    hapi.db_begin(args.carpeta)
    ids = {}
    for (m, i), v in hapi.ISO.items():
        ids.setdefault(m, []).append((i, v[hapi.ISO_INDEX["id"]]))
    proc = {"fecha": time.strftime("%Y-%m-%d %H:%M:%S"), "hitran_api": getattr(hapi, "HAPI_VERSION", "?"),
            "base_de_datos": "HITRAN (la vigente en hitran.org el dia de la descarga; HITRAN2020 en 10/2026)",
            "rango_cm-1": [NU_MIN, NU_MAX], "gases": {}}
    for nombre, m in GASES.items():
        datos = os.path.join(args.carpeta, nombre + ".data")
        if os.path.exists(datos) and os.path.getsize(datos) > 0:
            print(f"{nombre}: ya descargado, no se repite", flush=True)
            hapi.storage2cache(nombre)          # para contar sus lineas
        else:
            lista = sorted(g for _, g in ids[m])
            print(f"{nombre}: descargando {len(lista)} isotopologos de {NU_MIN:g} a {NU_MAX:g} cm^-1 "
                  "(puede tardar varios minutos y no muestra avance)...", flush=True)
            t0 = time.time()
            hapi.fetch_by_ids(nombre, lista, NU_MIN, NU_MAX)
            print(f"{nombre}: hecho en {time.time() - t0:.0f} s", flush=True)
        n = len(hapi.LOCAL_TABLE_CACHE[nombre]["data"]["nu"])
        proc["gases"][nombre] = {"molecula": m, "isotopologos": sorted(i for i, _ in ids[m]), "lineas": n,
                                 "sha256_data": sha256(datos),
                                 "sha256_header": sha256(os.path.join(args.carpeta, nombre + ".header"))}
        print(f"{nombre}: {n} lineas", flush=True)
    # funciones de particion (TIPS de hitran-api, sin red): tabla de datos para M3N (decision T3)
    tips = {"fuente": "TIPS de hitran-api (Gamache et al. 2017 y posteriores)", "T_K": list(T_TIPS), "Q": {}}
    for nombre, m in GASES.items():
        for i, _ in ids[m]:
            try:
                tips["Q"][f"{m}_{i}"] = [float(hapi.partitionSum(m, i, float(T))) for T in T_TIPS]
            except Exception as e:                   # algun isotopologo sin TIPS en esta version: se dice
                tips["Q"][f"{m}_{i}"] = None
                print(f"  aviso: sin funcion de particion para {nombre} isotopologo {i} ({e})", flush=True)
    with open(os.path.join(args.carpeta, "tips.json"), "w", encoding="utf-8") as f:
        json.dump(tips, f)
    proc["sha256_tips"] = sha256(os.path.join(args.carpeta, "tips.json"))
    # isotopologos: masa (para la anchura Doppler) y abundancia terrestre (la que llevan dentro las intensidades)
    iso = {f"{m}_{i}": {"id": v[hapi.ISO_INDEX["id"]], "nombre": v[hapi.ISO_INDEX["iso_name"]],
                        "abundancia": v[hapi.ISO_INDEX["abundance"]], "masa_u": v[hapi.ISO_INDEX["mass"]]}
           for (m, i), v in hapi.ISO.items() if m in GASES.values()}
    with open(os.path.join(args.carpeta, "isotopologos.json"), "w", encoding="utf-8") as f:
        json.dump(iso, f, indent=1, ensure_ascii=False)
    proc["sha256_isotopologos"] = sha256(os.path.join(args.carpeta, "isotopologos.json"))
    # el continuo del vapor MT_CKD 4.3 (AER), version fijada por su commit y comprobada por su huella
    destino = os.path.join(args.carpeta, "absco-ref_wv-mt-ckd.nc")
    if not (os.path.exists(destino) and sha256(destino) == MTCKD_SHA256):
        import urllib.request
        print("MT_CKD 4.3: descargando el continuo del vapor (83 kB)...", flush=True)
        urllib.request.urlretrieve(MTCKD_URL, destino + ".tmp")
        if sha256(destino + ".tmp") != MTCKD_SHA256:
            os.remove(destino + ".tmp")
            sys.exit("MT_CKD: la huella del archivo descargado no es la esperada; no se guarda.")
        os.replace(destino + ".tmp", destino)
    proc["mtckd"] = {"version": "4.3", "url": MTCKD_URL, "sha256": MTCKD_SHA256,
                     "licencia": "AER: uso cientifico y de investigacion; redistribuible con su aviso de copyright"}
    with open(os.path.join(args.carpeta, "procedencia.json"), "w", encoding="utf-8") as f:
        json.dump(proc, f, indent=1, ensure_ascii=False)
    total = sum(os.path.getsize(os.path.join(args.carpeta, x)) for x in os.listdir(args.carpeta))
    print(f"\nListo en {args.carpeta}: {sum(g['lineas'] for g in proc['gases'].values())} lineas, "
          f"{total / 1e6:.0f} MB. Pasale a Claude el archivo procedencia.json.")


if __name__ == "__main__":
    main()

# preparar_indices_refraccion.py -- Fase 5.2 de M3N (DISENO_FASE_5.2.md §7, decision tecnica T17): prepara el archivo
# datos/indices_refraccion_agua_hielo.npz con los INDICES DE REFRACCION del agua liquida y del hielo que usa
# fase52_indices.py (y, con el, la teoria de Mie de las gotas de las nubes). Se ejecuto UNA vez (13/10/2026) y el
# resultado esta en el repositorio; este programa queda para que cualquiera pueda repetirlo y comprobar que da lo mismo.
#
# Fuente de los datos: la base de datos de refractiveindex.info (Polyanskiy 2024, Sci. Data 11, 94), de DOMINIO
# PUBLICO (CC0 1.0), en github.com/polyanskiy/refractiveindex.info-database, commit fijado abajo. De ella, tres
# conjuntos, cada uno comprobado por su huella sha256:
#   - Hale y Querry 1973 (Appl. Opt. 12, 555): agua liquida a 25 °C, de 0,2 a 200 µm;
#   - Rowe, Fergoda y Neshyba 2020 (J. Geophys. Res. Atmos. 125, e2020JD032624): agua liquida a 240, 253, 263 y
#     273 K, de 1 a 15 000 cm^-1 (0,667 a 10 000 µm);
#   - Warren y Brandt 2008 (J. Geophys. Res. 113, D14220): hielo a 266 K, de 0,0443 µm a 2 m (comprobado frente a la
#     tabla impresa de los autores, atmos.washington.edu/ice_optical_constants, en 9 longitudes de onda: identico).
#
# Uso:
#     git clone https://github.com/polyanskiy/refractiveindex.info-database.git RUTA
#     git -C RUTA checkout c5c2f188e848453def5970e347399d653df2ffc2
#     python herramientas/preparar_indices_refraccion.py RUTA

import hashlib
import json
import os
import sys

import numpy as np

COMMIT = "c5c2f188e848453def5970e347399d653df2ffc2"
CARPETA = os.path.join("database", "data", "main", "H2O", "nk")
ARCHIVOS = {
    "Hale.yml": "df1af6b4352c3378cf81b149ac2280de30d80e2941d147786d5e7cd5044a8847",
    "Rowe-240K.yml": "0c1e808ecf0e30fa58ae0023882a2c396641f99b3e7d9c91d98f6c0e67e50663",
    "Rowe-253K.yml": "a8084501aaf81c47ae7badf696d99c252b8eef4b8c49148c8ee2f31f7273e76d",
    "Rowe-263K.yml": "0d7879587ed37c259ff1d6932a7892b6c797fb3490e49d2639f944a1286f157c",
    "Rowe-273K.yml": "c4bb53b5af850ae0488701087ff876701cec221e275a368db9d76051e06eb13f",
    "Warren-2008.yml": "a8577eaa1fabb9245957397344ce23efb20dd4c57cd7f442710d2d2373e0bbd6",
}
SALIDA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "datos",
                      "indices_refraccion_agua_hielo.npz")


def leer_tabla(ruta):
    """Las filas (longitud de onda en µm, n, k) del bloque "tabulated nk" de un archivo de refractiveindex.info.
    Lector propio y minimo (sin PyYAML): el formato es fijo, una fila de tres numeros por linea tras "data: |"."""
    filas, dentro = [], False
    with open(ruta, encoding="utf-8") as f:
        for linea in f:
            t = linea.strip()
            if t.startswith("- type:"):
                if "tabulated nk" not in t:
                    raise ValueError(f"{ruta}: se esperaba un bloque 'tabulated nk' y hay '{t}'")
                continue
            if t.startswith("data:"):
                dentro = True
                continue
            if dentro:
                partes = t.split()
                if len(partes) == 3:
                    try:
                        filas.append([float(x) for x in partes])
                        continue
                    except ValueError:
                        pass
                dentro = False
    a = np.array(filas, dtype=float)
    if a.ndim != 2 or a.shape[1] != 3 or len(a) < 10 or np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError(f"{ruta}: tabla vacia o longitudes de onda no crecientes")
    return a


def sha256(ruta):
    with open(ruta, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__ or "uso: python herramientas/preparar_indices_refraccion.py RUTA_DE_LA_BASE")
    base = os.path.join(sys.argv[1], CARPETA)
    t = {}
    for nombre, huella in ARCHIVOS.items():
        ruta = os.path.join(base, nombre)
        if sha256(ruta) != huella:
            sys.exit(f"{nombre}: la huella no es la esperada (¿otro commit de la base? se espera {COMMIT}).")
        t[nombre] = leer_tabla(ruta)
    rowe = [t[f"Rowe-{T}K.yml"] for T in (240, 253, 263, 273)]
    for r in rowe[1:]:
        if r.shape != rowe[0].shape or np.any(r[:, 0] != rowe[0][:, 0]):
            raise ValueError("las cuatro tablas de Rowe et al. no tienen la misma rejilla")
    procedencia = {
        "fuente": "refractiveindex.info-database (CC0 1.0), github.com/polyanskiy/refractiveindex.info-database",
        "commit": COMMIT, "sha256": ARCHIVOS,
        "agua_hale_querry_1973": "Appl. Opt. 12, 555, doi:10.1364/AO.12.000555; 25 °C (298 K)",
        "agua_rowe_2020": "J. Geophys. Res. Atmos. 125, e2020JD032624, doi:10.1029/2020JD032624; 240, 253, 263, 273 K",
        "hielo_warren_brandt_2008": "J. Geophys. Res. 113, D14220, doi:10.1029/2007JD009744; 266 K",
        "unidades": "longitud de onda en micrometros; m = n + i k",
    }
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "wb") as f:              # sin fecha dentro: repetirlo da el mismo archivo, byte a byte
        np.savez_compressed(
            f,
            hq_lam=t["Hale.yml"][:, 0], hq_n=t["Hale.yml"][:, 1], hq_k=t["Hale.yml"][:, 2],
            rowe_T=np.array([240.0, 253.0, 263.0, 273.0]), rowe_lam=rowe[0][:, 0],
            rowe_n=np.array([r[:, 1] for r in rowe]), rowe_k=np.array([r[:, 2] for r in rowe]),
            wb_lam=t["Warren-2008.yml"][:, 0], wb_n=t["Warren-2008.yml"][:, 1], wb_k=t["Warren-2008.yml"][:, 2],
            procedencia=np.array(json.dumps(procedencia, ensure_ascii=False, sort_keys=True)))
    print(f"Escrito {SALIDA} ({os.path.getsize(SALIDA) / 1e6:.2f} MB), sha256 {sha256(SALIDA)}")


if __name__ == "__main__":
    main()

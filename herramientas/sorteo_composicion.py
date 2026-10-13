# sorteo_composicion.py -- Fase 5.2 de M3N: el SORTEO de la composicion del aire de P3N (FICHA_COMPOSICION_P3N.md §6,
# aprobada por Carlos el 12/10/2026, 21:51: "Doy el visto bueno a todo").
#
# Reglas fijadas ANTES de sortear (ficha §6.1, §6.2 y §6.5):
#   - semilla: el entero de la huella sha256 (UTF-8) del texto "P3N · ficha de composición · 12/10/2026", modulo 2^63;
#     generador numpy.random.default_rng (PCG64);
#   - 1.er numero: la posicion s entre glaciacion (0) e interglaciar (1), uniforme en [0, 1); de ella salen los tres
#     gases a la vez, CO2 y CH4 en escala logaritmica y N2O en lineal, dentro de sus rangos naturales de los ultimos
#     800 000 años (CO2 172-300 ppm, Luthi et al. 2008; N2O 200-270 ppb, Schilt et al. 2010; CH4 350-800 ppb,
#     Loulergue et al. 2008 ⚠️ limites por comprobar: s no cambia al comprobarlos, solo el valor del CH4);
#   - 2.º numero: el oxigeno, uniforme en [18, 24] % (independiente).
# Se ejecuta tal cual; da siempre lo mismo. No forma parte del modelo.

import hashlib
import math

import numpy as np

TEXTO = "P3N · ficha de composición · 12/10/2026"
RANGOS = {"CO2_ppm": (172.0, 300.0, "log"), "CH4_ppb": (350.0, 800.0, "log"), "N2O_ppb": (200.0, 270.0, "lineal")}
O2_PCT = (18.0, 24.0)


def semilla():
    return int(hashlib.sha256(TEXTO.encode("utf-8")).hexdigest(), 16) % (2 ** 63)


def sortear():
    rng = np.random.default_rng(semilla())
    s = float(rng.uniform(0.0, 1.0))
    o2 = float(rng.uniform(*O2_PCT))
    gases = {}
    for k, (a, b, escala) in RANGOS.items():
        gases[k] = math.exp(math.log(a) + s * (math.log(b) - math.log(a))) if escala == "log" else a + s * (b - a)
    return {"semilla": semilla(), "s": s, **gases, "O2_pct": o2}


if __name__ == "__main__":
    r = sortear()
    print(f"texto de la semilla: {TEXTO!r}")
    print(f"semilla: {r['semilla']}")
    print(f"posicion glaciacion (0) - interglaciar (1): s = {r['s']:.6f}")
    print(f"CO2 {r['CO2_ppm']:.1f} ppm | CH4 {r['CH4_ppb']:.0f} ppb | N2O {r['N2O_ppb']:.1f} ppb | O2 {r['O2_pct']:.2f} %")

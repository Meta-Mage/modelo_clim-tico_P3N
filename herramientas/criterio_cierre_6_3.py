# herramientas/criterio_cierre_6_3.py -- aplica el criterio escrito de cierre de la Fase 6.3 (DISENO_FASE_6.3.md
# §6.22, aprobado por Carlos el 11/10/2026) a los textos de una simulacion larga del modo Tierra: el resumen.txt de
# clima_dinamico.py y el informe.txt de validar_i16.py. Solo LEE esos textos (no simula nada, no cambia nada).
#
# Uso (en la carpeta de M3N, con el venv activado):
#     python herramientas/criterio_cierre_6_3.py outputs/clima_dinamico/tierra_3.18.1/resumen.txt \
#            outputs/validacion_i16_3.18.1/informe.txt
#
# Salida: la tabla del criterio (A bloqueantes, B estructura con tolerancia de posicion de una celda = 5 grados,
# intensidades con su cociente frente a la referencia) y la lista de lo que necesita explicacion (regla de §6.22).
# Lo que no se puede comprobar con estos textos se marca "a mano" (por ejemplo, la exportacion y H3N, criterio A7).

import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fase6_clima import tendencia_santer          # noqa: E402  (la misma prueba de deriva de §6.17)

TOL_POSICION = 5.0          # grados: una celda de la rejilla de 5 (§6.22)
COCIENTE = (0.5, 2.0)       # intensidades: fuera de este intervalo, explicacion obligatoria (§6.22)
COTA_ENERGIA = 1e-7         # A4: cierre relativo de la energia
COTA_AGUA = 1e-11           # A4: kg/m2
UMBRAL_N = 0.2              # A5: W/m2, el del equilibrio (fase6_equilibrio)
NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"


def num(patron, texto, grupo=1):
    m = re.search(patron, texto, re.M)
    if not m:
        raise ValueError(f"no encuentro en el texto: {patron!r}")
    return float(m.group(grupo))


class Tabla:
    def __init__(self):
        self.filas, self.explicar = [], []

    def criterio(self, clave, texto, ok, detalle, explicable=True):
        estado = "a mano" if ok is None else ("✅" if ok else "❌")
        self.filas.append((clave, texto, estado, detalle))
        if ok is False and explicable:
            self.explicar.append(f"{clave} {texto}: {detalle}")

    def intensidad(self, clave, texto, m3n, obs, unidad=""):
        c = m3n / obs if obs else float("nan")
        fuera = not (COCIENTE[0] <= abs(c) <= COCIENTE[1]) or (m3n * obs < 0)
        signo = " (signo cambiado)" if m3n * obs < 0 else ""
        self.filas.append((clave, texto, "⚠️" if fuera else "·", f"M3N {m3n:g} / obs. {obs:g} {unidad} = {c:.2f}{signo}"))
        if fuera:
            self.explicar.append(f"{clave} {texto}: cociente {c:.2f}{signo}")


def seccion_a(t, res):
    # A1: un solo codigo en todos los tramos
    proc = res.split("PROCEDENCIA")[1].split("\n\n")[0]
    codigos = set(re.findall(r"codigo ([0-9a-f]{8,})", proc))
    t.criterio("A1", "cadena con un solo codigo", len(codigos) == 1, f"{len(codigos)} codigo(s): {', '.join(sorted(codigos))}",
               explicable=False)
    # A2: la cadena termino (el resumen existe y tiene la climatologia completa)
    m = re.search(r"CLIMATOLOGIA: (\d+) años registrados de (\d+)", res)
    t.criterio("A2", "estable, cadena completa", bool(m) and m.group(1) == m.group(2),
               f"{m.group(1)} de {m.group(2)} años registrados (recortes de vapor: en resultado.pkl, a mano)" if m else "sin climatologia",
               explicable=False)
    # A3: equilibrio
    conv = re.search(r"EQUILIBRIO: convergido (\w+) tras (\d+)", res)
    # 3.22.1: una climatologia continuada (clima_dinamico.py --continuar-desde) hereda el equilibrio de la de partida;
    # clima_dinamico.py solo deja continuar desde una climatologia que convergio y no tiene deriva
    heredado = re.search(r"EQUILIBRIO: el de (.+?) \((\d+) años de equilibrio", res)
    if conv:
        t.criterio("A3", "equilibrio", conv.group(1) == "True", f"convergido {conv.group(1)} tras {conv.group(2)} años",
                   explicable=False)
    elif heredado:
        t.criterio("A3", "equilibrio", True, f"heredado: continua desde el estado final de {heredado.group(1)} "
                   f"({heredado.group(2)} años de equilibrio; solo se puede continuar desde un equilibrio sin deriva)",
                   explicable=False)
    else:
        raise ValueError("el resumen no dice nada del equilibrio")
    # A4: cierres por año
    filas = re.findall(rf"^\s*\d+: cierre energia ({NUM}) \| error dinamica ({NUM}) W/m2 \| resto ({NUM}) \| agua atm ({NUM}) tierra ({NUM})",
                       res, re.M)
    v = np.array(filas, dtype=float)
    e_max, a_max = float(np.abs(v[:, 0]).max()), float(np.abs(v[:, 3:5]).max())
    t.criterio("A4", f"energia < {COTA_ENERGIA:g} y agua < {COTA_AGUA:g} kg/m2 cada año",
               e_max < COTA_ENERGIA and a_max < COTA_AGUA, f"max energia {e_max:.2e}, max agua {a_max:.2e} ({len(v)} años)",
               explicable=False)
    # A5: N del equilibrio y del año extra, y el error de la dinamica sin tendencia (Santer et al. 2000)
    b, sb, p, _ = tendencia_santer(v[:, 1])
    if heredado:                    # el N del equilibrio es el de la de partida (no el de esta)
        t.criterio("A5", "error de la dinamica controlado", p >= 0.05,
                   f"N del equilibrio: el de la de partida; error medio {v[:, 1].mean():.3f} W/m2, tendencia {b:+.4f} ± "
                   f"{sb:.4f} /año (p {p:.3f})", explicable=False)
    else:
        n_eq = num(rf"^\s*N\s+({NUM}) \(umbral", res)
        t.criterio("A5", "error de la dinamica controlado", abs(n_eq) < UMBRAL_N and p >= 0.05,
                   f"N del equilibrio {n_eq:+.3f} W/m2; error medio {v[:, 1].mean():.3f} W/m2, tendencia {b:+.4f} ± "
                   f"{sb:.4f} /año (p {p:.3f})", explicable=False)
    # A6: sin deriva
    sin = "sin deriva significativa" in res
    t.criterio("A6", "sin deriva significativa (Santer + BH)", sin,
               re.search(r"aire 2 m global: (.*?) K \|", res).group(1) + " K/año", explicable=False)


def seccion_a7(t, inf):
    partes = ("1. AIRE A 2 m", "2. PRECIPITACION", "3. TRANSPORTE", "4. BALANCE DE RADIACION", "5. HIELO MARINO",
              "6. PERFIL VERTICAL", "7. VIENTOS")
    faltan = [p for p in partes if p not in inf]
    t.criterio("A7", "validacion completa (exportacion y H3N: a mano)", not faltan,
               "informe con las 7 partes" if not faltan else f"faltan {faltan}", explicable=False)


def lat_de(txt):
    return float(txt)


def seccion_b(t, inf):
    # B1: vientos de superficie
    al_n = re.search(rf"Alisios norte.*?M3N\s+({NUM}) m/s a ({NUM})\s+obs\. 975 hPa\s+({NUM}) a ({NUM})", inf)
    al_s = re.search(rf"Alisios sur.*?M3N\s+({NUM}) m/s a ({NUM})\s+obs\. 975 hPa\s+({NUM}) a ({NUM})", inf)
    oe_n = re.search(rf"Oeste norte.*?M3N\s+({NUM}) m/s a ({NUM})\s+obs\. 975 hPa\s+({NUM}) a ({NUM})", inf)
    oe_s = re.search(rf"Oeste sur.*?M3N\s+({NUM}) m/s a ({NUM})\s+obs\. 975 hPa\s+({NUM}) a ({NUM})", inf)
    sentido = float(al_n.group(1)) < 0 and float(al_s.group(1)) < 0 and float(oe_n.group(1)) > 0 and float(oe_s.group(1)) > 0
    t.criterio("B1a", "alisios del este y oestes en los dos hemisferios", sentido,
               f"alisios {al_n.group(1)} / {al_s.group(1)} m/s; oestes {oe_n.group(1)} / {oe_s.group(1)} m/s")
    for nombre, m in (("norte", oe_n), ("sur", oe_s)):
        d = abs(float(m.group(2)) - float(m.group(4)))
        t.criterio("B1b", f"maximo del oeste de superficie ({nombre}) ±5°", d <= TOL_POSICION,
                   f"M3N {m.group(2)}° frente a {m.group(4)}° ({d:.1f}°)")
    # B2: celulas de Hadley
    had = {}
    for per, lado, val, la, bo, oval, ola, obo in re.findall(
            rf"^\s*(anual|dic-feb|jun-ago)\s+(norte|sur)\s+M3N\s+({NUM}) a\s+({NUM}), borde\s+({NUM})\s+obs\.\s+({NUM}) a\s+({NUM}), borde\s+({NUM})",
            inf, re.M):
        had[(per, lado)] = tuple(float(x) for x in (val, la, bo, oval, ola, obo))
    signos = all((v[0] > 0) == (lado == "norte") for (per, lado), v in had.items())
    t.criterio("B2a", "Hadley con el signo correcto (6 celdas)", signos and len(had) == 6, f"{len(had)} celdas leidas")
    inv_n = abs(had[("dic-feb", "norte")][0]) > abs(had[("jun-ago", "norte")][0])
    inv_s = abs(had[("jun-ago", "sur")][0]) > abs(had[("dic-feb", "sur")][0])
    t.criterio("B2b", "la de invierno mas fuerte que la de verano (los dos hemisferios)", inv_n and inv_s,
               f"norte {had[('dic-feb', 'norte')][0]:+.2f} frente a {had[('jun-ago', 'norte')][0]:+.2f}; "
               f"sur {had[('jun-ago', 'sur')][0]:+.2f} frente a {had[('dic-feb', 'sur')][0]:+.2f}")
    for (per, lado), v in had.items():
        d = abs(v[2] - v[5])
        t.criterio("B2c", f"borde de Hadley {per} {lado} ±5°", d <= TOL_POSICION, f"M3N {v[2]:+.1f}° frente a {v[5]:+.1f}° ({d:.1f}°)")
    # B3: chorro subtropical de invierno
    for etiqueta, nombre in (("Troposfera \\(norte, dic-feb\\)", "norte, dic-feb"), ("Troposfera \\(sur, jun-ago\\)", "sur, jun-ago")):
        m = re.search(rf"{etiqueta}\s+M3N\s+({NUM}) m/s a ({NUM}),.*?obs\.\s+({NUM}) m/s a ({NUM}),", inf)
        if m:
            d = abs(float(m.group(2)) - float(m.group(4)))
            t.criterio("B3", f"chorro subtropical ({nombre}) con nucleo ±5°", d <= TOL_POSICION,
                       f"M3N {m.group(1)} m/s a {m.group(2)}° frente a {m.group(3)} a {m.group(4)}° ({d:.1f}°)")
        else:
            t.criterio("B3", f"chorro subtropical ({nombre}) con nucleo", False, "sin nucleo")
    # B4: transporte de la atmosfera
    m = re.search(rf"Atmosfera, maximo\s+({NUM}) a ({NUM})\s+({NUM}) a ({NUM})\s+obs\. ({NUM}) a (\d+) N", inf)
    hacia_polos = float(m.group(1)) > 0 and float(m.group(3)) > 0
    t.criterio("B4a", "transporte hacia los polos en los dos hemisferios", hacia_polos, f"norte {m.group(1)} PW, sur {m.group(3)} PW")
    d_n = abs(float(m.group(2)) - float(m.group(6)))
    d_s = abs(abs(float(m.group(4))) - 40.0)
    t.criterio("B4b", "maximo de la atmosfera (norte) ±5° (obs. 43° N)", d_n <= TOL_POSICION, f"M3N {m.group(2)}° ({d_n:.1f}°)")
    t.criterio("B4b", "maximo de la atmosfera (sur) ±5° (obs. ~40° S)", d_s <= TOL_POSICION, f"M3N {m.group(4)}° ({d_s:.1f}°)")
    # B5: precipitacion
    m = re.search(rf"^\s*año\s+M3N\s+({NUM}) \(.*?obs\.\s+({NUM})", inf, re.M)
    d = abs(float(m.group(1)) - float(m.group(2)))
    t.criterio("B5a", "banda de lluvias (media anual) ±5°", d <= TOL_POSICION, f"M3N {m.group(1)}° frente a {m.group(2)}° ({d:.1f}°)")
    bandas = re.findall(rf"^\s*({NUM}) \|\s+({NUM}) \|\s+({NUM}) \|\s+({NUM}) \|\s+({NUM})\s*$",
                        inf.split("2. PRECIPITACION")[1].split("3. TRANSPORTE")[0], re.M)
    lat = np.array([float(b[0]) for b in bandas]); p = np.array([float(b[1]) for b in bandas])

    def hay_extremo(lo, hi, maximo):
        """¿Hay en cada hemisferio un extremo LOCAL (mayor o menor que sus dos vecinas, mirando el perfil entero)
        con |latitud| entre lo y hi? Devuelve {1: bool, -1: bool} y las latitudes encontradas."""
        res, donde = {}, {}
        for h in (1, -1):
            idx = sorted([i for i in range(len(lat)) if np.sign(lat[i]) == h], key=lambda i: abs(lat[i]))
            v = p[idx]; la = np.abs(lat[idx])
            ext = [float(la[k]) for k in range(1, len(v) - 1) if lo <= la[k] <= hi and
                   ((v[k] > v[k - 1] and v[k] > v[k + 1]) if maximo else (v[k] < v[k - 1] and v[k] < v[k + 1]))]
            res[h], donde[h] = bool(ext), ext
        return res, donde
    (mn, dmn), (mx, dmx) = hay_extremo(10, 45, False), hay_extremo(30, 70, True)
    t.criterio("B5b", "minimos subtropicales de lluvia en los dos hemisferios", mn[1] and mn[-1],
               f"norte en {dmn[1]}, sur en {dmn[-1]} (|lat|)")
    t.criterio("B5c", "maximos de latitudes medias en los dos hemisferios", mx[1] and mx[-1],
               f"norte en {dmx[1]}, sur en {dmx[-1]} (|lat|)")
    # B6: temperatura (cota rigurosa con la media anual y la amplitud del ciclo) y hielo
    tb = re.findall(rf"^\s*({NUM}) \|\s+({NUM}) \|\s+({NUM}) \|\s+({NUM}) \|\s+({NUM})\s+({NUM})\s*$",
                    inf.split("1. AIRE A 2 m")[1].split("2. PRECIPITACION")[0], re.M)
    T = {float(r[0]): (float(r[1]), float(r[4])) for r in tb}
    ecu_min = min(T[l][0] - T[l][1] for l in T if abs(l) < 10)            # el mes mas frio del ecuador, cota inferior
    pol_max = max(T[l][0] + T[l][1] for l in T if abs(l) > 80)            # el mes mas calido de los polos, cota superior
    t.criterio("B6a", "polos mas frios que el ecuador en todas las estaciones", pol_max < ecu_min,
               f"cota del polo mas calido {pol_max:.1f} °C < cota del ecuador mas frio {ecu_min:.1f} °C")
    hn = num(rf"Norte, maximo\s+({NUM})", inf); hs = num(rf"Sur, maximo\s+({NUM})", inf)
    t.criterio("B6b", "hielo marino en el invierno de los dos hemisferios", hn > 0 and hs > 0, f"norte {hn}, sur {hs} millones de km2")
    # B7: tropopausa
    filas = re.findall(rf"^\s*(\S+)\s+M3N\s+\|\s+{NUM} km\s+\|\s+({NUM}) km", inf, re.M)
    tro = {b: float(z) for b, z in filas}
    tropicos = min(tro["30..10"], tro["10..-10"], tro["-10..-30"]); polos = max(tro["90..60"], tro["-60..-90"])
    t.criterio("B7", "tropopausa mas alta en los tropicos que en los polos", tropicos > polos,
               f"tropicos >= {tropicos} km, polos <= {polos} km (interpolada)")
    return had


def intensidades(t, inf, had):
    for (per, lado), v in had.items():
        t.intensidad("I", f"Hadley {per} {lado}", v[0], v[3], "1e10 kg/s")
    for etiqueta in ("Troposfera \\(norte, dic-feb\\)", "Troposfera \\(sur, jun-ago\\)", "Troposfera \\(norte\\)",
                     "Troposfera \\(sur\\)", "Estratosfera \\(norte\\)", "Estratosfera \\(sur\\)"):
        m = re.search(rf"{etiqueta}\s+M3N\s+({NUM}) m/s.*?obs\.\s+({NUM}) m/s", inf)
        if m:
            t.intensidad("I", etiqueta.replace("\\", ""), float(m.group(1)), float(m.group(2)), "m/s")
    m = re.search(rf"Atmosfera, maximo\s+({NUM}) a {NUM}\s+({NUM}) a {NUM}\s+obs\. ({NUM})", inf)
    t.intensidad("I", "transporte de la atmosfera (norte)", float(m.group(1)), float(m.group(3)), "PW")
    t.intensidad("I", "transporte de la atmosfera (sur)", float(m.group(2)), float(m.group(3)), "PW")
    t.intensidad("I", "precipitacion global", num(rf"Media global\s+({NUM})\s+({NUM}) \(GPCP", inf),
                 num(rf"Media global\s+{NUM}\s+({NUM}) \(GPCP", inf), "mm/dia")
    t.intensidad("I", "agua precipitable", num(rf"Agua precipitable \(kg/m2\)\s+({NUM})", inf),
                 num(rf"Agua precipitable \(kg/m2\)\s+{NUM}\s+({NUM})", inf), "kg/m2")
    for nombre in ("Norte, maximo", "Norte, minimo", "Sur, maximo", "Sur, minimo"):
        m = re.search(rf"{nombre}\s+({NUM})\s+({NUM})", inf)
        t.intensidad("I", f"hielo marino {nombre.lower()}", float(m.group(1)), float(m.group(2)), "millones de km2")
    for nombre in ("Alisios norte", "Alisios sur", "Oeste norte", "Oeste sur"):
        m = re.search(rf"{nombre}.*?M3N\s+({NUM}) m/s.*?obs\. 975 hPa\s+({NUM})", inf)
        t.intensidad("I", f"viento de superficie: {nombre.lower()}", float(m.group(1)), float(m.group(2)), "m/s (975 hPa)")
    dT = num(rf"Media global\s+{NUM}\s+{NUM}\s+({NUM})", inf)
    t.filas.append(("I", "aire a 2 m global (diferencia, no cociente)", "·", f"M3N - obs. = {dT:+.2f} K (sin nubes: §6.14)"))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__ if __doc__ else "uso: criterio_cierre_6_3.py resumen.txt informe.txt")
    res = open(sys.argv[1], encoding="utf-8").read()
    inf = open(sys.argv[2], encoding="utf-8").read()
    t = Tabla()
    seccion_a(t, res)
    seccion_a7(t, inf)
    had = seccion_b(t, inf)
    intensidades(t, inf, had)
    print("CRITERIO DE CIERRE DE LA FASE 6.3 (DISENO_FASE_6.3.md §6.22)")
    print(f"resumen: {sys.argv[1]}\ninforme: {sys.argv[2]}\n")
    for clave, texto, estado, detalle in t.filas:
        print(f"{clave:4s} {estado:6s} {texto}\n            {detalle}")
    bloq = [f for f in t.filas if f[0].startswith("A") and f[2] == "❌"]
    print("\nBLOQUEANTES (A):", "todos cumplen" if not bloq else f"{len(bloq)} NO cumplen -> la 6.3 no se cierra")
    print(f"NECESITAN EXPLICACION CON FUENTE (B que falla o intensidad fuera de {COCIENTE[0]}-{COCIENTE[1]}): {len(t.explicar)}")
    for e in t.explicar:
        print("  -", e)


if __name__ == "__main__":
    main()

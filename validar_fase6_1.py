# validar_fase6_1.py -- validacion de la etapa 6.1 (aguas someras) con los casos de Williamson et al. (1992).
#
# Uso:  python3 validar_fase6_1.py [--rapido]
# Imprime una tabla por caso: errores normalizados l1, l2, l_inf (Williamson ✅), y la variacion
# relativa de la masa, la energia y la enstrofia potencial. --rapido omite el estudio de convergencia.

import math
import sys
import time
import numpy as np

from fase6_aguas_someras import (AguasSomeras, A_TIERRA, OMEGA_TIERRA, G_TIERRA, DIA, caso1, caso1_exacta,
                                 caso2, caso5, caso6, errores_normalizados, velocidad_fase_rh)

LAT_FILTRO = 60.0
DT_5GRADOS = 150.0            # s; con RAW (alfa 0,53) hace falta omega_max*dt <~ 0,45 (DISENO_FASE6_1.md §5)


def _modelo(filas=36, **kw):
    return AguasSomeras(A_TIERRA, OMEGA_TIERRA, G_TIERRA, filas=filas, lat_filtro=kw.pop("lat_filtro", LAT_FILTRO), **kw)


def _cambios(m, I0, estado):
    I = m.integrales(*estado)
    return {k: I[k] / I0[k] - 1 for k in I0}


def _fila(nombre, err, cambios, extra=""):
    print(f"  {nombre:<28} l1 {err[0]:.2e}  l2 {err[1]:.2e}  linf {err[2]:.2e}  | dM {cambios['masa']:+.1e}"
          f"  dE {cambios['energia']:+.1e}  dZ {cambios['enstrofia_potencial']:+.1e} {extra}")


def prueba_caso1(filas=36, dias=12):
    print(f"\nCaso 1: adveccion de la campana coseno, {dias} dias (una vuelta)")
    dt = DT_5GRADOS * 36 / filas * 2          # sin ondas de gravedad: el limite es el viento (~40 m/s)
    for nombre, alfa in (("alfa = 0", 0.0), ("alfa = pi/2 - 0,05", math.pi / 2 - 0.05), ("alfa = pi/2 (polos)", math.pi / 2)):
        h, u, v = caso1(alfa, filas=filas)
        m = _modelo(filas, viento_fijo=True)
        I0 = m.integrales(h, u, v)
        fin = m.integrar(h, u, v, dt, int(round(dias * DIA / dt)))
        ex = caso1_exacta(dias * DIA, alfa, filas=filas)
        c = _cambios(m, I0, fin)
        _fila(nombre, errores_normalizados(fin[0], ex, m.R.area), c, f" h_min {fin[0].min():.0f} m")


def prueba_caso2(filas=36, dias=5, imprimir=True):
    if imprimir:
        print(f"\nCaso 2: flujo zonal estacionario (solucion exacta), {dias} dias")
    dt = DT_5GRADOS * 36 / filas
    res = {}
    for nombre, alfa in (("alfa = 0", 0.0), ("alfa = 0,05", 0.05), ("alfa = pi/4", math.pi / 4),
                         ("alfa = pi/2 - 0,05", math.pi / 2 - 0.05), ("alfa = pi/2 (polos)", math.pi / 2)):
        h, u, v = caso2(alfa, filas=filas)
        m = _modelo(filas, alfa_coriolis=alfa)
        I0 = m.integrales(h, u, v)
        fin = m.integrar(h, u, v, dt, int(round(dias * DIA / dt)))
        err = errores_normalizados(fin[0], h, m.R.area)
        res[nombre] = err
        if imprimir:
            _fila(nombre, err, _cambios(m, I0, fin))
    return res


def prueba_caso5(filas=36, dias=15):
    print(f"\nCaso 5: flujo zonal sobre una montana aislada, {dias} dias (sin solucion exacta)")
    dt = DT_5GRADOS * 36 / filas
    h, u, v, hs = caso5(filas=filas)
    m = _modelo(filas, hs=hs)
    I0 = m.integrales(h, u, v)
    fin = m.integrar(h, u, v, dt, int(round(dias * DIA / dt)))
    c = _cambios(m, I0, fin)
    hh = fin[0] + hs
    print(f"  dM {c['masa']:+.1e}  dE {c['energia']:+.1e}  dZ {c['enstrofia_potencial']:+.1e}"
          f" | superficie libre h+hs: min {hh.min():.0f} m, max {hh.max():.0f} m")
    return fin, m


def prueba_caso6(filas=36, dias=14):
    print(f"\nCaso 6: onda de Rossby-Haurwitz de numero 4, {dias} dias")
    dt = DT_5GRADOS * 36 / filas
    h, u, v = caso6(filas=filas)
    m = _modelo(filas)
    I0 = m.integrales(h, u, v)
    fila50 = int(np.argmin(abs(m.R.phi_c - math.radians(50))))
    fases = []

    def reg(n, hh, uu, vv):
        c4 = np.fft.rfft(hh[fila50])[4]
        fases.append((n * dt, -np.angle(c4) / 4))
    fin = m.integrar(h, u, v, dt, int(round(dias * DIA / dt)), cada=int(DIA / dt), al_registrar=reg)
    t = np.array([x[0] for x in fases]); f = np.unwrap(np.array([x[1] for x in fases]) * 4) / 4
    nu = np.polyfit(t, f, 1)[0]
    c = _cambios(m, I0, fin)
    print(f"  dM {c['masa']:+.1e}  dE {c['energia']:+.1e}  dZ {c['enstrofia_potencial']:+.1e}")
    print(f"  velocidad angular de la onda: {math.degrees(nu) * DIA:.2f} grados/dia (formula de la vorticidad"
          f" barotropica no divergente: {math.degrees(velocidad_fase_rh()) * DIA:.2f}; en aguas someras algo menor)")
    print(f"  amplitud de la onda 4 en 50N: inicial {abs(np.fft.rfft(h[fila50])[4]) * 2 / h.shape[1]:.0f} m,"
          f" final {abs(np.fft.rfft(fin[0][fila50])[4]) * 2 / h.shape[1]:.0f} m")
    return fin, m


def convergencia():
    print("\nConvergencia del caso 2 (alfa = pi/4) al refinar la rejilla (orden 2 esperado):")
    prev = None
    for filas in (36, 72, 144):
        e = prueba_caso2(filas, imprimir=False)["alfa = pi/4"]
        orden = "" if prev is None else f"  orden l2 {math.log2(prev / e[1]):.2f}"
        print(f"  {180 / filas:.2f} grados: l2 {e[1]:.2e}{orden}")
        prev = e[1]


if __name__ == "__main__":
    t0 = time.time()
    prueba_caso1()
    prueba_caso2()
    prueba_caso5()
    prueba_caso6()
    if "--rapido" not in sys.argv:
        convergencia()
    print(f"\n(tiempo total {time.time() - t0:.0f} s)")

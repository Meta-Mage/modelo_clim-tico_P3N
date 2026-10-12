# fase6_zonal.py -- 3.21.0 (Fase 6.4; AUDITORIA_H3N_EXPLORAR P2.1, decision de Carlos del 12/10/2026, camino "a"):
# MEDIAS ZONALES POR NIVEL durante la climatologia con el nucleo dinamico (I16), para ver en H3N las celulas de la
# circulacion, los chorros y la temperatura en un corte latitud x altura.
#
# Que guarda: para cada dia del año, la media zonal (y entre todos los subpasos del nucleo de ese dia) de
#   u   (capas x filas)            viento zonal en los centros de celda, m/s
#   vps (capas x caras interiores) v * p_s en las caras entre filas, para la funcion de corriente
#   T   (capas x filas)            temperatura, K
#   ps  (filas)                    presion en superficie, Pa
# sumadas sobre los años REGISTRADOS (solo esos). La media de la funcion de corriente es la funcion de corriente de
# la media de v p_s (es lineal), asi que basta con acumular vps.
#
# Como se lee el nucleo: igual que el año extra de validar_i16.py (envolviendo NucleoSeco.avanzar y leyendo su
# resultado) -- SOLO SE LEE: el modelo no cambia (lo comprueba test_fase6_zonal.py, bit a bit). El dia y si el
# año es registrado salen del gancho AL_PASO de fase2b_atmosfera, que se encadena con el que ya hubiera.

import math

import numpy as np


def funcion_corriente(vps_zonal, sigma_semi, cos_caras, radio, gravedad):
    """Funcion de corriente de masa (kg/s) en las caras y las interfaces sigma, integrada desde el tope:
    psi_{k+1/2} = 2 pi a cos(phi) / g * sum_{j<=k} [v p_s]_j dsigma_j. La MISMA formula que
    validar_i16.funcion_corriente (lo comprueba la prueba); copiada para no importar un programa con argumentos."""
    dsig = np.diff(sigma_semi)
    acum = np.concatenate([np.zeros((1,) + vps_zonal.shape[1:]),
                           np.cumsum(vps_zonal * dsig.reshape((-1,) + (1,) * (vps_zonal.ndim - 1)), axis=0)])
    return 2 * math.pi * radio * cos_caras / gravedad * acum


class MediasZonales:
    """Sumas por dia del año de las medias zonales del nucleo durante UNA llamada a simular_fase2b (solo en el año
    registrado). Uso:
        z = MediasZonales(pasos_dia); z.instalar(); ... simular_fase2b(...) ...; z.quitar(); z.sumas()"""

    def __init__(self, pasos_dia):
        self.pasos_dia = int(pasos_dia)
        self.registrado = False
        self.dia = 0
        self.s = None             # dict de arrays (dias, ...) y "n" (dias,)
        self.geo = None
        self._orig_avanzar = None
        self._orig_al_paso = None

    def _al_paso(self, ano, p, total, registrado):
        self.registrado = bool(registrado)
        self.dia = p // self.pasos_dia
        self._dias = -(-total // self.pasos_dia)
        if self._orig_al_paso is not None:
            self._orig_al_paso(ano, p, total, registrado)

    def instalar(self):
        import fase2b_atmosfera as F
        import fase6_nucleo as N6
        self._F, self._N6 = F, N6
        self._orig_al_paso = F.AL_PASO
        self._orig_avanzar = N6.NucleoSeco.avanzar
        F.AL_PASO = self._al_paso
        original, yo = self._orig_avanzar, self

        def avanzar(nucleo, ant, act, forzamiento=None, **kw):
            res = original(nucleo, ant, act, forzamiento, **kw)
            if yo.registrado:
                yo._sumar(nucleo, res[1])
            return res
        N6.NucleoSeco.avanzar = avanzar

    def quitar(self):
        self._F.AL_PASO = self._orig_al_paso
        self._N6.NucleoSeco.avanzar = self._orig_avanzar

    def _sumar(self, nucleo, estado):
        ps, T, u, v = estado
        uc = 0.5 * (u + np.roll(u, -1, axis=-1))                       # u en el centro (como validar_i16)
        ps_cara = 0.5 * (ps[1:] + ps[:-1])
        campos = {"u": uc.mean(axis=2), "vps": (v * ps_cara[None]).mean(axis=2), "T": T.mean(axis=2),
                  "ps": ps.mean(axis=1)}
        if self.s is None:
            self.s = {k: np.zeros((self._dias,) + x.shape) for k, x in campos.items()}
            self.s["n"] = np.zeros(self._dias)
            Rj = nucleo.Rj
            self.geo = {"sigma_semi": np.asarray(nucleo.sh, dtype=float).copy(),
                        "lat_caras": np.degrees(Rj.phi_f[1:-1]), "cos_caras": np.cos(Rj.phi_f[1:-1]),
                        "lat": np.degrees(Rj.phi_c), "radio": float(Rj.a), "gravedad": float(nucleo.g)}
        d = min(self.dia, self._dias - 1)
        for k, x in campos.items():
            self.s[k][d] += x
        self.s["n"][d] += 1

    def sumas(self):
        """Las sumas por dia (o None si no hubo año registrado) y la geometria."""
        return self.s, self.geo


def combinar(lista):
    """Medias por dia de varias llamadas (años) a partir de sus sumas. lista: [(sumas, geo), ...]."""
    total, geo = None, None
    for s, g in lista:
        if s is None:
            continue
        geo = g
        if total is None:
            total = {k: v.copy() for k, v in s.items()}
        else:
            for k in total:
                total[k] += s[k]
    if total is None:
        return None
    n = np.maximum(total["n"], 1)
    medias = {k: total[k] / n.reshape((-1,) + (1,) * (total[k].ndim - 1)) for k in ("u", "vps", "T", "ps")}
    return {**medias, "subpasos_por_dia": total["n"], **geo}

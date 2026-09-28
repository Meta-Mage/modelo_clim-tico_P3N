# analisis_global.py -- Vista global del estado climatico de P3N.
#
# Herramienta de analisis (no cambia la fisica), pensada para ver "de un
# vistazo, pero completo" la situacion del planeta antes y despues de
# cada cambio del modelo (primer uso previsto: comparar Fase 2 con Fase 3).
#
# QUE CALCULA
#   - Estado global: temperatura media ponderada por area (año,
#     perihelio, afelio), hemisferio norte vs sur, y fraccion de la
#     superficie por debajo de 0 C -- global, tierra y agua por separado.
#   - Tabla por BANDAS de latitud: 60-90 N, bandas de 10 grados de 60 N a
#     60 S, y 60-90 S (14 bandas). Cada banda agrupa las filas de la
#     rejilla cuyo centro cae dentro de ella (la rejilla es de 5 grados,
#     con centros en 2.5, 7.5... -- por eso se usan bandas y no latitudes
#     sueltas). Tierra y agua SIEMPRE por separado, con el numero de
#     celdas de cada tipo en la banda.
#   - Para cada banda, tipo y dia de referencia: media ponderada por area
#     de la MINIMA diaria, de la MEDIA diaria y de la MAXIMA diaria de sus
#     celdas. Mas una fila "Año" con la media anual y los extremos
#     absolutos (celda y dia mas frio / mas calido de toda la banda).
#   - Grafico: dos paneles (tierra / agua), dia del año en horizontal,
#     latitud en vertical, color = temperatura media diaria de la franja
#     de 5 grados. Perihelio y afelio marcados.
#
# PONDERACION POR AREA: el area de una celda de la rejilla es
# proporcional a cos(latitud de su centro). Todas las medias de este
# archivo usan ese peso, para que una celda polar (pequeña) no cuente
# lo mismo que una ecuatorial (grande).
#
# ESTACIONES: los nombres (invierno, primavera...) se refieren al
# HEMISFERIO NORTE. En P3N importan poco (inclinacion axial 1.7 grados);
# la variacion anual la domina la distancia a S3N (excentricidad 0.046),
# por eso el perihelio y el afelio se muestran aparte.
#
# Usa la misma simulacion (y por tanto la misma cache) que
# analisis_latitudes.py, mapa_calor.py y consulta_punto.py: si ya se ha
# ejecutado cualquiera de ellas con el mismo mapa y parametros, este
# analisis sale en segundos.

import os
from datetime import date

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from parametros import *
from orbita import info_estaciones, anomalia_media
from rejilla import LATITUDES_GRADOS, FILAS, COLUMNAS
from fase1_geografia import TIERRA, AGUA, ALBEDO_POR_TIPO, INERCIA_POR_TIPO
from fase2_difusion import D_DIFUSION_REFERENCIA
from cache_simulacion import precalcular_orbita_cacheada, simular_rejilla_combinada_cacheada

# D de la difusion horizontal (igual que en el resto de herramientas).
D_GRID_ACTIVO = np.full((FILAS, COLUMNAS), D_DIFUSION_REFERENCIA)

# Dias de referencia por estacion, repartidos por igual en el tiempo
# real de cada estacion (inicio, y fracciones iguales a partir de ahi).
# 2 -> 8 dias (inicio y mitad de cada estacion, los de siempre).
# 3 -> 12 dias (inicio, un tercio y dos tercios de cada estacion).
DIAS_POR_ESTACION = 2

# Bandas de latitud (limite norte, limite sur), de norte a sur.
BANDAS = [(90, 60)] + [(lat, lat - 10) for lat in range(60, -60, -10)] + [(-60, -90)]

CARPETA_RESULTADOS = os.path.join("outputs", "analisis_global")

PESO_AREA = np.cos(np.radians(LATITUDES_GRADOS)).reshape(-1, 1) * np.ones((1, COLUMNAS))  # (FILAS, COLUMNAS)


# ---------------------------------------------------------------------
# Dias de referencia, perihelio y afelio
# ---------------------------------------------------------------------

def obtener_dias_referencia(num_dias, dias_por_estacion=DIAS_POR_ESTACION):
    """
    Devuelve una lista de (nombre, indice_de_registro). El indice i del
    registro corresponde al dia i+1 del año (mismo criterio que
    analisis_latitudes.obtener_dias_especiales, con el que coincide
    exactamente cuando dias_por_estacion = 2).
    """
    if dias_por_estacion == 2:
        etiquetas = ["Inicio", "Mitad"]
    elif dias_por_estacion == 3:
        etiquetas = ["Inicio", "1/3", "2/3"]
    else:
        etiquetas = ["Inicio"] + [f"{k}/{dias_por_estacion}" for k in range(1, dias_por_estacion)]

    dias = []
    for est in info_estaciones():
        for k in range(dias_por_estacion):
            dia_real = est["dia_inicio"] + est["duracion"] * k / dias_por_estacion
            indice = int(round(dia_real - 1)) % num_dias
            dias.append((f"{etiquetas[k]} de {est['estacion']}", indice))
    return dias


def dias_perihelio_afelio(num_dias):
    """
    Indices de registro del perihelio (anomalia media = 0) y del afelio
    (anomalia media = 180 grados), a partir de la anomalia media del
    instante inicial del año (dia 1, hora 0) -- misma fuente que usa
    todo el modelo, sin repetir constantes a mano.
    """
    am_inicio = anomalia_media(1, 0)  # radianes
    t_perihelio = ((2 * PI - am_inicio) % (2 * PI)) / (2 * PI) * ORBITA_PERIODO
    t_afelio = ((PI - am_inicio) % (2 * PI)) / (2 * PI) * ORBITA_PERIODO
    indice_perihelio = int(t_perihelio // 86400) % num_dias
    indice_afelio = int(t_afelio // 86400) % num_dias
    return indice_perihelio, indice_afelio


# ---------------------------------------------------------------------
# Medias ponderadas
# ---------------------------------------------------------------------

def media_ponderada(valores, mascara):
    """Media de valores (FILAS, COLUMNAS) ponderada por area, solo en las
    celdas donde mascara es True. NaN si la mascara esta vacia."""
    pesos = PESO_AREA * mascara
    total = pesos.sum()
    if total == 0:
        return np.nan
    return float((valores * pesos).sum() / total)


def fraccion_area(condicion, mascara):
    """Fraccion del area (de las celdas de la mascara) donde condicion es True."""
    pesos = PESO_AREA * mascara
    total = pesos.sum()
    if total == 0:
        return np.nan
    return float((pesos * condicion).sum() / total)


def mascara_banda(lat_norte, lat_sur):
    filas = (LATITUDES_GRADOS < lat_norte) & (LATITUDES_GRADOS > lat_sur)
    return filas.reshape(-1, 1) * np.ones((1, COLUMNAS), dtype=bool)


# ---------------------------------------------------------------------
# Texto
# ---------------------------------------------------------------------

def fmt(x, ancho=6):
    return f"{'--':>{ancho}}" if np.isnan(x) else f"{x:>{ancho}.1f}"


def pct(x):
    return "  --" if np.isnan(x) else f"{100 * x:>3.0f}%"


def nombre_banda(lat_norte, lat_sur):
    def g(lat):
        if lat == 0:
            return "0"
        return f"{abs(lat)}{'N' if lat > 0 else 'S'}"
    return f"{g(lat_norte)}-{g(lat_sur)}"


def resumen_global(registro_minima, registro_media, tipo_superficie, indice_perihelio, indice_afelio):
    todo = np.ones((FILAS, COLUMNAS), dtype=bool)
    tierra = tipo_superficie == TIERRA
    agua = tipo_superficie == AGUA
    norte = (LATITUDES_GRADOS > 0).reshape(-1, 1) * todo
    sur = (LATITUDES_GRADOS < 0).reshape(-1, 1) * todo

    media_anual = registro_media.mean(axis=0)
    minima_anual = registro_minima.min(axis=0)
    dia_perihelio = registro_media[indice_perihelio]
    dia_afelio = registro_media[indice_afelio]

    lineas = []
    lineas.append("ESTADO GLOBAL (medias ponderadas por area)")
    lineas.append(f"{'':<38}{'Global':>8}{'Tierra':>8}{'Agua':>8}")
    for etiqueta, campo in [
        ("Temperatura media anual (C)", media_anual),
        (f"Media diaria en perihelio, dia {indice_perihelio + 1}", dia_perihelio),
        (f"Media diaria en afelio, dia {indice_afelio + 1}", dia_afelio),
    ]:
        lineas.append(f"{etiqueta:<38}" + "".join(f"{fmt(media_ponderada(campo, m), 8)}" for m in (todo, tierra, agua)))
    lineas.append(f"{'Superficie con media anual < 0 C':<38}"
                  + "".join(f"{pct(fraccion_area(media_anual < 0, m)):>8}" for m in (todo, tierra, agua)))
    lineas.append(f"{'Superficie que baja de 0 C algun dia':<38}"
                  + "".join(f"{pct(fraccion_area(minima_anual < 0, m)):>8}" for m in (todo, tierra, agua)))
    lineas.append(f"{'Fraccion de la superficie':<38}{'':>8}"
                  + "".join(f"{pct(fraccion_area(m, todo)):>8}" for m in (tierra, agua)))
    lineas.append("")
    lineas.append(f"Hemisferio norte: media anual {fmt(media_ponderada(media_anual, norte), 0)} C | "
                  f"perihelio {fmt(media_ponderada(dia_perihelio, norte), 0)} C | "
                  f"afelio {fmt(media_ponderada(dia_afelio, norte), 0)} C")
    lineas.append(f"Hemisferio sur:   media anual {fmt(media_ponderada(media_anual, sur), 0)} C | "
                  f"perihelio {fmt(media_ponderada(dia_perihelio, sur), 0)} C | "
                  f"afelio {fmt(media_ponderada(dia_afelio, sur), 0)} C")
    return lineas


def resumen_bandas_anual(registro_minima, registro_media, registro_maxima, tipo_superficie):
    """Tabla corta (una linea por banda) con los valores anuales."""
    media_anual = registro_media.mean(axis=0)
    minima_abs = registro_minima.min(axis=0)
    maxima_abs = registro_maxima.max(axis=0)

    enc = (f"{'Banda':<9}| {'Tierra':^31} | {'Agua':^31}")
    sub = (f"{'':<9}| {'celdas':>6} {'media':>7} {'min.abs':>7} {'max.abs':>7} "
           f"| {'celdas':>6} {'media':>7} {'min.abs':>7} {'max.abs':>7}")
    lineas = ["RESUMEN ANUAL POR BANDAS (media anual ponderada; extremos absolutos del año)", enc, sub, "-" * len(sub)]
    for lat_n, lat_s in BANDAS:
        mb = mascara_banda(lat_n, lat_s)
        partes = []
        for tipo in (TIERRA, AGUA):
            m = mb & (tipo_superficie == tipo)
            n = int(m.sum())
            if n == 0:
                partes.append(f"{0:>6} {'--':>7} {'--':>7} {'--':>7}")
            else:
                partes.append(f"{n:>6} {fmt(media_ponderada(media_anual, m), 7)} "
                              f"{fmt(float(minima_abs[m].min()), 7)} {fmt(float(maxima_abs[m].max()), 7)}")
        lineas.append(f"{nombre_banda(lat_n, lat_s):<9}| {partes[0]} | {partes[1]}")
    return lineas


def tabla_detallada(registro_minima, registro_media, registro_maxima, tipo_superficie, dias_referencia):
    """Tabla completa: por banda, una linea por dia de referencia + una linea 'Año'."""
    enc = (f"{'Banda':<9}| {'Dia':<28}| {'Tierra: min / med / max':^22} | {'Agua: min / med / max':^22}")
    lineas = ["DETALLE POR BANDAS Y DIAS DE REFERENCIA",
              "(min/med/max = media ponderada por area de la minima, la media y la maxima diaria de las celdas;",
              " fila 'Año': media anual, y minima/maxima absolutas del año en la banda)",
              enc, "=" * len(enc)]
    for lat_n, lat_s in BANDAS:
        mb = mascara_banda(lat_n, lat_s)
        m_t = mb & (tipo_superficie == TIERRA)
        m_a = mb & (tipo_superficie == AGUA)
        banda = f"{nombre_banda(lat_n, lat_s)} (tierra {int(m_t.sum())}, agua {int(m_a.sum())} celdas)"
        lineas.append(banda)
        for nombre_dia, i in dias_referencia:
            celdas = []
            for m in (m_t, m_a):
                celdas.append(f"{fmt(media_ponderada(registro_minima[i], m))} {fmt(media_ponderada(registro_media[i], m))} "
                              f"{fmt(media_ponderada(registro_maxima[i], m))}   ")
            lineas.append(f"{'':<9}| {nombre_dia + ' (d' + str(i + 1) + ')':<28}| {celdas[0]:<22} | {celdas[1]:<22}")
        celdas = []
        for m in (m_t, m_a):
            if m.sum() == 0:
                celdas.append(f"{'--':>6} {'--':>6} {'--':>6}   ")
            else:
                celdas.append(f"{fmt(float(registro_minima[:, m].min()))} {fmt(media_ponderada(registro_media.mean(axis=0), m))} "
                              f"{fmt(float(registro_maxima[:, m].max()))}   ")
        lineas.append(f"{'':<9}| {'Año':<28}| {celdas[0]:<22} | {celdas[1]:<22}")
        lineas.append("-" * len(enc))
    return lineas


# ---------------------------------------------------------------------
# Grafico
# ---------------------------------------------------------------------

def media_zonal_por_tipo(registro_media, tipo_superficie, tipo):
    """(FILAS, dias): media de cada fila de 5 grados solo en las celdas del
    tipo dado (NaN si la fila no tiene ninguna). Dentro de una fila todas
    las celdas tienen la misma area, asi que la media simple es correcta."""
    mascara = (tipo_superficie == tipo)  # (FILAS, COLUMNAS)
    n = mascara.sum(axis=1)  # (FILAS,)
    suma = (registro_media * mascara).sum(axis=2)  # (dias, FILAS)
    with np.errstate(invalid="ignore", divide="ignore"):
        resultado = np.where(n > 0, suma / np.maximum(n, 1), np.nan)
    return resultado.T


def generar_grafico(registro_media, tipo_superficie, dias_referencia, indice_perihelio, indice_afelio, nombre_mapa, ruta):
    tierra = media_zonal_por_tipo(registro_media, tipo_superficie, TIERRA)
    agua = media_zonal_por_tipo(registro_media, tipo_superficie, AGUA)
    num_dias = registro_media.shape[0]

    validos = np.concatenate([tierra[~np.isnan(tierra)], agua[~np.isnan(agua)]])
    vmin, vmax = np.floor(validos.min()), np.ceil(validos.max())

    cmap = plt.get_cmap("turbo").copy()
    cmap.set_bad("#d0d0d0")

    fig, ejes = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
    for eje, datos, titulo in [(ejes[0], tierra, "Tierra"), (ejes[1], agua, "Agua")]:
        im = eje.imshow(np.ma.masked_invalid(datos), aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                        extent=[0.5, num_dias + 0.5, -90, 90], origin="upper", interpolation="nearest")
        for dia, etiqueta, estilo in [(indice_perihelio + 1, "Perihelio", "-"), (indice_afelio + 1, "Afelio", "-")]:
            eje.axvline(dia, color="white", linestyle=estilo, linewidth=1.4)
            a_la_izquierda = dia > num_dias * 0.85  # cerca del borde derecho: etiqueta a la izquierda de la linea
            eje.text(dia - 2 if a_la_izquierda else dia + 2, 84, etiqueta, color="white", fontsize=9, va="top",
                     ha="right" if a_la_izquierda else "left")
        for _, i in dias_referencia:
            eje.axvline(i + 1, color="black", linestyle=":", linewidth=0.7, alpha=0.6)
        for lat in (60, -60):
            eje.axhline(lat, color="black", linewidth=0.5, alpha=0.4)
        eje.axhline(0, color="black", linewidth=0.8, alpha=0.6)
        eje.set_title(titulo)
    ejes[0].set_ylabel("Latitud (grados)")
    ejes[0].set_yticks(range(-90, 91, 30))
    fig.colorbar(im, ax=ejes, label="Temperatura media diaria (C)  (gris: sin celdas de ese tipo)")
    fig.supxlabel("Dia del año  (punteado: dias de referencia; estaciones del hemisferio norte)")
    fig.suptitle(f"Evolucion anual por latitud, tierra y agua por separado -- {nombre_mapa}")
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------

if __name__ == "__main__":
    datos_orbita = precalcular_orbita_cacheada(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    from puente_c3n import cargar_mapa_activo_de_c3n
    tipo_superficie, altitud_metros, nombre_mapa = cargar_mapa_activo_de_c3n()

    print(f"Simulando mapa '{nombre_mapa}' (si ya esta en cache, es instantaneo)...")
    _, anos_convergencia, registro_minima, registro_media, registro_maxima = simular_rejilla_combinada_cacheada(
        datos_orbita, tipo_superficie, altitud_metros, EMISIVIDAD,
        ALBEDO_POR_TIPO, INERCIA_POR_TIPO, PROFUNDIDAD_OPTICA, D_GRID_ACTIVO,
        nombre_mapa=nombre_mapa,
    )
    num_dias = registro_media.shape[0]
    dias_referencia = obtener_dias_referencia(num_dias)
    indice_perihelio, indice_afelio = dias_perihelio_afelio(num_dias)

    cabecera = [
        f"ANALISIS GLOBAL -- mapa '{nombre_mapa}' -- {date.today().isoformat()}",
        f"D = {D_DIFUSION_REFERENCIA} W/m2/K | convergencia: {anos_convergencia} año(s) | año de {num_dias} dias",
    ]
    if anos_convergencia >= 50:
        cabecera.append("AVISO: se alcanzo el limite de 50 años sin confirmar convergencia -- puede no ser el equilibrio.")
    cabecera.append("")

    global_ = resumen_global(registro_minima, registro_media, tipo_superficie, indice_perihelio, indice_afelio)
    anual = resumen_bandas_anual(registro_minima, registro_media, registro_maxima, tipo_superficie)
    detalle = tabla_detallada(registro_minima, registro_media, registro_maxima, tipo_superficie, dias_referencia)

    # En la terminal: estado global + resumen anual por bandas. El detalle
    # por dias (largo) va solo al archivo.
    print("\n".join(cabecera + global_ + [""] + anual))

    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)
    base = os.path.join(CARPETA_RESULTADOS, f"analisis_global_{nombre_mapa}_{date.today().isoformat()}")
    with open(base + ".txt", "w", encoding="utf-8") as f:
        f.write("\n".join(cabecera + global_ + [""] + anual + [""] + detalle) + "\n")
    generar_grafico(registro_media, tipo_superficie, dias_referencia, indice_perihelio, indice_afelio,
                    nombre_mapa, base + ".png")
    print(f"\nTabla completa (con el detalle por dias) guardada en {base}.txt")
    print(f"Grafico guardado en {base}.png")

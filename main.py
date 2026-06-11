from parametros import *
from orbita import *
from geometria import *
from atmosfera import *
from radiacion import *
from temperatura import *
import math

# PARÁMETROS FIJOS
MASA = 0.97

# RANGOS DEL BARRIDO
INCLINACIONES = [round(i * 1.0, 1) for i in range(46)]   # 0° a 45°
LATITUDES = [round(-60 + i * 5, 0) for i in range(25)]   # -60° a 60°

# CRITERIOS DE HABITABILIDAD
T_SUPERVIVENCIA_MAX = 45
T_SUPERVIVENCIA_MIN = -50
T_HABITABLE_MAX = 40
T_HABITABLE_MIN = -10
T_CONFORT_MAX = 30
T_CONFORT_MIN = 0
T_OPTIMA_MAX = 20
T_OPTIMA_MIN = 10
T_PICO_MAX = 40
T_PICO_MIN = -20


def puntuacion_latitud(ciclo):
    media = sum(ciclo) / len(ciclo)
    maximo = max(ciclo)
    minimo = min(ciclo)

    if maximo > T_SUPERVIVENCIA_MAX or minimo < T_SUPERVIVENCIA_MIN:
        return 0, "descartada", f"máx {maximo:.1f}°C / mín {minimo:.1f}°C fuera de supervivencia"

    frac_hab = sum(1 for t in ciclo if T_HABITABLE_MIN <= t <= T_HABITABLE_MAX) / len(ciclo)
    puntos_hab = frac_hab * 50

    frac_conf = sum(1 for t in ciclo if T_CONFORT_MIN <= t <= T_CONFORT_MAX) / len(ciclo)
    puntos_conf = frac_conf * 25

    if T_OPTIMA_MIN <= media <= T_OPTIMA_MAX:
        puntos_media = 25
    elif media < T_OPTIMA_MIN:
        puntos_media = max(0, 25 - (T_OPTIMA_MIN - media) * 1.5)
    else:
        puntos_media = max(0, 25 - (media - T_OPTIMA_MAX) * 1.5)

    picos_calor = sum(1 for t in ciclo if t > T_PICO_MAX)
    picos_frio = sum(1 for t in ciclo if t < T_PICO_MIN)
    fraccion_picos = (picos_calor + picos_frio) / len(ciclo)
    penalizacion = fraccion_picos * 50

    puntuacion = max(0, puntos_hab + puntos_conf + puntos_media - penalizacion)

    if puntuacion >= 75:
        estado = "optima"
    elif puntuacion >= 50:
        estado = "habitable"
    elif puntuacion >= 25:
        estado = "marginal"
    else:
        estado = "no_habitable"

    motivo = (f"media {media:.1f}°C | "
              f"habitable {frac_hab*100:.0f}% | "
              f"confort {frac_conf*100:.0f}% | "
              f"picos {fraccion_picos*100:.1f}% | "
              f"punt {puntuacion:.1f}")

    return puntuacion, estado, motivo


def dias_confortables(ciclo, pasos_por_dia):
    dias = 0
    for d in range(270):
        inicio = d * pasos_por_dia
        fin = inicio + pasos_por_dia
        temps = ciclo[inicio:fin]
        media_dia = sum(temps) / len(temps)
        if T_CONFORT_MIN <= media_dia <= T_CONFORT_MAX:
            dias += 1
    return dias


def calcular_parametros_estelares(masa):
    radio = RADIO_SOL * (masa ** EXP_MASA_RADIO)
    temperatura = TEMPERATURA_SOL * (masa ** EXP_MASA_TEMPERATURA)
    luminosidad = 4 * PI * radio**2 * CONSTANTE_SB * temperatura**4
    masa_kg = masa * MASA_SOL
    semieje = (CONSTANTE_GRAVITACIONAL * masa_kg * ORBITA_PERIODO**2 / (4 * PI**2)) ** (1/3)
    return luminosidad, semieje


def barrido():
    resultados_optimos = []
    todas_combinaciones = []
    total = len(INCLINACIONES)
    combinacion = 0
    pasos_por_dia = int(86400 / PASO_TIEMPO)
    dias_muestra = [1, 28, 55, 82, 109, 136, 163, 190, 217, 244]

    luminosidad, semieje = calcular_parametros_estelares(MASA)

    html = open("resultados.html", "w", encoding="utf-8")
    html.write(f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Barrido Inclinación Axial P3N</title>
<style>
  body {{ font-family: monospace; background: #1a1a2e; color: #e0e0e0; padding: 20px; }}
  details {{ margin: 8px 0; border: 1px solid #444; border-radius: 4px; }}
  summary {{ padding: 10px; cursor: pointer; background: #16213e; }}
  summary:hover {{ background: #0f3460; }}
  .optima {{ color: #00ff88; }}
  .habitable {{ color: #ffcc00; }}
  .marginal {{ color: #ff9900; }}
  .descartada {{ color: #ff4444; }}
  .no_habitable {{ color: #aaaaaa; }}
  .contenido {{ padding: 10px 20px; }}
  table {{ border-collapse: collapse; margin: 8px 0; }}
  td, th {{ padding: 4px 10px; border: 1px solid #444; }}
  th {{ background: #0f3460; }}
  .resumen {{ background: #16213e; padding: 10px; margin: 10px 0; border-radius: 4px; }}
  .params {{ background: #0d1b2a; padding: 10px; margin-bottom: 20px; border-radius: 4px; color: #aaddff; }}
</style>
</head>
<body>
<h1>Barrido Inclinación Axial P3N</h1>
<div class="params">
<b>Parámetros fijos:</b><br>
Masa: {MASA} M☉ | Excentricidad: {ORBITA_EXCENTRICIDAD} | Albedo: {ALBEDO} |
Emisividad: {EMISIVIDAD} | Profundidad óptica: {PROFUNDIDAD_OPTICA} |
Inercia térmica: {INERCIA_TERMICA} SI | Exp. masa-T: {EXP_MASA_TEMPERATURA}<br>
<b>Rango barrido:</b> Inclinación axial {INCLINACIONES[0]}° a {INCLINACIONES[-1]}° cada 1°
</div>
""")

    for inclinacion in INCLINACIONES:
        combinacion += 1
        inclinacion_rad = inclinacion * PI / 180
        datos = precalcular_orbita(luminosidad, inclinacion_rad, semieje)

        print(f"[{combinacion}/{total}] Inclinación: {inclinacion}°...")

        latitudes_optimas = []
        latitudes_habitables = []
        latitudes_marginales = []
        latitudes_descartadas = []
        puntuaciones_validas = []
        contenido_lats = ""

        for latitud in LATITUDES:
            latitud_rad = latitud * PI / 180
            T_final, ciclo = simular(latitud_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
            punt, estado, motivo = puntuacion_latitud(ciclo)
            dc = dias_confortables(ciclo, pasos_por_dia)

            if estado == "optima":
                latitudes_optimas.append(latitud)
                puntuaciones_validas.append(punt)
            elif estado == "habitable":
                latitudes_habitables.append(latitud)
                puntuaciones_validas.append(punt)
            elif estado == "marginal":
                latitudes_marginales.append(latitud)
                puntuaciones_validas.append(punt)
            elif estado == "descartada":
                latitudes_descartadas.append(latitud)

            filas = ""
            for dia in dias_muestra:
                inicio = (dia - 1) * pasos_por_dia
                fin = min(inicio + pasos_por_dia, len(ciclo))
                temps = ciclo[inicio:fin]
                media_dia = sum(temps) / len(temps)
                max_dia = max(temps)
                min_dia = min(temps)
                filas += f"<tr><td>{dia}</td><td>{media_dia:.2f}</td><td>{max_dia:.2f}</td><td>{min_dia:.2f}</td></tr>"

            contenido_lats += f"""
            <details>
            <summary class="{estado}">Lat {latitud}° → {estado.upper()} | {motivo} | Días confort: {dc}</summary>
            <div class="contenido">
            <table>
            <tr><th>Día</th><th>Media</th><th>Máxima</th><th>Mínima</th></tr>
            {filas}
            </table>
            </div>
            </details>
            """

        punt_media = sum(puntuaciones_validas) / len(puntuaciones_validas) if puntuaciones_validas else 0
        color_resumen = "optima" if latitudes_optimas else ("habitable" if latitudes_habitables else ("marginal" if latitudes_marginales else "no_habitable"))

        html.write(f"""
<details>
<summary class="{color_resumen}">[{combinacion}/{total}] Inclinación: {inclinacion}° |
Óptimas: {len(latitudes_optimas)} | Habitables: {len(latitudes_habitables)} |
Marginales: {len(latitudes_marginales)} | Punt. media: {punt_media:.1f}/100</summary>
<div class="contenido">
<div class="resumen">
Latitudes óptimas: {latitudes_optimas}<br>
Latitudes habitables: {latitudes_habitables}<br>
Latitudes marginales: {latitudes_marginales}<br>
Latitudes descartadas: {latitudes_descartadas}
</div>
{contenido_lats}
</div>
</details>
""")

        todas_combinaciones.append({
            "inclinacion": inclinacion,
            "latitudes_optimas": latitudes_optimas,
            "latitudes_habitables": latitudes_habitables,
            "latitudes_marginales": latitudes_marginales,
            "latitudes_descartadas": latitudes_descartadas,
            "puntuacion_media": punt_media
        })

        if latitudes_optimas or latitudes_habitables:
            resultados_optimos.append(todas_combinaciones[-1])

    html.write("</body></html>")
    html.close()
    print("\nResultados guardados en resultados.html")

    return resultados_optimos, todas_combinaciones


def resumen(resultados_optimos, todas_combinaciones):
    print("\n" + "="*60)
    print("RESUMEN DEL BARRIDO")
    print("="*60)

    if resultados_optimos:
        print(f"Inclinaciones con latitudes óptimas o habitables: {len(resultados_optimos)}")
        mejor = max(resultados_optimos, key=lambda x: (len(x["latitudes_optimas"]), x["puntuacion_media"]))
        print(f"\nMejor inclinación:")
        print(f"  Inclinación: {mejor['inclinacion']}°")
        print(f"  Latitudes óptimas: {mejor['latitudes_optimas']}")
        print(f"  Latitudes habitables: {mejor['latitudes_habitables']}")
        print(f"  Puntuación media: {mejor['puntuacion_media']:.1f}/100")
    else:
        print("No se encontraron combinaciones óptimas o habitables.")

    print("\n--- TOP 5 POR PUNTUACIÓN MEDIA ---")
    top5 = sorted(todas_combinaciones, key=lambda x: x["puntuacion_media"], reverse=True)[:5]
    for i, comb in enumerate(top5):
        print(f"{i+1}. Inclinación {comb['inclinacion']}° → "
              f"Punt: {comb['puntuacion_media']:.1f}/100 | "
              f"Óptimas: {len(comb['latitudes_optimas'])} | "
              f"Habitables: {len(comb['latitudes_habitables'])} | "
              f"Marginales: {len(comb['latitudes_marginales'])}")


if __name__ == "__main__":
    print("Iniciando barrido de inclinación axial...")
    print(f"Masa fija: {MASA} M☉")
    print(f"Inclinaciones: {INCLINACIONES[0]}° a {INCLINACIONES[-1]}° cada 1°")
    print(f"Latitudes: {LATITUDES}")
    print(f"Total combinaciones: {len(INCLINACIONES)}")
    print("="*60)
    resultados_optimos, todas_combinaciones = barrido()
    resumen(resultados_optimos, todas_combinaciones)
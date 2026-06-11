from parametros import *
from orbita import *
from geometria import *
from atmosfera import *
from radiacion import *
import math
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\sala.AULASUC-214VNRO\ffmpeg\ffmpeg-master-latest-win64-gpl-shared\bin\ffmpeg.exe'

PASO_TIEMPO = 900
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


# ============================================================
# UTILIDADES
# ============================================================

def texto_parametros(dia=None, hora=None, latitud=None):
    lineas = [
        f"Masa: {S3N_MASAS_SOLARES} M☉ | Exp.T: {EXP_MASA_TEMPERATURA} | "
        f"Emis: {EMISIVIDAD} | Albedo: {ALBEDO} | Exc: {ORBITA_EXCENTRICIDAD} | "
        f"Inercia: {INERCIA_TERMICA} | Incl: {INCLINACION_AXIAL}°"
    ]
    if dia is not None:
        lineas.append(f"Día: {dia} | Hora: {hora}h | Latitud: {latitud}°")
    return "\n".join(lineas)


def nombre_archivo(nombre):
    return f"outputs/{nombre}_{TIMESTAMP}"


# ============================================================
# NÚCLEO DEL MODELO
# ============================================================

def precalcular_orbita(luminosidad, inclinacion_axial_rad, semieje):
    pasos = round(ORBITA_PERIODO / PASO_TIEMPO)
    datos = []
    for paso in range(pasos):
        t = paso * PASO_TIEMPO
        dia = int(t / 86400) + 1
        hora = (t % 86400) / 3600
        AM = anomalia_media(dia, hora)
        AE = anomalia_excentrica(AM)
        AV = anomalia_verdadera(AE)
        dist = distancia_S3N(AV, semieje)
        decl = declinacion_solar(AV, inclinacion_axial_rad)
        ang_h = angulo_horario(hora)
        toa = i_toa(dist, luminosidad)
        datos.append((toa, decl, ang_h))
    return datos


def t_eq(irradiancia_absorbida):
    return (irradiancia_absorbida / CONSTANTE_SB) ** 0.25


def simular(latitud_rad, datos_orbita, emisividad, inercia, albedo, profundidad_optica):
    C = inercia * math.sqrt(ROTACION_PERIODO / PI)
    toa_ini, decl_ini, ang_h_ini = datos_orbita[len(datos_orbita)//2]
    cenital_ini = angulo_cenital(latitud_rad, decl_ini, ang_h_ini)
    inst_ini = i_inst(toa_ini, cenital_ini)
    masa_ini = masa_aire(cenital_ini)
    if masa_ini is None:
        T = 273.15
    else:
        tra_ini = trans(masa_ini, profundidad_optica)
        atm_ini = i_atm(inst_ini, tra_ini)
        abs_ini = i_abs(atm_ini, albedo)
        T = t_eq(abs_ini)

    for año in range(50):
        T_inicio_año = T
        registro_año = []
        for toa, decl, ang_h in datos_orbita:
            cenital = angulo_cenital(latitud_rad, decl, ang_h)
            inst = i_inst(toa, cenital)
            masa = masa_aire(cenital)
            if masa is None:
                abs_local = 0
            else:
                tra = trans(masa, profundidad_optica)
                atm = i_atm(inst, tra)
                abs_local = i_abs(atm, albedo)
            dT = (PASO_TIEMPO / C) * (abs_local - (1 - emisividad/2) * CONSTANTE_SB * T**4)
            T = T + dT
            registro_año.append(T - 273.15)
        if abs(T - T_inicio_año) < 0.01:
            break

    return T - 273.15, registro_año


# ============================================================
# OPCIÓN 1 — MOMENTO CONCRETO
# ============================================================

def opcion_momento(datos):
    dia = int(input("Día del año (1-270): "))
    hora = int(input("Hora del día (0-23): "))
    latitud = float(input("Latitud (-90 a 90): "))
    latitud_rad = latitud * PI / 180

    print(f"\nConsulta: Día {dia}, Hora {hora}h, Latitud {latitud}°")
    print("="*60)

    T_final, ciclo = simular(latitud_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)

    pasos_por_dia = int(86400 / PASO_TIEMPO)
    inicio_dia = (dia - 1) * pasos_por_dia
    paso_hora = int(hora * 3600 / PASO_TIEMPO)
    temp_momento = ciclo[inicio_dia + paso_hora]

    # Estación
    AV_dia = anomalia_verdadera(anomalia_excentrica(anomalia_media(dia, hora)))
    est = estacion(AV_dia)

    # Distancia al Sol
    dist_dia = distancia_S3N(AV_dia, SEMIEJE_MAYOR)

    # Día o noche
    decl_dia = declinacion_solar(AV_dia, INCLINACION_AXIAL_RAD)
    ang_h_dia = angulo_horario(hora)
    cenital_dia = angulo_cenital(latitud_rad, decl_dia, ang_h_dia)
    momento_dia = "Día" if cenital_dia < PI/2 else "Noche"

    print(f"Temperatura en ese momento: {temp_momento:.2f}°C")
    print(f"Estación: {est}")
    print(f"Distancia a S3N: {dist_dia/1e9:.3f} millones de km")
    print(f"Momento: {momento_dia}")

    # Ciclo del día
    print(f"\nCiclo completo del día {dia}:")
    for i in range(pasos_por_dia):
        if i % 4 == 0:
            h = i / 4
            temp = ciclo[inicio_dia + i]
            print(f"  Hora {h:.0f}: {temp:.2f}°C")

    # Climograma anual para esa latitud
    dias_año = list(range(1, 271))
    medias_año, maximas_año, minimas_año = [], [], []
    for d in dias_año:
        inicio = (d - 1) * pasos_por_dia
        fin = inicio + pasos_por_dia
        temps = ciclo[inicio:fin]
        medias_año.append(sum(temps) / len(temps))
        maximas_año.append(max(temps))
        minimas_año.append(min(temps))

    plt.figure(figsize=(14, 6))
    plt.plot(dias_año, maximas_año, color='tomato', linewidth=1.5, label='Máxima diaria')
    plt.plot(dias_año, medias_año, color='gold', linewidth=2, label='Media diaria')
    plt.plot(dias_año, minimas_año, color='steelblue', linewidth=1.5, label='Mínima diaria')
    plt.fill_between(dias_año, minimas_año, maximas_año, alpha=0.15, color='gray')
    plt.axvline(dia, color='white', linestyle='--', linewidth=1.5, label=f'Día {dia}')
    plt.axhline(0, color='white', linestyle='--', linewidth=0.8)
    plt.axhline(35, color='orange', linestyle=':', linewidth=0.8, label='Límite tolerable')
    plt.axhline(-5, color='cyan', linestyle=':', linewidth=0.8)
    plt.xlabel("Día del año")
    plt.ylabel("Temperatura (°C)")
    plt.title(f"Climograma anual — Latitud {latitud}°")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.figtext(0.5, -0.05, texto_parametros(dia=dia, hora=hora, latitud=latitud),
                ha='center', fontsize=7, color='gray')
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(f"{nombre_archivo(f'climograma_lat{latitud}_dia{dia}')}.png")
    plt.show()
    plt.close()

    # Ciclo diario
    horas = [i * 15 / 60 for i in range(pasos_por_dia)]
    temps_dia = ciclo[inicio_dia:inicio_dia + pasos_por_dia]

    plt.figure(figsize=(10, 5))
    plt.plot(horas, temps_dia, color='tomato', linewidth=2)
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axvline(hora, color='white', linestyle='--', linewidth=1.5, label=f'Hora {hora}h')
    plt.xlabel("Hora del día")
    plt.ylabel("Temperatura (°C)")
    plt.title(f"Ciclo diario — Latitud {latitud}°, Día {dia} ({est})")
    plt.xticks(range(0, 25, 2))
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.figtext(0.5, -0.05, texto_parametros(dia=dia, hora=hora, latitud=latitud),
                ha='center', fontsize=7, color='gray')
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(f"{nombre_archivo(f'ciclo_diario_lat{latitud}_dia{dia}')}.png")
    plt.show()
    plt.close()

    # Comparativo de latitudes — ambos hemisferios
    latitudes_comp = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
    plt.figure(figsize=(12, 6))
    colores = ['#4fc3f7', '#81d4fa', '#b3e5fc', '#e0f7fa', '#ffffff',
               '#ffcc80', '#ffa726', '#f57c00', '#e65100']

    for i, lat in enumerate(latitudes_comp):
        lat_rad = lat * PI / 180
        _, ciclo_lat = simular(lat_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        inicio = (dia - 1) * pasos_por_dia
        temps = ciclo_lat[inicio:inicio + pasos_por_dia]
        estilo = '-' if lat >= 0 else '--'
        plt.plot(horas, temps, color=colores[i], linewidth=2,
                 linestyle=estilo, label=f"Lat {lat}°")

    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axhline(35, color='orange', linestyle=':', linewidth=0.8, label='Límite superior tolerable')
    plt.axhline(-5, color='cyan', linestyle=':', linewidth=0.8, label='Límite inferior tolerable')
    plt.axvline(hora, color='white', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlabel("Hora del día")
    plt.ylabel("Temperatura (°C)")
    plt.title(f"Ciclo diario comparativo — Día {dia} (— hemisferio N, -- hemisferio S)")
    plt.xticks(range(0, 25, 2))
    plt.legend(fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.figtext(0.5, -0.05, texto_parametros(dia=dia, hora=hora, latitud=latitud),
                ha='center', fontsize=7, color='gray')
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(f"{nombre_archivo(f'ciclo_comparativo_dia{dia}')}.png")
    plt.show()
    plt.close()

    # Vídeo mapa de calor 24 horas
    import matplotlib.animation as animation
    import numpy as np

    print("Generando vídeo mapa de calor...")
    lats = list(range(-90, 91, 10))
    lons = list(range(-180, 181, 10))
    pasos_hora = int(3600 / PASO_TIEMPO)
    hora_inicio_pasos = int(hora * 3600 / PASO_TIEMPO)
    mapa = np.zeros((len(lats), len(lons), 24))

    for i, lat in enumerate(lats):
        lat_rad = lat * PI / 180
        _, ciclo_lat = simular(lat_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        inicio_dia_v = (dia - 1) * int(86400 / PASO_TIEMPO)
        for j, lon in enumerate(lons):
            for h in range(24):
                desfase = int(lon / 360 * 86400 / PASO_TIEMPO)
                paso = inicio_dia_v + hora_inicio_pasos + h * pasos_hora + desfase
                paso = paso % len(ciclo_lat)
                mapa[i, j, h] = ciclo_lat[paso]

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    im = ax.imshow(mapa[:, :, 0], aspect='auto', origin='lower',
                   extent=[-180, 180, -90, 90], vmin=-40, vmax=60, cmap='turbo')
    plt.colorbar(im, ax=ax, label='Temperatura (°C)')
    ax.set_xlabel("Longitud (°)", color='white')
    ax.set_ylabel("Latitud (°)", color='white')
    ax.tick_params(colors='white')
    titulo = ax.set_title(f"P3N — Día {dia}, Hora {hora:.0f}h", color='white')
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    ax.grid(True, alpha=0.3, color='white', linewidth=0.5)

    def actualizar(h):
        im.set_data(mapa[:, :, h])
        hora_actual = (hora + h) % 24
        titulo.set_text(f"P3N — Día {dia}, Hora {hora_actual:.0f}h")
        return [im, titulo]

    ani = animation.FuncAnimation(fig, actualizar, frames=24, interval=200, blit=True)
    writer = animation.FFMpegWriter(fps=4, bitrate=1800)
    ani.save(f"{nombre_archivo(f'mapa_calor_dia{dia}_hora{hora:.0f}')}.mp4",
             writer=writer, dpi=100, savefig_kwargs={'facecolor': 'black'})
    plt.close()
    print("Vídeo diario guardado.")


# ============================================================
# OPCIÓN 2 — VISIÓN ANUAL GLOBAL
# ============================================================

def opcion_anual(datos):
    import matplotlib.animation as animation
    import numpy as np

    latitudes_consulta = list(range(-60, 65, 5))
    pasos_por_dia = int(86400 / PASO_TIEMPO)

    # Calcular ciclos para todas las latitudes
    print("Calculando ciclos para todas las latitudes...")
    ciclos_por_latitud = {}
    stats_por_latitud = {}

    for lat in latitudes_consulta:
        lat_rad = lat * PI / 180
        T_final, ciclo = simular(lat_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        ciclos_por_latitud[lat] = ciclo

        media = sum(ciclo) / len(ciclo)
        maximo = max(ciclo)
        minimo = min(ciclo)

        # Días confortables (media diaria entre 0 y 30°C)
        dias_confort = 0
        for d in range(270):
            inicio = d * pasos_por_dia
            fin = inicio + pasos_por_dia
            temps = ciclo[inicio:fin]
            media_dia = sum(temps) / len(temps)
            if 0 <= media_dia <= 30:
                dias_confort += 1

        # Puntuación
        if maximo > 60 or minimo < -50:
            punt, estado = 0, "DESCARTADA"
        else:
            frac_hab = sum(1 for t in ciclo if -10 <= t <= 40) / len(ciclo)
            frac_conf = sum(1 for t in ciclo if 0 <= t <= 30) / len(ciclo)
            if 10 <= media <= 20:
                p_media = 25
            elif media < 10:
                p_media = max(0, 25 - (10 - media) * 1.5)
            else:
                p_media = max(0, 25 - (media - 20) * 1.5)
            punt = frac_hab * 50 + frac_conf * 25 + p_media
            if punt >= 75:
                estado = "ÓPTIMA"
            elif punt >= 50:
                estado = "HABITABLE"
            elif punt >= 25:
                estado = "MARGINAL"
            else:
                estado = "no_hab"

        stats_por_latitud[lat] = {
            "media": media, "maximo": maximo, "minimo": minimo,
            "punt": punt, "estado": estado, "dias_confort": dias_confort
        }

    # ---- RESUMEN ----
    print("\n" + "="*75)
    print("RESUMEN GLOBAL")
    print("="*75)
    print(f"{'Lat':>5} | {'Media':>7} | {'Máx':>7} | {'Mín':>7} | {'Punt':>6} | {'Días confort':>12} | Estado")
    print("-"*75)
    for lat in latitudes_consulta:
        s = stats_por_latitud[lat]
        print(f"{lat:>5} | {s['media']:>7.2f} | {s['maximo']:>7.2f} | {s['minimo']:>7.2f} | "
              f"{s['punt']:>6.1f} | {s['dias_confort']:>12} | {s['estado']}")
    print("="*75)

    # ---- DETALLE ----
    print("\n" + "="*75)
    print("DETALLE POR LATITUD")
    print("="*75)
    for lat in latitudes_consulta:
        s = stats_por_latitud[lat]
        ciclo = ciclos_por_latitud[lat]
        print(f"\n--- Latitud {lat}° | {s['estado']} | Puntuación: {s['punt']:.1f} | "
              f"Media: {s['media']:.2f}°C | Máx: {s['maximo']:.2f}°C | Mín: {s['minimo']:.2f}°C | "
              f"Días confort: {s['dias_confort']} ---")
        print(f"  {'Día':>4} | {'Media':>7} | {'Máxima':>7} | {'Mínima':>7}")
        print(f"  {'-'*33}")
        for dia in range(1, 271, 10):
            inicio = (dia - 1) * pasos_por_dia
            fin = min(inicio + pasos_por_dia * 10, len(ciclo))
            temps = ciclo[inicio:fin]
            m = sum(temps) / len(temps)
            mx = max(temps)
            mn = min(temps)
            print(f"  {dia:>4} | {m:>7.2f} | {mx:>7.2f} | {mn:>7.2f}")

    # ---- GRÁFICOS ----
    dias_año = list(range(1, 271))
    colores_n = ['#ffcc80', '#ffa726', '#f57c00', '#e65100', '#bf360c',
                 '#ff8a65', '#ff7043', '#ff5722', '#dd2c00', '#b71c1c', '#ff1744', '#d50000', '#ff6d00']
    colores_s = ['#4fc3f7', '#29b6f6', '#039be5', '#0288d1', '#0277bd',
                 '#81d4fa', '#b3e5fc', '#e1f5fe', '#80deea', '#4dd0e1', '#26c6da', '#00bcd4', '#00acc1']

    plt.figure(figsize=(16, 8))
    lats_positivas = [l for l in latitudes_consulta if l >= 0]
    lats_negativas = [l for l in latitudes_consulta if l < 0]

    for i, lat in enumerate(lats_positivas):
        ciclo = ciclos_por_latitud[lat]
        medias_lat = [sum(ciclo[(d-1)*pasos_por_dia:(d-1)*pasos_por_dia+pasos_por_dia]) /
                      pasos_por_dia for d in dias_año]
        plt.plot(dias_año, medias_lat, color=colores_n[i % len(colores_n)],
                 linewidth=2, label=f"N {lat}°")

    for i, lat in enumerate(lats_negativas):
        ciclo = ciclos_por_latitud[lat]
        medias_lat = [sum(ciclo[(d-1)*pasos_por_dia:(d-1)*pasos_por_dia+pasos_por_dia]) /
                      pasos_por_dia for d in dias_año]
        plt.plot(dias_año, medias_lat, color=colores_s[i % len(colores_s)],
                 linewidth=2, linestyle='--', label=f"S {abs(lat)}°")

    plt.axhline(0, color='white', linestyle='--', linewidth=0.8)
    plt.axhline(30, color='orange', linestyle=':', linewidth=0.8, label='Límite confort superior')
    plt.axhline(-10, color='cyan', linestyle=':', linewidth=0.8, label='Límite confort inferior')
    plt.xlabel("Día del año")
    plt.ylabel("Temperatura media (°C)")
    plt.title("Comparativo anual de latitudes — hemisferio N (—) y S (--)")
    plt.legend(fontsize=7, ncol=3)
    plt.grid(True, alpha=0.3)
    plt.figtext(0.5, -0.05, texto_parametros(), ha='center', fontsize=7, color='gray')
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(f"{nombre_archivo('comparativo_anual')}.png")
    plt.show()
    plt.close()

    # Vídeo mapa de calor anual
    print("Generando vídeo mapa de calor anual...")
    lats_v = list(range(-90, 91, 10))
    lons_v = list(range(-180, 181, 10))
    mapa_anual = np.zeros((len(lats_v), len(lons_v), 270))

    for i, lat in enumerate(lats_v):
        lat_rad = lat * PI / 180
        _, ciclo_lat = simular(lat_rad, datos, EMISIVIDAD, INERCIA_TERMICA, ALBEDO, PROFUNDIDAD_OPTICA)
        for j, lon in enumerate(lons_v):
            desfase = int(lon / 360 * 86400 / PASO_TIEMPO)
            for d in range(270):
                inicio = d * pasos_por_dia + desfase
                inicio = inicio % len(ciclo_lat)
                fin = min(inicio + pasos_por_dia, len(ciclo_lat))
                temps = ciclo_lat[inicio:fin]
                mapa_anual[i, j, d] = sum(temps) / len(temps)

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    im = ax.imshow(mapa_anual[:, :, 0], aspect='auto', origin='lower',
                   extent=[-180, 180, -90, 90], vmin=-40, vmax=60, cmap='turbo')
    plt.colorbar(im, ax=ax, label='Temperatura media (°C)')
    ax.set_xlabel("Longitud (°)", color='white')
    ax.set_ylabel("Latitud (°)", color='white')
    ax.tick_params(colors='white')
    titulo_v = ax.set_title("P3N — Día 1", color='white')
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    ax.grid(True, alpha=0.3, color='white', linewidth=0.5)

    def actualizar_anual(d):
        im.set_data(mapa_anual[:, :, d])
        titulo_v.set_text(f"P3N — Día {d + 1}")
        return [im, titulo_v]

    ani = animation.FuncAnimation(fig, actualizar_anual, frames=270, interval=100, blit=True)
    writer = animation.FFMpegWriter(fps=10, bitrate=1800)
    ani.save(f"{nombre_archivo('mapa_calor_anual')}.mp4",
             writer=writer, dpi=100, savefig_kwargs={'facecolor': 'black'})
    plt.close()
    print("Vídeo anual guardado.")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    # Tabla de estaciones
    print("\n--- ESTACIONES ---")
    print(f"{'Estación':12} | {'AV inicio':10} | {'Día inicio':10} | {'Duración':10}")
    print("-" * 52)
    for est in info_estaciones():
        print(f"{est['estacion']:12} | {est['AV_grados']:>8.2f}°  | "
              f"{est['dia_inicio']:>10.1f} | {est['duracion']:>8.1f} días")
    print("-" * 52)

    # Parámetros
    print("="*60)
    print("PARÁMETROS DE LA SIMULACIÓN")
    print("="*60)
    print(f"Masa estelar:               {S3N_MASAS_SOLARES} M☉")
    print(f"Exponente masa-temperatura: {EXP_MASA_TEMPERATURA}")
    print(f"Excentricidad orbital:      {ORBITA_EXCENTRICIDAD}")
    print(f"Emisividad:                 {EMISIVIDAD}")
    print(f"Albedo:                     {ALBEDO}")
    print(f"Profundidad óptica:         {PROFUNDIDAD_OPTICA}")
    print(f"Inercia térmica:            {INERCIA_TERMICA} SI")
    print(f"Inclinación axial:          {INCLINACION_AXIAL}°")
    print("="*60)

    datos = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)

    print("¿Qué quieres calcular?")
    print("1. Temperatura en un momento y lugar concreto")
    print("2. Visión anual global por latitudes")
    opcion = input(":")

    if opcion == "1":
        opcion_momento(datos)
    elif opcion == "2":
        opcion_anual(datos)
    else:
        print("Opción no válida")
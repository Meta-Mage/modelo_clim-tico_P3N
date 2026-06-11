from parametros import *
from orbita import *
from geometria import *
from atmosfera import *
from radiacion import *
from temperatura import *
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\sala.AULASUC-214VNRO\ffmpeg\bin\ffmpeg.exe'

def grafico_ciclo_diario(ciclo, dia, latitud):
    pasos_por_dia = int(86400 / PASO_TIEMPO)
    inicio = (dia - 1) * pasos_por_dia
    horas = [i * 15 / 60 for i in range(pasos_por_dia)]
    temps = ciclo[inicio:inicio + pasos_por_dia]

    plt.figure(figsize=(10, 5))
    plt.plot(horas, temps, color='tomato', linewidth=2)
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.xlabel("Hora del día")
    plt.ylabel("Temperatura (°C)")
    plt.title(f"Ciclo diario — Latitud {latitud}°, Día {dia}")
    plt.xticks(range(0, 25, 2))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"ciclo_diario_lat{latitud}_dia{dia}.png")
    plt.show()
    print(f"Gráfico guardado como ciclo_diario_lat{latitud}_dia{dia}.png")

if __name__ == "__main__":
    datos = precalcular_orbita(S3N_LUMINOSIDAD, INCLINACION_AXIAL_RAD, SEMIEJE_MAYOR)
    latitud = 15
    latitud_rad = latitud * PI / 180
    T_final, ciclo = simular(latitud_rad, datos, EMISIVIDAD, INERCIA_TERMICA)
    grafico_ciclo_diario(ciclo, 98, latitud)


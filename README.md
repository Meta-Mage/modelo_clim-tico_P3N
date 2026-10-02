# M3N — Modelo climático de P3N (proyecto B3N)

Modelo físico del clima de **P3N**, el planeta del proyecto de worldbuilding **B3N**, que orbita la estrella **S3N**. Calcula la temperatura en cualquier punto del planeta a lo largo del año, a partir de la órbita, la geografía y la atmósfera.

Desarrollado por Carlos (Meta-Mage) y Ozan. Forma parte del ecosistema técnico de B3N:

- **M3N**: este modelo climático.
- **C3N**: el editor de mapas (app web local), en su propio repositorio.
- **H3N**: la futura herramienta que unirá los dos.

**Versión actual: `v2.2b`** (Fase 2b, octubre de 2026).

---

## Qué hace el modelo (en una página)

- **Rejilla** de 72 × 36 celdas (5° × 5°). Cada celda es tierra o agua y tiene una altitud. La geografía la dibuja C3N.
- **Órbita kepleriana** de P3N alrededor de S3N: excentricidad, estaciones y solsticio calculado.
- **Luz** de la estrella según la hora, el día y la latitud, atenuada por la atmósfera con el reparto medido en la Tierra sin nubes.
- **Suelo** con una columna de 6 capas que conduce el calor hacia abajo; **océano** con una capa de mezcla de ~48 m.
- **Atmósfera de dos capas** (Fase 2b):
  - capa límite (el primer kilómetro);
  - troposfera;
  - efecto invernadero de dos capas;
  - calor que el aire le quita al suelo (calor sensible);
  - convección.
- **Transporte de calor** hacia los polos por el aire y el océano (difusión).
- **Resultado principal:** temperatura del **aire a 2 m**, la de un parte meteorológico. También se guarda la del suelo.

Lo que **todavía no** tiene (ver el plan de fases): hielo, evaporación, nubes, vientos, corrientes oceánicas y biomas.

---

## Instalación (Linux, probado en Pop!_OS)

Necesitas `git` y Python 3.12 o superior (Pop!_OS los trae).

```
git clone https://github.com/Meta-Mage/modelo_clim-tico_P3N.git M3N
cd M3N
python3 -m venv venv
source venv/bin/activate
pip install numpy scipy matplotlib pytest
```

- `venv` es un entorno virtual: una carpeta con las librerías solo de este proyecto. **Cada vez que abras una terminal nueva**, actívalo con `source venv/bin/activate` (verás `(venv)` al principio de la línea).
- Editor recomendado: **VSCodium** (`codium .` dentro de la carpeta).
- Opcional, solo para los vídeos del modelo antiguo de un punto: `sudo apt install ffmpeg`.

### El mapa: puente con C3N

M3N lee el mapa activo de C3N desde:

```
~/Documentos/B3N/mapa_activo/mapa_activo_m3n.json
```

C3N lo escribe solo al **activar un mapa**. Si el archivo no existe, las herramientas avisan con un error que lo explica.

---

## Cómo usarlo

Todas las herramientas usan **la misma simulación** y la guardan en caché (`outputs/cache/`). La primera vez con un mapa o unos parámetros nuevos tarda unos 6 minutos; después, segundos. Si cambias cualquier parámetro o el mapa, la caché lo detecta sola y recalcula.

| Comando | Qué da |
|---|---|
| `python analisis_global.py` | Estado de todo el planeta: medias globales, tierra/agua, perihelio/afelio, tabla por bandas de latitud y gráfico anual |
| `python consulta_punto.py LAT LON` | Mínima, media y máxima de un punto concreto a lo largo del año (sin argumentos, las pregunta) |
| `python analisis_latitudes.py` | Tabla y gráfico por latitud |
| `python mapa_calor.py` | Mapa de calor del planeta |

Los resultados se guardan en `outputs/` (no se suben a GitHub).

---

## Archivos

### Núcleo actual

| Archivo | Qué hace |
|---|---|
| `parametros.py` | Parámetros físicos: estrella, planeta, órbita. **Aquí se cambian** |
| `orbita.py` | Órbita: anomalías, distancia a S3N, estaciones |
| `geometria.py` | Geometría solar: declinación, ángulo horario, ángulo cenital |
| `rejilla.py` | La rejilla de 72 × 36 celdas |
| `fase1_geografia.py` | Tierra/agua: albedo, inercia, luz absorbida por celda |
| `fase2_inercia_multicapa.py` | Columna de suelo de varias capas (conducción) |
| `fase2_difusion.py` | Transporte horizontal de calor (difusión) |
| `fase2_combinado.py` | Modelo de la Fase 2 (sin atmósfera con cuerpo). Se usa como referencia de validación |
| `fase2b_atmosfera.py` | **Modelo actual (Fase 2b)**: atmósfera de dos capas, con 7 interruptores |
| `cache_simulacion.py` | Caché de resultados y punto de entrada de las herramientas |
| `puente_c3n.py` | Lectura del mapa activo de C3N |

### Herramientas

`analisis_global.py`, `consulta_punto.py`, `analisis_latitudes.py`, `mapa_calor.py`.

### Documentación de diseño

| Archivo | Qué contiene |
|---|---|
| `NOTAS_DISENO_FASE0.md` | Decisiones de la reescritura vectorizada |
| `DISENO_FASE2B.md` | Diseño, fuentes, calibración y validación de la Fase 2b |

### Histórico

`temperatura.py`, `main.py`, `atmosfera.py`, `radiacion.py` y los `fase0_*.py` son el **modelo original de un solo punto** y su paso a rejilla. Se conservan porque el modelo actual reutiliza parte de sus funciones y porque sirven de referencia, pero **no** incluyen la geografía ni la atmósfera actuales.

---

## Plan de fases

| Versión | Fase | Estado |
|---|---|---|
| v2.0 | Fase 0: reescritura vectorizada (NumPy) | ✅ |
| v2.1 | Fase 1: geografía real (tierra/agua, altitud) | ✅ |
| v2.2 | Fase 2: transporte horizontal + columna de suelo | ✅ |
| v2.2b | Fase 2b: atmósfera de dos capas, aire a 2 m | ✅ |
| v2.3 | Fase 3: hielo y nieve (retroalimentación hielo-albedo, calor latente) | ⏳ |
| v2.4 | Fase 4: exportación de datos estructurados | ⏳ |
| v2.5 | Fase 5: humedad, evaporación, nubes, precipitación | ⏳ |
| v2.6 / v2.7 | Fase 6 / 6b: circulación atmosférica / corrientes oceánicas | ⏳ |
| v2.8 / v2.9 | Fase 7 / 7b: biomas / biomas → albedo e inercia | ⏳ |
| v3.0 | Fase 8: atmósfera de varias capas | ⏳ |

Al terminar M3N está previsto un **barrido final de parámetros** (masa de la estrella, órbita, presión, emisividades) para fijar el mundo definitivo. Los valores actuales son provisionales.

---

## Trabajar con Git

```
git pull                               # bajar lo que haya subido el otro
git add archivo1.py archivo2.py        # añadir solo lo que has cambiado
git commit -m "Qué has cambiado"
git push
```

Mejor añadir los archivos por nombre que con `git add .`, para no subir resultados o archivos temporales por error.

---

## Contacto

- Carlos: **Meta-Mage** en GitHub.

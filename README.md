# M3N — Modelo climático de P3N (proyecto B3N)

Modelo físico del clima de **P3N**, el planeta del proyecto de worldbuilding **B3N**, que orbita la estrella **S3N**. Calcula la temperatura en cualquier punto del planeta a lo largo del año a partir de la órbita, la geografía y la atmósfera.

Desarrollado por Carlos (Meta-Mage) y Ozan. Forma parte del ecosistema técnico de B3N:

- **M3N**: este modelo climático.
- **C3N**: el editor de mapas (app web local), en su propio repositorio.
- **H3N**: el espacio de trabajo que une los dos, en su propio repositorio.

**Versión actual: `v2.4.2`** (duración del día generalizada, octubre de 2026).

---

## Qué hace el modelo (en una página)

- **Rejilla** de 72 × 36 celdas (5° × 5°). Cada celda es tierra o agua y tiene una altitud. La geografía la dibuja C3N.
- **Órbita kepleriana** de P3N alrededor de S3N: excentricidad, estaciones y solsticio calculado.
- **Luz** de la estrella según la hora, el día y la latitud, atenuada por la atmósfera. El albedo del océano depende de la altura del sol.
- **Suelo** con una columna de 6 capas que conduce el calor hacia abajo; **océano** con una capa de mezcla de unos 48 m.
- **Atmósfera de dos capas** (Fase 2b):
  - capa límite (el primer kilómetro) y troposfera;
  - efecto invernadero de dos capas;
  - calor sensible entre la superficie y el aire;
  - ajuste convectivo.
- **Transporte de calor** hacia los polos por el aire y por el océano (difusión).
- **Hielo marino** con espesor real (Fase 3): se forma, crece y se funde con energía real, y refleja más luz que el agua. La nieve en tierra es solo un diagnóstico hasta que haya precipitación.
- **Resultado principal:** temperatura del **aire a 2 m**, la de un parte meteorológico. También se guardan la de la superficie, la del aire, el hielo y el **registro hora a hora** del año (Fase 4).

Lo que **todavía no** tiene (ver el plan de fases): evaporación, humedad, nubes, lluvia, vientos, corrientes oceánicas y biomas.

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

### Los puentes con C3N y H3N

M3N lee el mapa activo de C3N y escribe el clima para H3N en dos carpetas compartidas:

```
~/Documentos/B3N/mapa_activo/mapa_activo_m3n.json     C3N → M3N  (geografía; C3N lo escribe al activar o guardar el mapa activo)
~/Documentos/B3N/clima_activo/clima_activo_m3n.json   M3N → H3N  (clima; lo escribe python exportar_clima.py)
```

Si el mapa puente no existe, las herramientas avisan con un error que lo explica.

---

## Cómo usarlo

Todas las herramientas usan **la misma simulación** y la guardan en caché (`outputs/cache/`).

- La primera vez con un mapa, unos parámetros o un código nuevos tarda **unos 20 minutos**; después, segundos.
- La caché detecta sola cualquier cambio: del mapa, de los parámetros y, desde la v2.4.1, también del **código** de la física.

| Comando | Qué da |
|---|---|
| `python analisis_global.py` | Estado de todo el planeta: medias globales, tierra y agua, perihelio y afelio, hielo marino, tabla por bandas de latitud y gráfico anual |
| `python exportar_clima.py` | Exporta el clima del mapa activo para H3N (Fase 4) |
| `python consulta_punto.py LAT LON` | Clima de un punto a lo largo del año: aire a 2 m, superficie y hielo. Sin argumentos, pregunta las coordenadas |
| `python analisis_latitudes.py` | Tabla y gráfico por latitud |
| `python mapa_calor.py` | Mapa de la temperatura media anual, con la costa y el hielo marino |
| `python -m pytest` | Pruebas automáticas (ver más abajo) |

Los resultados se guardan en `outputs/`, que **no se sube a GitHub**.

---

## Archivos

### Núcleo actual

| Archivo | Qué hace |
|---|---|
| `parametros.py` | Parámetros físicos: estrella, planeta, órbita. **Aquí se cambian** |
| `orbita.py` | Órbita: anomalías, distancia a S3N, estaciones |
| `geometria.py` | Geometría solar: declinación, ángulo horario, ángulo cenital |
| `rejilla.py` | La rejilla de 72 × 36 celdas (fila 0 = norte) |
| `fase1_geografia.py` | Tierra y agua: albedo, inercia, luz absorbida por celda |
| `fase2_inercia_multicapa.py` | Columna de suelo de varias capas (conducción) |
| `fase2_difusion.py` | Transporte horizontal de calor (difusión) |
| `fase2b_atmosfera.py` | **Modelo actual**: atmósfera de dos capas, hielo marino y registro horario, con 9 interruptores |
| `fase2_combinado.py` | Modelo de la Fase 2 (sin atmósfera con cuerpo). Se conserva como referencia de validación |
| `cache_simulacion.py` | Caché de resultados y punto de entrada de las herramientas |
| `puente_c3n.py` | Lectura del mapa activo de C3N |
| `exportar_clima.py` | Exportación del clima para H3N (formato `m3n-clima`, ver `DISENO_FASE4.md`) |

### Herramientas

`analisis_global.py`, `consulta_punto.py`, `analisis_latitudes.py`, `mapa_calor.py`, `exportar_clima.py`.

### Pruebas

| Archivo | Qué comprueba | Tiempo |
|---|---|---|
| `test_modelo.py` | Formato de exportación, calendario y estaciones, clave de la caché, puente con C3N, conservación de la energía del océano profundo | Segundos |
| `test_fase0.py` | La rejilla de la Fase 0 frente al modelo original de un punto | Medio minuto |

Las validaciones de cada fase que exigen simular un clima completo (por ejemplo, "con todo apagado, el modelo da lo mismo que la versión anterior") están documentadas, con sus resultados, en los `DISENO_FASE*.md`.

### Documentación de diseño

| Archivo | Qué contiene |
|---|---|
| `NOTAS_DISENO_FASE0.md` | Decisiones de la reescritura vectorizada |
| `DISENO_FASE2B.md` | Atmósfera de dos capas: diseño, fuentes, calibración, validación y **valores vigentes** (sección 13) |
| `DISENO_FASE3.md` | Hielo marino: física, fuentes, validación |
| `DISENO_FASE4.md` | Exportación de datos: formato `m3n-clima` y validación |

### Histórico

`temperatura.py`, `main.py`, `atmosfera.py`, `radiacion.py`, `fase0_*.py`, `fase1_altitud.py`, `fase1_mapa_mixto.py` y `fase2_inercia_estacional.py` son el **modelo original de un solo punto**, su paso a rejilla y pruebas de entonces. Se conservan porque el modelo actual reutiliza parte de sus funciones (`temperatura.py` aporta el paso de tiempo y el precálculo de la órbita) y porque sirven de referencia. **No** incluyen la geografía ni la atmósfera actuales.

---

## Plan de fases

| Versión | Fase | Estado |
|---|---|---|
| v2.0 | Fase 0: reescritura vectorizada (NumPy) | ✅ |
| v2.1 | Fase 1: geografía real (tierra/agua, altitud) | ✅ |
| v2.2 | Fase 2: transporte horizontal + columna de suelo | ✅ |
| v2.2b / v2.2c | Fase 2b: atmósfera de dos capas, aire a 2 m; auditoría y recalibración | ✅ |
| v2.3 | Fase 3: hielo marino (retroalimentación hielo-albedo, calor latente de fusión) | ✅ |
| v2.4 / v2.4.1 | Fase 4: exportación de datos estructurados; revisión general | ✅ |
| v2.4.2 | Duración del día generalizada (sin "24 h" escondidas); diseño de la Fase 5a en borrador | ✅ |
| v2.5 | Fase 5: humedad, evaporación, nubes, precipitación | ⏳ |
| v2.6 / v2.7 | Fase 6 / 6b: circulación atmosférica / corrientes oceánicas | ⏳ |
| v2.8 / v2.9 | Fase 7 / 7b: biomas / biomas → albedo e inercia | ⏳ |
| v3.0 | Fase 8: atmósfera de varias capas | ⏳ |

Al terminar M3N está previsto un **barrido final de parámetros** (masa de la estrella, órbita, presión, emisividades) para fijar el mundo definitivo. Los valores actuales son provisionales.

## Limitaciones conocidas

Están documentadas para que nadie las tome por resultados:

- **Duración del día (v2.4.2).** `ROTACION_PERIODO` (en `parametros.py`) es el día **solar** de P3N y ya se puede cambiar: la posición del sol, las capas del suelo, los registros y la exportación se adaptan solos. Condiciones: debe ser múltiplo exacto del paso de tiempo (900 s) y, para tener registro horario, de 3600 s. **Pendiente:** los coeficientes de difusión D están calibrados para un día de 24 h y dependen de la rotación (D ∝ 1/Ω², Williams y Kasting 1997); se escalarán en la Fase 5a. Los archivos históricos (`temperatura.py` en sus funciones de gráficas, `main.py`) siguen suponiendo 24 h y 270 días: no los usa el modelo actual.

- **La superficie sale más caliente que el aire, sobre todo en el océano** (unos 13 °C en el ecuador). Falta la evaporación: en la Tierra, la superficie pierde unos 82 W/m² evaporando agua (Wild et al. 2019). Llega con la Fase 5.
- **La presión no varía con la altitud** y el aire de la capa límite sobre las montañas se trata como si estuviera a nivel del mar en el infrarrojo y la convección. Pendiente.
- **El aire a 2 m se interpola sin tener en cuenta la estabilidad** (perfil neutro), aunque `DISENO_FASE2B.md` §2.7 la contemplaba. Pendiente.
- **Cada simulación en caché ocupa unos 185 MB** (`outputs/cache/`), porque guarda también el registro hora a hora. Se pueden borrar sin miedo: se recalculan cuando hagan falta.
- **El hielo polar es muy grueso** (unos 15–20 m), porque los polos de P3N casi no tienen verano y el modelo no tiene dinámica del hielo. Ver `DISENO_FASE3.md` §4.1.

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

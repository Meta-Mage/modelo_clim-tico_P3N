# M3N — Modelo climático de P3N (proyecto B3N)

Modelo físico del clima de **P3N**, el planeta del proyecto de worldbuilding **B3N**, que orbita la estrella **S3N**. Calcula la temperatura en cualquier punto del planeta a lo largo del año a partir de la órbita, la geografía y la atmósfera.

Desarrollado por Carlos (Meta-Mage) y Ozan. Forma parte del ecosistema técnico de B3N:

- **M3N**: este modelo climático.
- **C3N**: el editor de mapas (app web local), en su propio repositorio.
- **H3N**: el espacio de trabajo que une los dos, en su propio repositorio.

**Versión actual: `3.12.1`** (octubre de 2026). Versiones con [SemVer](https://semver.org/lang/es/) y fases con numeración jerárquica: la norma está en `NOMENCLATURA.md` y la historia en `CHANGELOG.md`.

Lo que se exporta es el **modelo de 2 capas** (fases 0–4; v2.4.3 con la corrección de la capa de mezcla, `DISENO_FASE_2.3.md` §14). Encima hay tres bloques programados, probados y **apagados** hasta validarlos, que no cambian el clima exportado:
- la atmósfera de 20 capas (Fase 2.3, interruptor I10);
- el ciclo del agua (Fase 5.1, I11–I15);
- el **núcleo dinámico propio** (Fase 6, I16), con la Fase 6.3 casi cerrada. Su primera climatología del modo Tierra (30 años) está en `DISENO_FASE_6.3.md` §6.14 y se valida con `validar_i16.py`.

---

## Qué hace el modelo (en una página)

- **Rejilla** de 72 × 36 celdas (5° × 5°). Cada celda es tierra o agua y tiene una altitud. La geografía la dibuja C3N.
- **Órbita kepleriana** de P3N alrededor de S3N: excentricidad, estaciones y solsticio calculado.
- **Luz** de la estrella según la hora, el día y la latitud, atenuada por la atmósfera. El albedo del océano depende de la altura del sol.
- **Suelo** con una columna de 6 capas que conduce el calor hacia abajo; **océano** con una capa de mezcla de unos 48 m.
- **Atmósfera de dos capas** (Fase 2.2):
  - capa límite (el primer kilómetro) y troposfera;
  - efecto invernadero de dos capas;
  - calor sensible entre la superficie y el aire;
  - ajuste convectivo.
- **Transporte de calor** hacia los polos por el aire y por el océano (difusión).
- **Hielo marino** con espesor real (Fase 3): se forma, crece y se funde con energía real, y refleja más luz que el agua. La nieve en tierra es solo un diagnóstico hasta que haya precipitación.
- **Resultado principal:** temperatura del **aire a 2 m**, la de un parte meteorológico. También se guardan la de la superficie, la del aire, el hielo y el **registro hora a hora** del año (Fase 4).

Lo que **todavía no** tiene el clima exportado (ver el plan de fases): evaporación, humedad, nubes, lluvia, vientos, corrientes oceánicas y biomas. La evaporación, la humedad, la lluvia (v3.1) y los vientos (núcleo dinámico de la Fase 6) están programados y probados, pero apagados hasta validarlos.

---

## Instalación (Linux, probado en Pop!_OS)

Necesitas `git` y Python 3.12 o superior (Pop!_OS los trae).

```
git clone https://github.com/Meta-Mage/modelo_clim-tico_P3N.git M3N
cd M3N
python3 -m venv venv
source venv/bin/activate
pip install numpy scipy matplotlib pytest numba global-land-mask
```

- `numba` compila el ajuste convectivo de la v3.0 (unas 20 veces más rápido) y, desde la v3.6.0, las tendencias del núcleo dinámico. Sin él funciona igual (mismo resultado), más lento. `global-land-mask` da el mapa de la Tierra del modo Tierra.

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

- La primera vez con un mapa, unos parámetros o un código nuevos tarda **unos 6 minutos** en el PC de Carlos (unos 25–30 en el entorno de la IA); después, segundos.
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
| `fase2b_atmosfera.py` | **Modelo actual**: atmósfera de dos capas, hielo marino y registro horario, con 16 interruptores (I1–I9 encendidos; I10–I16 apagados: N capas, ciclo del agua, vapor radiativo y núcleo dinámico) |
| `fase30_multicapa.py` | v3.0: atmósfera de N capas en coordenada sigma (infrarrojo gris de dos flujos, ajuste convectivo exacto, diagnóstico del transporte) |
| `fase31_agua.py` | v3.1: física del agua (saturación, convección húmeda de Betts-Miller, condensación, lluvia/nieve, albedo según la estrella, suelo según su agua) |
| `modo_tierra.py` | v3.0: mapa de la Tierra para el **modo Tierra** (`M3N_MODO=tierra`, ver `parametros.py`) |
| `fase2_combinado.py` | Modelo de la Fase 2 (sin atmósfera con cuerpo). Se conserva como referencia de validación |
| `cache_simulacion.py` | Caché de resultados y punto de entrada de las herramientas |
| `puente_c3n.py` | Lectura del mapa activo de C3N |
| `exportar_clima.py` | Exportación del clima para H3N (formato `m3n-clima`, ver `DISENO_FASE_4.md`) |
| `fase6_aguas_someras.py` | Fase 6.1: aguas someras en la esfera (rejilla C, filtro polar, leapfrog + RAW) y casos de Williamson et al. (1992) |
| `fase6_nucleo.py` | Fase 6.2: núcleo dinámico seco de 20 capas (Simmons y Burridge, semiimplícito, hiperdifusión, puntos de control). Aún sin conectar a la física |
| `fase6_nucleo_nb.py` | v3.6.0: las tendencias del núcleo compiladas con numba, idénticas bit a bit a las de numpy (~1,4 veces más rápido el paso completo) |
| `fase6_forzamientos.py` | Fase 6.2: estado de Jablonowski y Williamson (2006) y forzamiento de Held y Suarez (1994) |
| `fase6_trazadores.py` | Fase 6.3: transporte del vapor positivo y conservativo (forma de flujo, Courant zonal arbitrario) |
| `fase6_acoplamiento.py` | Fase 6.3: piezas del acoplamiento física–núcleo (viento de las caras a los centros y vuelta; diagnóstico del calor por rozamiento). Aún sin conectar |
| `fase6_superficie.py` | Fase 6.3: intercambio con la superficie con viento real (rugosidad del océano del ECMWF, coeficientes a la altura de la capa baja, albedo de Cox y Munk según el viento). Aún sin conectar |
| `fase6_capa_limite.py` | Fase 6.3: capa límite en columna (arrastre de Louis a la altura real de la capa más baja, mezcla vertical de Frierson et al. 2006 implícita, calor por rozamiento exacto y positivo). Aún sin conectar |

### Herramientas

`analisis_global.py`, `consulta_punto.py`, `analisis_latitudes.py`, `mapa_calor.py`, `exportar_clima.py`, `prueba_paso.py` (sensibilidad al paso de tiempo), `calibrar_v30.py` (v3.0: calibración de la difusión en modo Tierra y convergencia en el número de capas; se lanza con `M3N_MODO=tierra python calibrar_v30.py`), `validar_v30.py` (v3.0: validación en modo Tierra frente a observaciones, o en P3N frente a la v2.4.3; con `--v31`, la v3.1 con el ciclo del agua). `calibrar_v30.py` acepta `--d-atm`, `--d-oc` y `--v31`. `calibrar_tau_vapor.py` (v3.1, prototipo I15: calibración del infrarrojo con vapor; necesita `pip install pyrtlib`). `validar_fase6_1.py` (Fase 6.1: casos de Williamson; unos minutos, `--rapido` sin el estudio de convergencia). `clima_dinamico.py` (v3.10.0: simulación larga con el núcleo dinámico, I16: arranque caliente, equilibrio y climatología de N años; barra de progreso y punto de control; `--tierra` para el modo Tierra, `--arranque v31` como alternativa, `--dias N` solo para probar). `validar_i16.py` (v3.12.0: validación del modo Tierra con I16 frente a la Tierra real, después de `clima_dinamico.py --tierra`; ~15 min por el año extra; `--sin-ano-extra` solo la parte rápida). En `herramientas/`: `diagnostico_tierra.py` (repite el año 1 del modo Tierra con I16 y guarda los últimos subpasos si algo se rompe) y `etiquetas_canonicas.sh` (crea las etiquetas de versión canónicas, `NOMENCLATURA.md` §4.1). `held_suarez.py` (Fase 6.2–6.3: prueba de Held y Suarez con el núcleo; barra de progreso y punto de control; opciones `--tau`, `--dias`, `--filas 72` para 2,5° y `--niveles m3n` para los niveles de M3N) y `analizar_held_suarez.py` (su informe y sus figuras).

### Pruebas

| Archivo | Qué comprueba | Tiempo |
|---|---|---|
| `test_modelo.py` | Formato de exportación, calendario y estaciones, clave de la caché, puente con C3N, conservación de la energía del océano profundo | Segundos |
| `test_fase0.py` | La rejilla de la Fase 0 frente al modelo original de un punto | Medio minuto |
| `test_v30.py` | v3.0: conservación de la energía del infrarrojo y del ajuste convectivo, estabilidad tras el ajuste, cierre del diagnóstico de transporte, I10 apagado, valores del modo Tierra | Segundos |
| `test_v31.py` | v3.1: conservación de energía y agua (convección, condensación, simulación corta con todo encendido), saturación, espectro, suelo | Medio minuto |
| `test_fase6_1.py` | Fase 6.1: masa y energía exactas, casos 2 y 6 de Williamson, orden de convergencia, filtro polar, estabilidad de RAW | ~15 s |
| `test_fase6_2.py` | Fase 6.2: energía exacta, reposo sobre montañas, onda de Lamb, equilibrio de JW06, semiimplícito, hiperdifusión, puntos de control | ~10 s |
| `test_fase6_3.py` | Fase 6.3: transporte del vapor (positivo, conservativo, polos) | Segundos |
| `test_fase6_acoplamiento.py` | Fase 6.3: columna con p_s variable, presión de capa de SB81, caché de la convección, conversión caras–centros, calor por rozamiento | ~10 s |
| `test_fase6_superficie.py` | Fase 6.3: rugosidad del océano, coeficientes, albedo con viento (la primera vez construye una tabla, ~10–40 s) | Segundos |
| `test_fase6_capa_limite.py` | Fase 6.3: capa límite (Louis, Monin-Obukhov, energía, momento y vapor exactos, calor por rozamiento positivo, perfil logarítmico) | ~15 s |
| `test_fase6_i16.py` | Fase 6.3: modelo acoplado (I16): agua y energía exactas, corrector de energía, punto de control, reposo sobre montañas, criterio de equilibrio, arranque caliente, climatología, corrección de presión | ~2 min |
| `test_validar_i16.py` | v3.12.0: validación del modo Tierra (tablas de referencia, calendario, transporte, función de corriente, validación de prueba completa) | ~1–2 min |
| `test_nomenclatura.py` | La norma de `NOMENCLATURA.md`: versión única, CHANGELOG ordenado, documentos y referencias, nombres antiguos, interruptores | Segundos |

Las validaciones de cada fase que exigen simular un clima completo (por ejemplo, "con todo apagado, el modelo da lo mismo que la versión anterior") están documentadas, con sus resultados, en los `DISENO_FASE*.md`.

### Documentación de diseño

Uno por fase o etapa (`NOMENCLATURA.md` §7). Además: `NOMENCLATURA.md` (la norma de nombres y versiones) y `CHANGELOG.md` (qué trae cada versión).

| Archivo | Qué contiene |
|---|---|
| `DISENO_FASE_0.md` | Fase 0: decisiones de la reescritura vectorizada |
| `DISENO_FASE_2.2.md` | Fase 2.2: atmósfera de dos capas: diseño, fuentes, calibración, validación y **valores vigentes** (§13) |
| `DISENO_FASE_2.3.md` | Fase 2.3: atmósfera de N capas: diseño, pruebas, modo Tierra, calibración, validación y revisión del código |
| `DISENO_FASE_3.md` | Fase 3: hielo marino: física, fuentes, validación |
| `DISENO_FASE_4.md` | Fase 4: exportación, formato `m3n-clima` y validación |
| `DISENO_FASE_5.md` | Fase 5: plan del agua y las nubes (borrador de la 5.1 y división 5.1/5.2) |
| `DISENO_FASE_5.1.md` | Fase 5.1: ciclo del agua sobre las N capas: implementación y validación |
| `DISENO_FASE_6.md` | Fase 6: plan del núcleo dinámico propio, por etapas validadas |
| `DISENO_FASE_6.1.md` | Fase 6.1: aguas someras, validación con Williamson et al., límite de estabilidad de RAW |
| `DISENO_FASE_6.2.md` | Fase 6.2: núcleo seco: Simmons y Burridge, semiimplícito, pruebas, hiperdifusión, Held y Suarez |
| `DISENO_FASE_6.3.md` | Fase 6.3: acoplamiento con la física, equilibrio y climatología, modo Tierra, corrector de energía y validación (§6) |

### Histórico

`temperatura.py`, `main.py`, `atmosfera.py`, `radiacion.py`, `fase0_*.py`, `fase1_altitud.py`, `fase1_mapa_mixto.py` y `fase2_inercia_estacional.py` son el **modelo original de un solo punto**, su paso a rejilla y pruebas de entonces. Se conservan porque el modelo actual reutiliza parte de sus funciones (`temperatura.py` aporta el paso de tiempo y el precálculo de la órbita) y porque sirven de referencia. **No** incluyen la geografía ni la atmósfera actuales.

---

## Plan de fases

Numeración y reglas: `NOMENCLATURA.md` §2. Qué versión trae cada cosa: `CHANGELOG.md`.

| Fase | Contenido | Estado | Versiones |
|---|---|---|---|
| 0 | Rejilla vectorizada (72 × 36) | ✅ en uso | v2.0.0 |
| 1 | Geografía real | ✅ en uso | v2.1.0 |
| 2.1 | Transporte horizontal de calor y suelo | ✅ en uso | v2.2.0 |
| 2.2 | Atmósfera de 2 capas, aire a 2 m | ✅ en uso | v2.2.1–v2.2.2 |
| 2.3 | Atmósfera de 20 capas (antigua Fase 8, adelantada) | 🔧 apagada (I10): sin circulación, polos demasiado cálidos (`DISENO_FASE_2.3.md` §13 bis) | v3.0.0 |
| 3 | Hielo marino | ✅ en uso | v2.3.0 |
| 4 | Exportación `m3n-clima` y registro horario | ✅ en uso | v2.4.0–v2.4.1 |
| — | Día de 19,84 h, hora de P3N, paso de 992 s (parámetros, sin fase) | ✅ en uso | v2.4.2–v2.4.3 |
| 5.1 | Ciclo del agua | 🔧 apagada (I11–I15), espera a la circulación (`DISENO_FASE_5.1.md` §10 bis) | v3.1.0 |
| 6.1 | Núcleo: una capa sobre la esfera | 🔧 cerrada | v3.1.0 |
| 6.2 | Núcleo seco de 20 capas | 🔧 cerrada | v3.1.0–v3.3.0 |
| 6.3 | Núcleo + física de M3N (I16) y validación en modo Tierra | 🟡 en curso | v3.3.0–v3.12.1 |
| 6.4 | P3N con el núcleo (clima oficial; será la 4.0.0) | ⏳ | — |
| 5.2 | Nubes y radiación; después, recalibrar la estrella | ⏳ | — |
| 6.5 | Corrientes oceánicas | ⏳ | — |
| 7.1 | Clasificación del clima y biomas | ⏳ | — |
| 7.2 | Efecto del bioma en la superficie | ⏳ | — |
| 9 | Calibración final del mundo (barrido de parámetros) | ⏳ | — |

La tabla sigue el orden de trabajo acordado; los números son identificadores y no se cambian al reordenar. La Fase 8 está retirada: se absorbió en la 2.3 y, la radiación por bandas, en la 5.2.

Al terminar M3N está previsto el **barrido final de parámetros** (Fase 9: masa de la estrella, órbita, presión, emisividades) para fijar el mundo definitivo. Los valores actuales son provisionales.

## Limitaciones conocidas

Están documentadas para que nadie las tome por resultados:

- **Duración del día (v2.4.2 / v2.4.3).** `ROTACION_PERIODO` (en `parametros.py`) es el día **solar** de P3N: **19,84 h** desde la v2.4.3 (día sideral 19 h 46 min 46 s; año de 326,75 días). La **hora de P3N** es 1/24 del día (2976 s) y el registro horario y la exportación van en horas de P3N. El paso de tiempo es de **992 s** (72 por día, 3 por hora). **Prueba de sensibilidad del paso superada** (06/10/2026, PC de Carlos, mapa pizarra1: 14,451 °C con 992 s frente a 14,436 °C con 496 s; diferencia +0,015 °C de media global y 0,022 °C como máximo en una banda, con límites de 0,1 y 0,3 °C). El coeficiente de difusión de la atmósfera se escala con la rotación (D ∝ 1/Ω², Williams y Kasting 1997; factor 0,683); el del océano no. Los archivos históricos (`temperatura.py` en sus funciones de gráficas, `main.py`) siguen suponiendo 24 h y 270 días: no los usa el modelo actual.

- **La superficie sale más caliente que el aire, sobre todo en el océano** (unos 13 °C en el ecuador). Falta la evaporación: en la Tierra, la superficie pierde unos 82 W/m² evaporando agua (Wild et al. 2019). Llega con la Fase 5.1 (programada y apagada).
- **La presión no varía con la altitud** y el aire de la capa límite sobre las montañas se trata como si estuviera a nivel del mar en el infrarrojo y la convección. Pendiente.
- **El aire a 2 m se interpola sin tener en cuenta la estabilidad** (perfil neutro), aunque `DISENO_FASE_2.2.md` §2.7 la contemplaba. Pendiente.
- **Cada simulación en caché ocupa unos 185 MB** (`outputs/cache/`), porque guarda también el registro hora a hora. Se pueden borrar sin miedo: se recalculan cuando hagan falta.
- **El hielo polar es muy grueso** (unos 15–20 m), porque los polos de P3N casi no tienen verano y el modelo no tiene dinámica del hielo. Ver `DISENO_FASE_3.md` §4.1.

---

## Trabajar con Git

```
git pull                               # bajar lo que haya subido el otro
git add archivo1.py archivo2.py        # añadir solo lo que has cambiado
git commit -m "Qué has cambiado"
git push
```

Mejor añadir los archivos por nombre que con `git add .`, para no subir resultados o archivos temporales por error.

Cada versión nueva: commit `M3N X.Y.Z: <resumen>`, etiqueta anotada `vX.Y.Z` y entrada en `CHANGELOG.md` (`NOMENCLATURA.md` §4–5).

---

## Contacto

- Carlos: **Meta-Mage** en GitHub.

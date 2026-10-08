# Registro de cambios de M3N

Formato [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/), con versiones [SemVer](https://semver.org/lang/es/). Las reglas y las correspondencias con los nombres antiguos están en `NOMENCLATURA.md`. Va de la versión más reciente a la más antigua.

En cada entrada:
- **Fase** es la del plan de trabajo (NOMENCLATURA §2).
- **Clima exportado** dice si cambia el clima que M3N exporta por defecto para H3N.
- **Pruebas** es el número de pruebas automáticas que pasan; "no consta" si no está documentado.

## [3.12.1] — 2026-10-08

- Fase: ninguna (nomenclatura).
- Clima exportado: no cambia (el modelo es idéntico).
- Pruebas: 113.

### Añadido
- `NOMENCLATURA.md`: la norma de fases (numeración jerárquica sin letras), versiones (SemVer), etiquetas, commits, paquetes, interruptores y documentos, con todas las correspondencias.
- `CHANGELOG.md` (este archivo), con toda la historia.
- `test_nomenclatura.py`: comprueba que el repositorio cumple la norma.
- `herramientas/etiquetas_canonicas.sh`: crea las etiquetas canónicas como alias de las antiguas.
- `herramientas/diagnostico_tierra.py`: el diagnóstico del 07/10, antes suelto en Descargas.

### Cambiado
- Documentos de diseño renombrados a `DISENO_FASE_<id>.md` (tabla en NOMENCLATURA §7), con las referencias actualizadas en todo el repositorio.
- Fases: 2b → 2.2; la antigua 8 (atmósfera de N capas) → 2.3; 5a → 5.1; 5b → 5.2; 6b → 6.5; 7 y 7b → 7.1 y 7.2; el barrido final → 9.
- Versiones citadas en documentos y comentarios con su número canónico (v3.1-pre12 → v3.10.0, etc.).
- `VERSION` sin la "v" (formato SemVer); README actualizado (plan de fases, documentos, herramientas y pruebas).

### Corregido
- La cabecera de `fase2b_atmosfera.py` hablaba de "diez interruptores"; son dieciséis.
- El factor de rotación de D aparecía como 0,684 en tres sitios; el valor es 0,683 (0,68294).

## [3.12.0] — 2026-10-08 (antes: v3.1-pre16)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 106.

### Añadido
- `validar_i16.py`: validación del modo Tierra con I16 frente a la Tierra real (temperatura CRU 1961–1990, precipitación GPCP, radiación CERES, transporte de Trenberth y Caron, hielo NSIDC, vientos y células de Hadley del modelo).
- `referencias_tierra.py`: los valores observados, cada uno con su fuente.
- `DISENO_FASE_6.3.md` §6.14: primera climatología del modo Tierra (equilibrio en el año 30, 30 años).

### Cambiado
- La versión de las cabeceras de `clima_dinamico.py` sale de una sola constante.
- Documentación: las decisiones de Carlos que constaban como abiertas (estrella hasta las nubes; rejilla de 5° por ahora) y el README al día.

## [3.11.0] — 2026-10-07 (antes: v3.1-pre15)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 101.

### Añadido
- Corrector global de energía de la dinámica (como CAM; Lauritzen y Williamson 2019), `CORRECTOR_ENERGIA_I16`. Medido en modo Tierra: la dinámica perdía −1,31 ± 0,29 W/m² en equilibrio (§6.13).

## [3.10.2] — 2026-10-07 (antes: v3.1-pre14)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 100.

### Corregido
- Inestabilidad en la meseta antártica durante la noche polar. La corrección de presión de la hiperdifusión pasa a ser una interpolación vertical sin extrapolar, con peso ≤ 0,5 (§6.12). Sustituye a la 3.10.1.

## [3.10.1] — 2026-10-07 (antes: v3.1-pre13)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 100.

### Corregido
- Inestabilidad junto a los escalones de la Antártida: corrección de presión limitada a media capa (§6.11).

## [3.10.0] — 2026-10-07 (antes: v3.1-pre12)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 99.

### Añadido
- Equilibrio (`fase6_equilibrio.py`), arranque caliente y climatología de N años (`fase6_clima.py`).
- `clima_dinamico.py`: la cadena completa, con barra de progreso y punto de control (§6.10).

## [3.9.0] — 2026-10-07 (antes: v3.1-pre11)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 95.

### Añadido
- Modelo acoplado física–núcleo: interruptor I16 `nucleo_dinamico`, apagado (§6.9).

## [3.8.0] — 2026-10-07 (antes: v3.1-pre10)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 91.

### Añadido
- Opciones del núcleo (calor de rozamiento, corrección de presión).
- `fase6_superficie.py`: rugosidad del océano del ECMWF y albedo de Cox-Munk con viento local, sin conectar (§6.8).

## [3.7.0] — 2026-10-07 (antes: v3.1-pre9)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 85.

### Añadido
- Columna con presión en superficie variable y presión de capa de Simmons y Burridge.
- `fase6_acoplamiento.py`.

### Corregido
- Caché de la geometría de la convección (§6.7).

## [3.6.0] — 2026-10-07 (antes: v3.1-pre8)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 79.

### Añadido
- Núcleo compilado con numba (`fase6_nucleo_nb.py`): idéntico bit a bit y 1,4 veces más rápido (§6.6).

## [3.5.0] — 2026-10-07 (antes: v3.1-pre7)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 76.

### Añadido
- Capa límite en columna (`fase6_capa_limite.py`: arrastre de Louis y mezcla de Frierson), sin conectar.
- Perfil del coste del núcleo (§6.5–6.6).

## [3.4.0] — 2026-10-06/07 (antes: v3.1-pre6)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 69.

### Añadido
- `held_suarez.py --niveles m3n`. Con los niveles de M3N no hace falta capa esponja (§6.4).
- Decisiones §1.1–1.8 de la 6.3, con sus fuentes (§6).

## [3.3.0] — 2026-10-06 (antes: v3.1-pre5)

- Fase: 6.2 (cierre) y 6.3 (inicio).
- Clima exportado: no cambia.
- Pruebas: 69.

### Añadido
- Cierre de la 6.2: Held y Suarez a 5° y 2,5°.
- Borrador de la 6.3 y transporte del vapor positivo y conservativo (`fase6_trazadores.py`).
- ⚠️ El reparto exacto entre las versiones 3.1.0, 3.2.0 y 3.3.0 no consta entrega a entrega en los documentos. Se comprueba con `git show --stat v3.1.0 v3.2.0 v3.3.0`.

## [3.2.0] — 2026-10-05/06 (antes: v3.1-pre4)

- Fase: 6.2.
- Clima exportado: no cambia.
- Pruebas: no consta.

### Añadido
- `held_suarez.py --filas 72`, para la prueba a 2,5° (DISENO_FASE_6.2 §10).

## [3.1.0] — 2026-10-05 (antes: v3.1-pre3; incluye los borradores v3.0-pre2, v3.1-pre1 y v3.1-pre2, nunca etiquetados)

- Fase: 2.3 (revisión), 5.1, 6.1 y 6.2 (inicio).
- Clima exportado: **sí**, por la corrección de la capa de mezcla del océano: medía 43,3 m sin querer desde la v2.4.3 y ahora ~48,7 m fijos (DISENO_FASE_2.3 §14). ⚠️ Según los documentos, la corrección entró con la revisión del 05/10 (borrador v3.0-pre2), incluida aquí.
- Pruebas: no consta.

### Añadido
- Ciclo del agua (I11–I14) y prototipo de vapor radiativo (I15), apagados (DISENO_FASE_5.1).
- Fase 6.1 (aguas someras, Williamson et al. 1992) y núcleo seco de la 6.2.

### Corregido
- Revisión del código de la 3.0.0 (DISENO_FASE_2.3 §14): la capa de mezcla, la huella de la caché, la escritura atómica y diagnósticos.

## [3.0.0] — 2026-10-05 (antes: v3.0-pre1)

- Fase: 2.3 (antes "Fase 8", adelantada).
- Clima exportado: no cambia (I10 apagado).
- Pruebas: no consta.

### Añadido
- Atmósfera de N capas en sigma (I10), modo Tierra, `calibrar_v30.py` y `validar_v30.py` (DISENO_FASE_2.3).

## [2.4.3] — 2026-10-04

- Fase: ninguna (parámetros del mundo).
- Clima exportado: **sí**.

### Cambiado
- Día solar de 19,84 h (Carlos y Ozan); hora de P3N = 1/24 del día; paso de 992 s; D de la atmósfera escalado con la rotación (factor 0,683).

## [2.4.2] — 2026-10-04

- Fase: ninguna.
- Clima exportado: no cambia (la duración del día se generaliza sin cambiar su valor).

### Cambiado
- Sin "24 h" escondidas en el código; borrador del diseño del agua (DISENO_FASE_5).

## [2.4.1] — 2026-10-03

- Fase: 4.

### Corregido
- Revisión general: hitos del año, codificación, agua a 0 m, huella del código en la caché, océano profundo y salto del hielo (DISENO_FASE_4, DISENO_FASE_3 §2.3).

## [2.4.0] — 2026-10-03 (sin etiqueta propia: incluida en la v2.4.1)

- Fase: 4.

### Añadido
- Exportación `m3n-clima` v1 y registro horario (DISENO_FASE_4).

## [2.3.0] — 2026-10-02 (antes: v2.3)

- Fase: 3.
- Clima exportado: **sí** (15,30 °C con hielo).

### Añadido
- Hielo marino termodinámico (I9) (DISENO_FASE_3).

## [2.2.2] — 2026-10-02 (antes: v2.2c; sin etiqueta propia, incluida en la v2.3)

- Fase: 2.2.
- Clima exportado: **sí**.

### Corregido
- Auditoría:
  - umbral convectivo de 40,5 K;
  - albedo del océano según el sol (I8);
  - D_atm 2,4 y D_oc 0,12;
  - masa de S3N 0,874 (DISENO_FASE_2.2 §12).

## [2.2.1] — 2026-10-02 (antes: v2.2b)

- Fase: 2.2.
- Clima exportado: **sí**.

### Añadido
- Atmósfera de 2 capas, calor sensible, convección y aire a 2 m (I1–I7) (DISENO_FASE_2.2).

## [2.2.0] — (antes: v2.2; commit a0b5059)

- Fase: 2.1.

### Añadido
- Transporte horizontal de calor (difusión) y columna de suelo.

## [2.1.0] — (antes: v2.1)

- Fase: 1.

### Añadido
- Geografía real del mapa de C3N.

## [2.0.0] — 2026-09-20 (antes: v2.0)

- Fase: 0.

### Añadido
- Reescritura vectorizada en una rejilla de 72 × 36 a partir del modelo original de un punto (DISENO_FASE_0).

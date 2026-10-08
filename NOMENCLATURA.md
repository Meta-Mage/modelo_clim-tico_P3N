# Nomenclatura de M3N: fases, versiones, entregas y documentos

**Norma del proyecto desde la versión 3.12.1 (08/10/2026), aprobada por Carlos:** "a partir de ahora pasa a ser muy prioritario que seamos rigurosos y estrictos con esto". Las pruebas automáticas (`test_nomenclatura.py`) comprueban que el repositorio la cumple.

Aquí están las reglas, las tablas de correspondencia con los nombres anteriores y la lista de comprobación para añadir cosas nuevas. El historial de versiones está en `CHANGELOG.md`.

---

## 1. Principio: el plan y el código se numeran por separado

- Una **fase** es un bloque del **plan de trabajo**: qué física o qué capacidad se construye.
- Una **versión** es un estado del **código**, aplicado, probado y etiquetado en git.
- Una fase puede ocupar muchas versiones, y una versión puede no pertenecer a ninguna fase, como un arreglo.
- Cada versión indica en `CHANGELOG.md` a qué fase pertenece.
- El número de versión no dice nada de la fase, y el número de fase no dice nada de la versión.

---

## 2. Fases: numeración jerárquica (desglose del trabajo)

Es el esquema estándar de gestión de proyectos (WBS, *work breakdown structure*):

- **Fase N** es un área de trabajo: 0, 1, 2…
- **Etapa N.M** es una parte de esa fase: 6.1, 6.2…
- **N.M.K** solo si una etapa necesita partes.

**Reglas:**

1. **Sin letras.** Nada de "2b", "5a" ni "6b": siempre decimales.
2. **El número identifica, no ordena.** El orden de trabajo lo marca la hoja de ruta (tabla de abajo y README). Adelantar o retrasar una fase no cambia su número.
3. **Un número asignado no cambia nunca ni se reutiliza.** Una fase descartada o absorbida por otra queda en la tabla como *retirada*, con una nota.
4. **Lo nuevo:**
   - una etapa nueva dentro de una fase es el siguiente decimal libre (6.6, 6.7…);
   - un área nueva es el siguiente entero libre (10, 11…).
5. **Estados:**
   - ✅ hecha y en uso (en el clima exportado);
   - 🔧 hecha y apagada (programada y probada, detrás de un interruptor);
   - 🟡 en curso;
   - ⏳ pendiente;
   - ✖ retirada.

### 2.1 Fases de M3N

| Fase | Nombre | Estado | Versiones | Documento | Nombre anterior |
|---|---|---|---|---|---|
| 0 | Rejilla vectorizada (72 × 36, 5°) | ✅ | v2.0.0 | `DISENO_FASE_0.md` | Fase 0 |
| 1 | Geografía real (tierra/agua, altitud, albedo e inercia por celda) | ✅ | v2.1.0 | — (README) | Fase 1 |
| 2 | Atmósfera y transporte de calor | — | — | — | — |
| 2.1 | Transporte horizontal de calor (difusión) y suelo de varias capas | ✅ | v2.2.0 | — (README) | Fase 2 |
| 2.2 | Atmósfera de 2 capas, calor sensible, convección, aire a 2 m | ✅ | v2.2.1, v2.2.2 | `DISENO_FASE_2.2.md` | Fase 2b |
| 2.3 | Atmósfera de N (20) capas en sigma | 🔧 (I10) | v3.0.0 | `DISENO_FASE_2.3.md` | Fase 8 (adelantada), "v3.0" |
| 3 | Hielo marino termodinámico | ✅ | v2.3.0 | `DISENO_FASE_3.md` | Fase 3 |
| 4 | Exportación del clima (`m3n-clima`) y registro horario | ✅ | v2.4.0, v2.4.1 | `DISENO_FASE_4.md` | Fase 4 |
| 5 | Agua y nubes (plan general) | — | — | `DISENO_FASE_5.md` | Fase 5 |
| 5.1 | Ciclo del agua (vapor, lluvia, nieve, suelo, convección húmeda) | 🔧 (I11–I14; I15 prototipo) | v3.1.0 | `DISENO_FASE_5.1.md` | Fase 5a, "v3.1" |
| 5.2 | Nubes y radiación (con el vapor radiativo); después, recalibrar la estrella | ⏳ | — | — | Fase 5b (y la radiación por bandas de la antigua Fase 8) |
| 6 | Circulación (plan general del motor de vientos propio) | — | — | `DISENO_FASE_6.md` | Fase 6 |
| 6.1 | Una capa de fluido sobre la esfera (aguas someras) | 🔧 cerrada | v3.1.0 | `DISENO_FASE_6.1.md` | 6.1 |
| 6.2 | Núcleo dinámico seco de 20 capas | 🔧 cerrada | v3.1.0 – v3.3.0 | `DISENO_FASE_6.2.md` | 6.2 |
| 6.3 | Núcleo + física de M3N (I16) y validación en modo Tierra | 🟡 | v3.3.0 – v3.15.0 (salvo la v3.12.1, sin fase) | `DISENO_FASE_6.3.md` | 6.3 |
| 6.4 | P3N con el núcleo dinámico (clima oficial) | ⏳ | — | — | 6.4 |
| 6.5 | Corrientes oceánicas | ⏳ | — | — | Fase 6b |
| 6.6 | Arrastre de ondas de gravedad que conserva el momento (después de la 5.2) | ⏳ | — | — | (nueva, 09/10/2026; sustituye a la capa esponja prevista en `DISENO_FASE_6.3.md` §6.4) |
| 7 | Biomas | — | — | — | Fase 7 |
| 7.1 | Clasificación del clima y bioma de cada celda | ⏳ | — | — | Fase 7 |
| 7.2 | Efecto del bioma en la superficie (albedo, rugosidad, inercia, agua) | ⏳ | — | — | Fase 7b |
| 8 | (retirada: absorbida en la 2.3 y, la radiación por bandas, en la 5.2) | ✖ | — | — | Fase 8 "varios niveles" |
| 9 | Calibración final del mundo (barrido de parámetros) | ⏳ | — | — | "barrido final" |

**Cambios de parámetros del mundo sin fase propia** (v2.4.2 y v2.4.3: día de 19,84 h, hora de P3N, paso de 992 s): solo son versiones, registradas en `CHANGELOG.md`.

**Orden de trabajo acordado** (hoja de ruta): 6.3 → 6.4 (con ella, la 4.0.0: clima oficial provisional, sin nubes; decisión del 09/10/2026) → 5.2 (después, la estrella) → 6.6 (después de la 5.2, decisión del 09/10/2026) → 6.5 → 7.1 → 7.2 → 9. El orden entre la 6.6 y la 6.5 es una propuesta de la IA 🔶. La rejilla de 2,5° se decide cuando el planeta esté definido (decisión del 06/10/2026).

---

## 3. Versiones: Versionado Semántico (SemVer 2.0.0)

Formato **MAYOR.MENOR.ARREGLO** ([semver.org](https://semver.org/lang/es/)), adaptado a un modelo científico:

| Sube | Cuándo en M3N | Ejemplo |
|---|---|---|
| **MAYOR** | Cambia el modelo que produce el **clima oficial exportado**, de forma que los resultados dejan de ser comparables, o el formato de exportación cambia de forma incompatible | encender el núcleo dinámico para P3N → 4.0.0 |
| **MENOR** | Capacidad nueva (física, herramienta, opción), aunque esté apagada | 3.12.0 → 3.13.0 |
| **ARREGLO** | Corrección, rendimiento sin cambio de resultados, documentación o nomenclatura, sin capacidades nuevas | 3.12.0 → 3.12.1 |

**Reglas:**

1. **Toda entrega aplicada es una versión.** No hay entregas sin número.
2. **El clima exportado.** Si una versión cambia el clima exportado por defecto, `CHANGELOG.md` lo dice en su línea "¿Cambia el clima exportado?". Si el cambio es pequeño y justificado (un arreglo de física) basta con MENOR; si deja de ser comparable, MAYOR.
3. **Versiones previas:** solo como candidatas antes de una MAYOR, con la forma estándar `4.0.0-rc.1`, `4.0.0-rc.2`… Nunca "preN": en SemVer, "pre16" se ordena antes que "pre6", y en Python (PEP 440) "pre" significa *release candidate*.
4. **Dónde se escribe la versión:**
   - en el código, `VERSION` en `clima_dinamico.py`;
   - en `README.md`, la línea "Versión actual: `X.Y.Z`";
   - en `CHANGELOG.md`, la primera entrada.

   Sin la "v". Las tres coinciden siempre (lo comprueba `test_nomenclatura.py`).
5. **El formato de datos lleva su propia versión,** independiente de la del modelo: `m3n-clima` va por su versión 1 (`DISENO_FASE_4.md`).
6. **Historia anterior a la norma.** Las versiones v2.x se escribieron antes de adoptarla: sus números canónicos solo **normalizan la forma** a tres cifras, no reinterpretan su significado. Desde la v3.0.0 se aplican las reglas.
7. **Formas cortas en los documentos.** Por costumbre, «la v3.0» y «la v3.1» designan las variantes del modelo con la Fase 2.3 y con la 5.1 encendidas, y «la v2.3», «la v2.2», versiones antiguas. Se admiten como nombres de variante o de versión de la serie (3.0.x, 3.1.x…). En texto nuevo se escribe la versión completa (v3.12.1) o el número de fase.
8. **C3N y H3N** siguen la misma norma en sus repositorios (ya usan SemVer 0.y.z: C3N 0.9.5, H3N 0.2.1).

### 3.1 Correspondencia de versiones anteriores

| Etiqueta antigua (se conserva) | Versión canónica | Commit | Fase |
|---|---|---|---|
| v2.0 (nunca se etiquetó: comprobado el 08/10/2026) | v2.0.0 (sin etiqueta) | — | 0 |
| v2.1 | v2.1.0 | 92e3b9c | 1 |
| v2.2 | v2.2.0 | c67729d (a0b5059 más el arreglo del solsticio) | 2.1 |
| v2.2b | v2.2.1 | f7f7d1c | 2.2 |
| v2.2c (sin etiqueta: incluida en la v2.3) | v2.2.2 | — | 2.2 |
| v2.3 | v2.3.0 | 79d98b7 | 3 |
| v2.4 (sin etiqueta: incluida en la v2.4.1) | v2.4.0 | — | 4 |
| v2.4.1, v2.4.2, v2.4.3 | v2.4.1, v2.4.2, v2.4.3 (sin cambio) | 0165ba7, 9ada011, bdfc2ab | 4; parámetros |
| v3.0-pre1 | v3.0.0 | b879bd4 | 2.3 |
| v3.0-pre2, v3.1-pre1, v3.1-pre2 (nunca etiquetadas) | incluidas en la v3.1.0 | — | 2.3, 5.1 |
| v3.1-pre3 | v3.1.0 | fd3e84f | 5.1, 6.1, 6.2 |
| v3.1-pre4 | v3.2.0 | 5bd92ac | 6.2 |
| v3.1-pre5 | v3.3.0 | 37e8b87 | 6.2, 6.3 |
| v3.1-pre6 | v3.4.0 | 216f535 | 6.3 |
| v3.1-pre7 | v3.5.0 | 573e83f | 6.3 |
| v3.1-pre8 | v3.6.0 | bce0b5b | 6.3 |
| v3.1-pre9 | v3.7.0 | d8b4995 | 6.3 |
| v3.1-pre10 | v3.8.0 | 73aebab | 6.3 |
| v3.1-pre11 | v3.9.0 | 948aa36 | 6.3 |
| v3.1-pre12 | v3.10.0 | 17d31a6 | 6.3 |
| v3.1-pre13 | v3.10.1 | 84e6f18 | 6.3 |
| v3.1-pre14 | v3.10.2 | 5ab05f8 | 6.3 |
| v3.1-pre15 | v3.11.0 | 5275548 | 6.3 |
| v3.1-pre16 | v3.12.0 | 627d938 | 6.3 |
| — | v3.12.1 (primera con esta norma) | 1165eca | — |
| — | v3.13.0 | 6d1f720 | 6.3 |
| — | v3.13.1 | 6f45433 | 6.3 |
| — | v3.14.0 | 371ca04 | 6.3 |
| — | v3.15.0 | — | 6.3 |

"—" en Commit significa que el hash no consta en los documentos. Se obtiene con `git rev-list -n1 <etiqueta>`.

---

## 4. Etiquetas, commits y paquetes

### 4.1 Etiquetas de git

- **Formato:** etiqueta **anotada** `vX.Y.Z` en el commit de esa versión: `git tag -a v3.13.0 -m "M3N 3.13.0: <resumen>"`.
- **Inmutables:** una etiqueta publicada nunca se borra ni se mueve.
- **Historial anterior:** las etiquetas antiguas (v3.1-pre6…) se conservan. `herramientas/etiquetas_canonicas.sh` crea las canónicas como **alias** del mismo commit. Es idempotente y no toca nada si ya existen bien.

### 4.2 Commits

- Mensaje del commit de una versión: `M3N X.Y.Z: <resumen en una línea>`.
- Todo va a la rama `master`.

### 4.3 Paquetes (scripts de parche)

- **Nombre:** `NNN_M3N_X.Y.Z_desde_A.B.C.sh`.
  - `NNN` es el orden de entrega, con **tres cifras** (021, 022…), para que se ordene bien también por encima de 99.
  - `desde` es la versión que debe tener el repo; el script lo comprueba con la línea "Versión actual" del README.
- **Comprobaciones del script:**
  - termina en `OK: M3N X.Y.Z aplicada correctamente`;
  - no cambia nada si algo no cuadra;
  - se verifica byte a byte sobre una copia exacta de la base antes de entregarlo.
- **Orden de aplicación:**
  1. aplicar el script;
  2. pasar las pruebas (`python -m pytest -q`);
  3. hacer el commit;
  4. crear la etiqueta;
  5. subir (`git push && git push origin vX.Y.Z`). Las etiquetas se suben **una a una**: el 08/10/2026 GitHub rechazó en bloque, sin dar el motivo, un `git push origin --tags` con 20 etiquetas nuevas, y aceptó sin problema la v3.12.1 enviada sola.
- **Paquetes anteriores (01–20, con dos cifras y nombres de versión antiguos):**

| Paquete | Versión antigua | Versión canónica |
|---|---|---|
| 01 | v2.4.2 | v2.4.2 |
| 02–07 | v2.4.3 … v3.1-pre3 (el 04–05/10; el detalle uno a uno no consta en los documentos) | v2.4.3 … v3.1.0 |
| 08 | v3.1-pre4 | v3.2.0 |
| 09 | v3.1-pre5 | v3.3.0 |
| 10 – 20 | v3.1-pre6 – v3.1-pre16 | v3.4.0 – v3.12.0 (tabla 3.1) |
| 021 | — | v3.12.1 (esta norma) |
| 022 | — | v3.13.0 |
| 023 | — | v3.13.1 |
| 024 | — | v3.14.0 |
| 025 | — | v3.15.0 |

---

## 5. CHANGELOG.md

Formato estándar [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/), de la más reciente a la más antigua. Cada versión lleva:

- `## [X.Y.Z] — AAAA-MM-DD` (y `(antes: <etiqueta antigua>)` si la tuvo);
- **Fase:** a qué fase o fases pertenece;
- **¿Cambia el clima exportado?** sí o no (y cuánto, si se sabe);
- **Pruebas:** número de pruebas automáticas que pasan;
- **Apartados:** Añadido / Cambiado / Corregido / Retirado (solo los que tengan contenido).

---

## 6. Interruptores

- **Formato:** `I<n>` más una clave descriptiva en `INTERRUPTORES_FASE2B` (`fase2b_atmosfera.py`).
- **Inmutables:** un número no cambia ni se reutiliza.
- **Nuevos:** I17, I18… Un interruptor nuevo nace **apagado**, y con él apagado el resultado tiene que ser idéntico bit a bit al anterior (prueba obligatoria).
- Esta tabla debe coincidir con el código (lo comprueba `test_nomenclatura.py`):

| Interruptor | Clave | Qué hace | Por defecto | Fase |
|---|---|---|---|---|
| I1 | `luz_absorbida` | Reparto de la luz (Wild et al. 2019) | encendido | 2.2 |
| I2 | `capacidad_atmosfera` | Capa límite y troposfera con temperatura y capacidad propias | encendido | 2.2 |
| I3 | `calor_sensible` | Calor sensible suelo–aire (Louis, Tiedtke y Geleyn 1982) | encendido | 2.2 |
| I4 | `capa_limite_radiativa` | La capa límite absorbe y emite infrarrojo | encendido | 2.2 |
| I5 | `ajuste_convectivo` | Ajuste convectivo entre capa límite y troposfera | encendido | 2.2 |
| I6 | `difusion_reubicada` | Difusión en la troposfera y el océano, no en la piel del suelo | encendido | 2.2 |
| I7 | `altitud_en_fisica` | Altitud dentro de la física | encendido | 2.2 |
| I8 | `albedo_oceano_solar` | Albedo del agua según la altura del sol | encendido | 2.2 (v2.2.2) |
| I9 | `hielo_marino` | Hielo marino termodinámico | encendido | 3 |
| I10 | `atmosfera_multicapa` | Atmósfera de N capas en sigma | apagado | 2.3 |
| I11 | `ciclo_agua` | Vapor, evaporación, condensación, lluvia y nieve, suelo | apagado | 5.1 |
| I12 | `conveccion_humeda` | Convección húmeda (Betts-Miller simplificado) | apagado | 5.1 |
| I13 | `suelo_termico_agua` | Propiedades térmicas del suelo según su agua | apagado | 5.1 |
| I14 | `albedo_espectral` | Albedo de la nieve y del hielo con el espectro de S3N | apagado | 5.1 |
| I15 | `vapor_radiativo` | Infrarrojo según el vapor del modelo (prototipo) | apagado | 5.1 (→ 5.2) |
| I16 | `nucleo_dinamico` | Núcleo dinámico propio (viento real) | apagado | 6.3 |

---

## 7. Documentos de diseño

- **Uno por fase o etapa:** `DISENO_FASE_<id>.md` (por ejemplo, `DISENO_FASE_6.3.md`). Nunca un documento por versión.
- **Dentro de cada documento:**
  - secciones numeradas (§N, §N.M);
  - las decisiones con fecha;
  - el resultado de cada versión como una sección nueva, sin reescribir las anteriores.
- **Todos tienen que figurar** en la tabla del README.
- **Correspondencia con los nombres anteriores:**

| Ahora | Antes |
|---|---|
| `DISENO_FASE_0.md` | `NOTAS_DISENO_FASE0.md` |
| `DISENO_FASE_2.2.md` | `DISENO_FASE2B.md` |
| `DISENO_FASE_2.3.md` | `DISENO_V3.0.md` |
| `DISENO_FASE_3.md` | `DISENO_FASE3.md` |
| `DISENO_FASE_4.md` | `DISENO_FASE4.md` |
| `DISENO_FASE_5.md` | `DISENO_FASE5A.md` |
| `DISENO_FASE_5.1.md` | `DISENO_V3.1.md` |
| `DISENO_FASE_6.md` | `DISENO_FASE6_BORRADOR.md` |
| `DISENO_FASE_6.1.md` | `DISENO_FASE6_1.md` |
| `DISENO_FASE_6.2.md` | `DISENO_FASE6_2.md` |
| `DISENO_FASE_6.3.md` | `DISENO_FASE6_3.md` |

---

## 8. Lo que NO se renombra, y por qué

- **Los módulos de Python conservan su nombre:** `fase2b_atmosfera.py`, `fase30_multicapa.py`, `fase31_agua.py`, `fase6_*.py`, `INTERRUPTORES_FASE2B`…
  - Renombrarlos rompería los `import` y los puntos de control guardados, que referencian clases por su módulo (por ejemplo, `fase6_clima.Acumulador` dentro de `climatologia.pkl`).
  - Son identificadores de código, no nombres de fase.
- **Las etiquetas antiguas de git** (ver 4.1).
- **Las versiones v2.x** (ver 3, regla 6).

---

## 9. Lista de comprobación para añadir algo nuevo

- **Una etapa o fase nueva:**
  1. darle el siguiente número libre (2, regla 4);
  2. añadirla a la tabla 2.1 y a la hoja de ruta del README;
  3. crear `DISENO_FASE_<id>.md` cuando empiece el diseño.
- **Una versión nueva:**
  1. decidir MAYOR, MENOR o ARREGLO (tabla 3);
  2. actualizar `VERSION`, la línea del README y una entrada nueva arriba en `CHANGELOG.md`;
  3. crear el paquete `NNN_M3N_X.Y.Z_desde_A.B.C.sh`;
  4. hacer el commit y la etiqueta como en 4.1–4.2.
- **Un interruptor nuevo:**
  1. darle el siguiente número (I17…), apagado;
  2. hacer la prueba "apagado = idéntico bit a bit";
  3. añadirlo a la tabla 6.
- **Una herramienta nueva:**
  1. si es parte del flujo normal, va en la carpeta del repo; si es auxiliar, en `herramientas/`;
  2. añadirla al README.

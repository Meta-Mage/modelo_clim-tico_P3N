# M3N — Diseño de la Fase 2b: atmósfera con cuerpo

**Estado:** aprobado por Carlos (02/10/2026) e **implementado** en `fase2b_atmosfera.py` (02/10/2026). Ver la sección 11 para lo que cambió al implementarlo y los resultados de la validación.
**Fecha:** 01/10/2026 · **Versión prevista:** `v2.2b` · **Parte de:** `v2.2` (commit `a0b5059` + arreglo del solsticio)

---

## 0. Resumen en una página

Hoy, en M3N, el suelo es lo único que tiene temperatura. La atmósfera es un filtro sin memoria:

- tira al espacio el 40 % de la luz que atenúa (τ = 0,3);
- devuelve calor al suelo **al instante**, calculado con la temperatura del propio suelo: `(1 − ε/2)·σT⁴`.

La Fase 2b la convierte en dos capas con temperatura y capacidad propias:

| Capa | Qué representa | Capacidad (1 bar) | Responde en |
|---|---|---|---|
| **Capa límite (CL)** | primer ~0,85 km de aire (1000→900 hPa) | 1,11×10⁶ J/m²/K | horas |
| **Troposfera (TR)** | el resto de la columna (900→0 hPa) | 9,98×10⁶ J/m²/K | semanas |

Y añade cinco piezas de física:

1. **Reparto de la luz** con valores medidos de la Tierra sin nubes (Wild et al. 2019).
2. **Radiación infrarroja de dos capas grises**, calibrada con el balance terrestre sin nubes.
3. **Calor sensible** suelo↔aire con la fórmula bulk y estabilidad de Louis (1979).
4. **Ajuste convectivo** entre la CL y la troposfera (Manabe y Wetherald 1967): de día el aire caliente sube; de noche las capas se desacoplan.
5. **Difusión horizontal reubicada**: la parte atmosférica actúa sobre la troposfera y la oceánica sobre la capa de mezcla del océano (reparto guiado por Trenberth y Caron 2001).

Además:

- **Temperatura de referencia nueva:** la del **aire a 2 m**, la de un parte meteorológico. La del suelo se sigue guardando.
- **La altitud entra en la física** (hoy solo se resta al final, en los resultados).
- **Masa de la estrella recalibrada** para mantener la media global en torno a los valores actuales (decisión de worldbuilding de Carlos).

**Principio de validación:** con todos los interruptores apagados, el modelo nuevo debe reproducir **exactamente** `v2.2`.

---

## 1. Por qué hace falta (diagnóstico)

Comprobado sobre el código y los resultados de `analisis_global.py` (mapa `prueba1`, 01/10/2026):

| Síntoma | Causa | Se resuelve en |
|---|---|---|
| Tierra ecuatorial con ~35 °C de diferencia día/noche (como un desierto de roca) | Solo frenan la oscilación la radiación (~3,8 W/m²/K) y el subsuelo (~21 W/m²/K); falta el aire | 2b (calor sensible + atmósfera con capacidad) |
| La radiación nocturna hacia el suelo cae con el suelo | La atmósfera no tiene temperatura propia | 2b |
| El 40 % de la luz atenuada desaparece del balance | `i_atm()` multiplica por `exp(−τ·m)` y el resto no va a ningún sitio | 2b (reparto de la luz) |
| τ = 0,3 deja llegar al suelo solo el 60 % de la luz (medido en la Tierra sin nubes: 72,6 %) | τ demasiado alto | 2b (τ ≈ 0,18) |
| La difusión borra el contraste tierra/mar | D actúa sobre la piel del suelo (1–2 cm de roca) | 2b (difusión reubicada) |
| Una montaña solo "parece" fría en las tablas | `correccion_altitud()` se aplica después de simular | 2b (altitud en la física) |
| Océano ecuatorial a 32–35 °C | Falta la evaporación | **Fase 5** (no en 2b) |
| Océano polar líquido a −10 °C | Falta el hielo | **Fase 3** (no en 2b) |

### Qué estaba ya en el plan y qué no

| Mecanismo | ¿Estaba? | Dónde queda |
|---|---|---|
| Calor sensible | No | 2b |
| Evaporación (calor latente) | Sí, Fase 5 | Fase 5 |
| Atmósfera que almacena calor | Solo en parte (Fase 8: "varios niveles") | 2b pone 2 capas; la Fase 8 añade más |
| Luz absorbida por la atmósfera | No | 2b |

---

## 2. Componentes de la física

Notación: `T_s` suelo/superficie (capa 0 de la columna), `T_b` capa límite, `T_t` troposfera. Todas en K dentro del modelo.

### 2.1 Onda corta (luz de la estrella)

**Referencia:** Wild et al. (2019), balance global **sin nubes** (cielo despejado), media anual. De 340 W/m² en lo alto de la atmósfera:

| Flujo | W/m² | Fracción |
|---|---|---|
| Llega al suelo (directa + difusa) | 247 ± 3 | 72,6 % |
| Absorbido por la atmósfera | 73 | 21,5 % |
| Reflejado al espacio (atmósfera + suelo) | 53 ± 2 | 15,6 % |
| Absorbido por el suelo | 214 | — |
| Reflejado por el suelo | 33 | — |

**Esquema nuevo:**

1. **Transmisión hacia el suelo:** se mantiene la forma actual `exp(−τ·m)`, con la masa de aire de Pickering (2002) que ya usa M3N, pero interpretada como luz **directa + difusa** que llega al suelo. Con **τ = 0,18** (calculado: con la fórmula de masa de aire de M3N, la media sobre el planeta da 0,7265), se reproduce el 72,6 % medido.
2. **Lo atenuado en el camino de bajada** (27,4 %) se reparte así:
   - **73 % lo absorbe la atmósfera** (calienta la troposfera);
   - **27 % vuelve al espacio** (dispersión).

   Se obtiene restando a los totales de Wild la parte que corresponde a la luz reflejada por el suelo en su camino de subida (~5–6 W/m² absorbidos, ~27–28 escapan).
3. **Luz reflejada por el suelo:** al subir atraviesa la atmósfera. Por coherencia se aplica la misma transmisión y la parte absorbida va a la troposfera. Efecto pequeño (~1–2 % del total), pero cierra el balance.

**Limitación honesta:** son valores de la atmósfera terrestre. Quien absorbe la mayor parte de la luz es el vapor de agua. Asumimos para P3N una atmósfera tipo Tierra a 1 bar. En la Fase 5, con humedad variable, la absorción dependerá de ella.

**Interruptor I1 (luz absorbida):** apagado → esquema antiguo exacto (τ = 0,3 y lo atenuado se pierde).

### 2.2 Onda larga (infrarrojo): dos capas grises

Cada capa absorbe una fracción `ε` del infrarrojo que la atraviesa y emite `ε·σT⁴` hacia arriba y hacia abajo (modelo gris clásico; Pierrehumbert 2010, cap. 4).

```
Hacia arriba:
  suelo:        F0 = σ T_s⁴
  sobre la CL:  F1 = (1 − ε_b)·F0 + ε_b·σ T_b⁴
  al espacio:   OLR = (1 − ε_t)·F1 + ε_t·σ T_t⁴

Hacia abajo:
  bajo la TR:   D1 = ε_t·σ T_t⁴
  al suelo:     DLR = (1 − ε_b)·D1 + ε_b·σ T_b⁴
```

Cada capa gana lo que absorbe y pierde lo que emite por arriba y por abajo.

**Calibración de `ε_b` y `ε_t`** (no se eligen a ojo). En un "modo Tierra" (constante solar 1361 W/m², órbita y eje terrestres), deben reproducir los valores sin nubes de Wild et al. (2019):

- infrarrojo que baja al suelo **DLR ≈ 314 ± 3 W/m²**;
- infrarrojo que sale al espacio **OLR ≈ 267 ± 3 W/m²**.

**Comprobación del caso límite:** con `ε_b = 0`, `ε_t = 0,77` y la troposfera sin capacidad, su balance da `T_t⁴ = T_s⁴/2` y la pérdida neta del suelo vuelve a ser exactamente `(1 − ε/2)·σT_s⁴`. Es decir: el modelo actual es un caso particular del nuevo.

**Interruptores:** I2 (capacidad de la atmósfera; apagado → capas instantáneas) e I4 (`ε_b`; apagado → CL transparente al infrarrojo).

### 2.3 Calor sensible (suelo ↔ capa límite)

Fórmula bulk estándar (Garratt 1992):

```
H = ρ · c_p · C_H · U · (θ_s − θ_b)
```

- `ρ = p/(R·T)`: densidad del aire (~1,2 kg/m³ a 1 bar).
- `c_p = 1004 J/kg/K`.
- `U`: viento. **Provisional: 5 m/s constante** hasta la Fase 6.
- `θ`: temperatura potencial (temperatura corregida por la altura a la que está el aire). Así se compara el suelo con la CL "como si estuvieran a la misma altura": `θ_b = T_b + Γ_d·z_b`, con `Γ_d = g/c_p = 9,0 K/km` en P3N y `z_b = 0,41 km` (altura del centro de la CL).
- `C_H`: coeficiente de intercambio. Tiene dos partes:
  - **Neutro (sin efectos de estabilidad):**
    - **Océano:** `C_H ≈ 1,1×10⁻³` (Large y Pond 1982; Smith 1988).
    - **Tierra:** depende de la rugosidad del terreno `z0` mediante `C_HN = κ² / [ln(z/z0m)·ln(z/z0h)]` (κ = 0,4, z = 10 m, z0h = z0m/10). Con **suelo desnudo (z0m = 0,01 m, decisión de Carlos)** sale **2,5×10⁻³**. La vegetación (z0m mayor) llegará por bioma en la Fase 7b.
  - **Estabilidad:** funciones de **Louis (1979)**, que dependen del número de Richardson (mide si el aire de abajo está más caliente que el de arriba):
    - suelo más caliente que el aire (día): `C_H` aumenta, típicamente ×2–3;
    - suelo más frío (noche): `C_H` cae mucho. Es la inversión nocturna: el aire frío pegado al suelo apenas se mezcla.

  Es más riguroso que un factor fijo día/noche y es el esquema clásico de los modelos globales. El valor se calcula con las temperaturas del paso anterior, lo que permite resolver el sistema de forma implícita (sección 4).

**Interruptor I3 (calor sensible):** apagado → H = 0.

### 2.4 Ajuste convectivo (capa límite ↔ troposfera)

En la Tierra, de día el aire calentado por el suelo sube y reparte el calor hacia arriba. De noche la capa baja queda estable y desacoplada. El método clásico es el **ajuste convectivo** (Manabe y Strickler 1964; Manabe y Wetherald 1967):

- Si la diferencia `T_b − T_t` supera la que permite el gradiente crítico **Γ_c = 6,5 K/km** entre las alturas de las dos capas, se redistribuye calor entre ellas hasta dejarla justo en el límite, conservando la energía total.
- Si no la supera (aire estable, típico de noche), no hay intercambio convectivo, solo radiativo.

Alturas efectivas en P3N, por equilibrio hidrostático con escala de altura `H = R·T/g ≈ 8,1 km`:

- centro de la CL (950 hPa) ≈ 0,41 km;
- centro de la troposfera (550 hPa) ≈ 4,83 km;
- diferencia crítica: `ΔT_c ≈ 6,5 × 4,42 ≈ 28,7 K`.

**Interruptor I5 (ajuste convectivo):** apagado → sin intercambio convectivo.

### 2.5 Difusión horizontal reubicada

**Ahora:** D = 0,55 W/m²/K actúa sobre la capa 0 de la columna (piel del suelo y capa del océano).

**Referencia:** Trenberth y Caron (2001). La atmósfera transporta la mayor parte del calor hacia los polos (pico ~5,0 PW a 43° N). El océano solo domina entre 0° y 17° N. A 35° la atmósfera lleva el 78 % del total en el hemisferio norte y el 92 % en el sur.

**Nuevo:**

```
D_atm = (1 − f)·D   → actúa sobre T_t (troposfera)
D_oc  = f·D         → actúa sobre la capa de mezcla del océano
```

- `D = 0,55` sigue siendo el total (Williams y Kasting 1997; decisión del 28/09).
- `f` se **calibra**, no se elige: se diagnostica el transporte meridiano de cada componente y se ajusta `f` para que el océano lleve en torno al **15–25 %** del total en latitudes medias, como en la Tierra.
- **Limitación:** con D constante no se reproduce que el océano domine en el trópico profundo. Eso queda para la Fase 6b (corrientes), que sustituirá `D_oc`.
- **Sobre tierra no hay difusión de superficie.** La tierra recibe el calor transportado a través del aire (calor sensible y radiación). Esto es lo que devolverá el contraste tierra/mar.

**Interruptor I6 (difusión reubicada):** apagado → D sobre la capa 0, como ahora.

### 2.6 Altitud dentro de la física

Hoy `correccion_altitud()` (−6,5 °C/km) se aplica solo a los resultados. Propuesta para 2b:

- La temperatura del aire de referencia en una celda con altitud `z_s` se obtiene restando `Γ = 6,5 K/km × z_s`, y **esa** es la que entra en el calor sensible y en el aire a 2 m. Así, una meseta intercambia calor con aire más frío y se enfría de verdad.
- Se elimina la corrección a posteriori (se haría dos veces).
- **No** se incluye todavía la reducción de presión y de masa de aire con la altura (afectaría a la capacidad y al efecto invernadero locales). Es un refinamiento posible, anotado como pendiente.

**Altitud máxima del mapa:** 5000 m (la escala de `C3N/app.py`, `ALTITUD_MAXIMA_METROS = 5000`). Confirmado por Carlos el 02/10 como valor válido; anula el 8848 m que recogía el briefing del 20/09.

**Interruptor I7 (altitud en la física):** apagado → corrección a posteriori, como ahora.

### 2.7 Temperatura del aire a 2 m (la "temperatura del sitio")

Es la que da un parte meteorológico: aire a 1,25–2 m, a la sombra, ventilado (estándar de la OMM). **Es la que M3N dará por defecto.**

Diagnóstico: se interpola entre la superficie (`T_s`) y el aire de la CL con el perfil logarítmico de la capa superficial, con las mismas funciones de estabilidad del calor sensible. Es el método de los modelos operativos (ECMWF, IFS Part IV, cap. 3: variables diagnósticas de capa límite).

- **De día** (CL bien mezclada): el aire a 2 m queda unos 4 K por encima de `T_b`, por el gradiente adiabático entre 0,41 km y el suelo, y algo más cerca de `T_s`.
- **De noche** (inversión): el aire a 2 m se acerca a `T_s`, más frío que la CL.

**Limitación:** con una sola capa límite, la inversión nocturna está aproximada, no resuelta en vertical.

---

## 3. Parámetros

| Parámetro | Valor | Fuente / criterio | Estado |
|---|---|---|---|
| Presión en superficie | 1 bar | Decisión de Carlos (01/10) | Fijo (hasta el barrido final) |
| Gravedad P3N | 9,053 m/s² | `parametros.py` | Derivado |
| Capacidad CL / TR | 1,11×10⁶ / 9,98×10⁶ J/m²/K | `Δp/g·c_p` | Derivado |
| τ (luz que llega al suelo) | ≈ 0,18 | Wild et al. 2019 (72,6 %) | Calibrado |
| Reparto de lo atenuado | 73 % absorbido / 27 % al espacio | Wild et al. 2019 | Calibrado |
| `ε_b`, `ε_t` | por calibrar | DLR 314 y OLR 267 (Wild et al. 2019) | Se calibra |
| Viento U | 5 m/s | Provisional hasta la Fase 6 | Provisional |
| `C_H` océano (neutro) | 1,1×10⁻³ | Large y Pond 1982; Smith 1988 | Fijo |
| `z0m` tierra | 0,01 m (`C_HN` ≈ 2,5×10⁻³) | Suelo desnudo (decisión de Carlos, 02/10); Garratt 1992 | Fijo hasta la Fase 7b |
| Estabilidad | Louis 1979 | — | Fijo |
| Γ_c (ajuste convectivo) | 6,5 K/km | Manabe y Wetherald 1967 | Fijo |
| D total | 0,55 W/m²/K | Williams y Kasting 1997 | Fijo (decisión 28/09) |
| Fracción oceánica `f` | por calibrar (océano ~15–25 % en latitudes medias) | Trenberth y Caron 2001 | Se calibra |
| Masa de S3N | ~0,85 M☉ (estimación) | Calibración: media global del aire a 2 m = **15,5 °C** (decisión de Carlos, 02/10) | Se calibra al final |
| Inercia tierra | 2500 (sin cambios) | Roca densa; rango real 400–2500 | Pendiente para la Fase 7b |
| Inercia agua | 1,2×10⁶ (≈ 48 m de capa de mezcla) | Sin cambios | Fijo |

---

## 4. Método numérico

**Idea central:** la capa límite y la troposfera se añaden **encima** de la columna de suelo que ya existe (`fase2_inercia_multicapa.py`). La columna pasa a ser `[TR, CL, suelo_0, …, suelo_N]`.

- **Implícito** (estable, sistema tridiagonal ya existente, ampliado): conducción en el suelo, calor sensible suelo↔CL con el `C_H` del paso anterior, y difusión horizontal (troposfera y océano).
- **Explícito** (como ahora): onda corta y onda larga. La onda larga es no lineal (T⁴) y se calcula con las temperaturas del paso.
- **Ajuste convectivo:** se aplica al final de cada paso. Es instantáneo y barato, y conserva la energía exactamente.
- **Paso de tiempo:** 900 s, sin cambios. Comprobación prevista de que la CL (que responde en horas) es estable con ese paso: con capacidad 1,1×10⁶ y acoplamiento radiativo de ~5 W/m²/K, el tiempo propio es de ~2,5 días, así que no hay problema. El calor sensible va en la parte implícita.

---

## 5. Optimización

Objetivo: que una simulación completa de 2b **no tarde más** que la de `v2.2` (~8 min con `prueba1`).

1. **Medir primero.** Perfilar la simulación actual (`cProfile`) para ver dónde se va el tiempo. No se optimiza nada sin medir.
2. **Arranque cercano al equilibrio.** Hoy se tarda ~16 años en converger, sobre todo por el océano. Opciones, por orden de preferencia y a decidir según las medidas:
   - partir del último equilibrio guardado en la caché con parámetros parecidos;
   - extrapolar la deriva anual del océano tras los primeros años (acelerar la convergencia sin cambiar el resultado final).
3. **Mismo sistema tridiagonal:** las dos capas nuevas añaden 2 incógnitas por columna a un sistema que ya se resuelve vectorizado. Coste previsto: pequeño.
4. **Caché:** se mantiene, con **todos** los parámetros nuevos en la huella (lección del solsticio).
5. **Último recurso:** compilar el bucle interno (por ejemplo, con `numba`). Solo si las medidas lo justifican, porque añade una dependencia.

---

## 6. Validación

| Prueba | Criterio |
|---|---|
| **V0. Caso límite** | Todos los interruptores apagados (I1–I7) → mismo resultado que `v2.2`, celda a celda (diferencia máxima < 1×10⁻⁶ K) |
| **V1. Energía en lo alto de la atmósfera** | Media anual global: luz absorbida (suelo + atmósfera) = infrarrojo emitido al espacio, con diferencia < 0,1 %. Sustituye al test actual, que solo mira la superficie |
| **V2. Energía por capa** | Cada capa cierra su balance anual (lo que entra = lo que sale) |
| **V3. Calibración modo Tierra** | Luz al suelo, absorbida y reflejada ≈ 247/73/53; DLR ≈ 314; OLR ≈ 267 (Wild et al. 2019) |
| **V4. Comprobaciones físicas en P3N** | Oscilación día/noche del suelo en tierra ecuatorial: de ~35 °C a ~15–20 °C; la radiación infrarroja nocturna hacia el suelo, casi igual que la diurna; reaparece el contraste tierra/mar; ajuste convectivo activo de día y casi inactivo de noche |
| **V5. Sensibilidad al paso de tiempo** | 900 s frente a 450 s: diferencias pequeñas (como en las Fases 0–2) |

---

## 7. Resultados que cambian

- **Media global:** sin recalibrar la estrella, subiría de ~16 °C a ~40–44 °C (estimación con un modelo simple), porque deja de perderse luz. La masa de S3N se calibra para que la media global del aire a 2 m sea **15,5 °C**. Estimación: ~0,85 M☉ (la estrella pasa de ~5.700 K a ~5.350 K, algo más anaranjada); el valor exacto se fija simulando.
- **Todas las simulaciones anteriores dejan de ser comparables.** La caché las invalida sola.
- **Herramientas** (`analisis_global.py`, `consulta_punto.py`, `analisis_latitudes.py`, `mapa_calor.py`): se adaptan. Por defecto muestran el **aire a 2 m** y además la temperatura del suelo. El registro guarda la mínima, media y máxima diarias de ambas, y la media diaria de la CL y la troposfera.
- **Puente C3N → M3N:** sin cambios.

---

## 8. Conexiones con otras fases

- **Fase 3 (hielo):** usará mínimas realistas. El hielo marino irá con la formulación por entalpía de Wagner y Eisenman (2015) y la nieve en tierra con la temperatura del suelo. Los problemas de la conversación perdida (todo el continente congelado, −225 °C sin difusión, 8,9 % de desequilibrio) se reexaminan con V1.
- **Fase 5:** la evaporación será otro flujo como H (misma fórmula bulk con humedad). Las nubes cambiarán el reparto de la luz de 2.1. El vapor de agua hará que `ε` dependa de la humedad.
- **Fase 6:** sustituye el viento constante. **Fase 6b:** sustituye `D_oc`.
- **Fase 7 (biomas):** Köppen usa la temperatura del aire a 2 m, que existirá desde 2b. **Fase 7b:** rugosidad e inercia por bioma.
- **Fase 8:** amplía las dos capas a varias.
- **Barrido final de parámetros** (propuesta de Carlos, 01/10): al terminar M3N, nuevo barrido de masa estelar, órbita, presión y ε para fijar el mundo definitivo.

---

## 9. Plan de trabajo

Cada pieza lleva su interruptor y su validación antes de pasar a la siguiente.

1. **Perfilado** de `v2.2` (medir tiempos).
2. **Esqueleto:** configuración de interruptores y columna ampliada. Validar V0 con todo apagado.
3. **Onda corta** (I1).
4. **Dos capas con capacidad + onda larga** (I2, I4) y calibración de `ε_b` y `ε_t` en modo Tierra (V3).
5. **Calor sensible con Louis** (I3).
6. **Ajuste convectivo** (I5).
7. **Difusión reubicada** (I6) y calibración de `f`.
8. **Altitud en la física** (I7).
9. **Aire a 2 m**, registro y adaptación de herramientas.
10. **Test de energía nuevo** (V1, V2) y sensibilidad (V5).
11. **Optimización** según el perfilado.
12. **Calibración de la masa de S3N.**
13. **Comparación antes/después** con `analisis_global.py`, documentación y tag `v2.2b`.

---

## 10. Decisiones

**Tomadas por Carlos (01/10/2026):**

- Fase 2b antes del hielo.
- Presión 1 bar.
- Dos capas (CL + troposfera) en lugar de una capa con desfase fijo.
- Viento constante provisional y estabilidad día/noche (concretada como Louis 1979).
- Difusión: lo más riguroso (reparto atmósfera/océano calibrado).
- Bajar la luminosidad mediante la masa de S3N, con objetivo de media global del aire a 2 m = 15,5 °C (02/10).
- Rugosidad de la tierra: suelo desnudo, z0m = 0,01 m; la vegetación llegará por bioma (02/10).
- Altitud: basta con la propuesta de 2.6 por ahora (02/10).
- Altitud máxima del mapa: 5000 m, no 8848 m (02/10).
- Barrido final de parámetros al terminar M3N.
- Versión `v2.2b`.
- Se permite ejecutar simulaciones y tests para validar.

**Pendientes:** ninguna para empezar a programar.

---

## 11. Implementación y validación (02/10/2026)

Código: `fase2b_atmosfera.py` (física, interruptores) y `cache_simulacion.py` (`simular_fase2b_cacheada`, `simular_modelo_cacheado`). Las cuatro herramientas usan ya la Fase 2b; sus tablas dan la **temperatura del aire a 2 m**.

### Cambios respecto al diseño (decididos al implementar)

1. **Calibración del infrarrojo sin "modo Tierra".** Se descartó simular una Tierra sin nubes para calibrar: en equilibrio, una Tierra sin nubes emitiría al espacio lo que absorbe (~287 W/m²), no los 267 medidos. Los 267 corresponden a la atmósfera real, calentada también por las nubes. Se calibró directamente con las temperaturas medias observadas (Atmósfera Estándar US 1976, promediada en masa en cada capa): suelo 289 K, CL 285,3 K, TR 246,3 K. Con ellas, infrarrojo hacia el suelo = 314 y hacia el espacio = 267 dan **ε_b = 0,740** y **ε_t = 0,661**. Absorción total del infrarrojo: 0,912 (antes 0,77 en una capa, calibrada para la Tierra con nubes).
2. **Aire a 2 m:** interpolación logarítmica **neutra** entre el suelo y el aire de la CL (fracción ln(2/z0h)/ln(10/z0h): 0,83 en tierra, 0,88 en océano). No se aplican aún las funciones de estabilidad en la interpolación. Efecto esperado: de noche, el aire a 2 m sale algo más templado que el real, con la inversión infravalorada. Refinamiento pendiente.
3. **Calor sensible:** se resuelve como un sistema implícito de 2×2 (suelo y CL) por celda en cada paso, en lugar de dentro del sistema tridiagonal de la columna. Conserva la energía exactamente y permite que `C_H` cambie en cada paso sin refactorizar nada.
4. **Aceleración de la convergencia:** lo último en converger es un único modo lento, el océano polar, que de noche casi no intercambia calor con el aire estable (cada año cambia ~61 % de lo que cambió el anterior). Cuando ese ritmo se mantiene tres años seguidos, se salta lo que falta (suma de la serie geométrica). Después se siguen simulando años normales con el criterio de siempre. Se probó también arrancar con el océano de capacidad reducida, y se descartó: sesgaba el estado medio.
5. **Fracción oceánica del transporte:** `f` = **0,13**. Con ese valor, el océano lleva el 17 % (norte) y el 13 % (sur) del transporte a 35° (Trenberth y Caron 2001: 22 % y 8 %).
6. **Masa de S3N:** **0,884 M☉** (antes 0,97). La luz que llega a P3N baja de 1875 a 1441 W/m², y la estrella pasa de ~5700 K a ~5450 K. Calibrada simulando: 0,85 → 7,9 °C; 0,885 → 15,7 °C; **0,884 → 15,45 °C**.

### Validación

| Prueba | Resultado |
|---|---|
| V0. Interruptores apagados = `v2.2` | **Idéntico bit a bit** (diferencia 0,0 en mínimas, medias y máximas) |
| V1. Energía en lo alto de la atmósfera | Diferencia 3,8×10⁻⁵ (criterio < 10⁻³) |
| Aceleración frente a simulación sin acelerar | Diferencia máxima 0,02 °C, media 0,005 °C (dentro de la tolerancia de convergencia) |
| V5. Paso de 450 s frente a 900 s | Aire a 2 m: máx. 0,17 °C, media 0,02 °C; suelo (máximas): máx. 0,84 °C, media 0,12 °C |
| Combinaciones de interruptores | Todas finitas y estables; las combinaciones no permitidas se rechazan con un error explicativo |

### Resultados en `prueba1` (masa 0,884)

- **Media global del aire a 2 m: 15,5 °C.** Tierra 23,2 °C, océano 12,9 °C. Ciclo anual global de 14,4 a 17,0 °C (excentricidad).
- **Tierra ecuatorial:** aire a 2 m con mínima 23,4, media 27,2 y máxima 31,4 °C. **Oscilación día/noche del aire ~8 °C** (realista para los trópicos) y del suelo ~20 °C (suelo desnudo, sin evaporación: típico de zonas semiáridas). Antes, 35 °C en el suelo y nada que separara el aire del suelo.
- **Océano ecuatorial:** ~33 °C (sigue alto: falta la evaporación, Fase 5). **Océano polar:** −28 a −35 °C y líquido (falta el hielo, Fase 3).
- **Troposfera:** de ~−2 °C (ecuador) a ~−35 °C (polos), media en masa.
- **Tiempo:** ~6 min por simulación completa (11 años con 2 saltos, ~32 s por año), frente a ~8,5 min de `v2.2` antes de optimizar y ~2 min después. La Fase 2b cuesta más por año (dos capas de aire, dos sistemas de difusión).

### Pendientes conocidos

- Estabilidad en la interpolación del aire a 2 m (punto 2).
- Presión y masa de aire menores con la altitud (sección 2.6).
- Evaporación (Fase 5), que bajará el océano tropical y la oscilación del suelo.
- Barrido final de parámetros al terminar M3N (masa, órbita, presión, emisividades).

---

## Referencias

- Charnock, H. (1955). Wind stress on a water surface. *Q. J. R. Meteorol. Soc.* 81.
- ECMWF (2023). *IFS Documentation CY48R1, Part IV: Physical Processes*, cap. 3.
- Garratt, J. R. (1992). *The Atmospheric Boundary Layer*. Cambridge University Press.
- Large, W. G. y Pond, S. (1982). Sensible and latent heat flux measurements over the ocean. *J. Phys. Oceanogr.* 12.
- Louis, J.-F. (1979). A parametric model of vertical eddy fluxes in the atmosphere. *Boundary-Layer Meteorol.* 17, 187–202.
- Manabe, S. y Strickler, R. F. (1964). Thermal equilibrium of the atmosphere with a convective adjustment. *J. Atmos. Sci.* 21.
- Manabe, S. y Wetherald, R. T. (1967). Thermal equilibrium of the atmosphere with a given distribution of relative humidity. *J. Atmos. Sci.* 24.
- Pickering, K. A. (2002). The southern limits of the ancient star catalog. *DIO* 12 (fórmula de masa de aire ya usada en M3N).
- Pierrehumbert, R. T. (2010). *Principles of Planetary Climate*. Cambridge University Press.
- Smith, S. D. (1988). Coefficients for sea surface wind stress, heat flux, and wind profiles. *J. Geophys. Res.* 93.
- Trenberth, K. E. y Caron, J. M. (2001). Estimates of meridional atmosphere and ocean heat transports. *J. Climate* 14, 3433–3443.
- Wagner, T. J. W. y Eisenman, I. (2015). How climate model complexity influences sea ice stability. *J. Climate* 28.
- Wild, M. et al. (2019). The cloud-free global energy balance and inferred cloud radiative effects. *Clim. Dyn.* 52.
- Williams, D. M. y Kasting, J. F. (1997). Habitable planets with high obliquities. *Icarus* 129.

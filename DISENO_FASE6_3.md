# Fase 6.3 de M3N: acoplar el núcleo dinámico con la física de M3N (BORRADOR v1)

**Estado: v2 (06–07/10/2026).** Las decisiones del §1 están revisadas con fuentes y **aprobadas por Carlos** (§6). Hechos: el transporte del vapor (§5), el paso 2 del §3 (núcleo con los niveles de M3N: **sin capa esponja**, §6.4) y la capa límite en columna (paso 4, §6.5). La 6.2 está cerrada (`DISENO_FASE6_2.md`).

Etiquetas:
- ✅ verificado;
- ⚠️ pendiente;
- 🔶 propuesta de la IA, que decide Carlos.

---

## 0. Qué cambia

Ahora mismo, en `paso_agua` (v3.1), el aire se mueve con una **difusión horizontal** de calor (X = T + gz/c_p) y de vapor (`fact_atm30`, `fact_vapor`). En la 6.3:
- **Esa difusión desaparece de la atmósfera:** la sustituye `NucleoSeco`, que da a cada paso el viento, la presión en superficie y las temperaturas transportadas.
- **El océano mantiene su difusión** hasta la Fase 6b.
- **El vapor** lo transportan los flujos de masa del núcleo, con un esquema que no da valores negativos (§3).
- **La física de columna sigue igual y en el mismo orden:** radiación, superficie, calor sensible, evaporación, convección, condensación, lluvia y nieve, suelo e hielo.

## 1. Cuestiones de diseño (cada una, con su recomendación)

### 1.1 Reparto del paso de tiempo 🔶

- **La física** sigue con su paso: 900 s en la Tierra y 992 s en P3N.
- **La dinámica** da 2 pasos por cada paso de la física (450 s y 496 s; DISENO_FASE6_2 §4.4).
- **Acoplamiento "en serie"** (process splitting): primero la dinámica con los 2 pasos, después la física sobre el estado resultante, con su tendencia aplicada como ajuste. Es lo que hacen la mayoría de los modelos ⚠️ (por citar).
- **Detalle con el leapfrog:** la física modifica el estado actual y el filtro RAW sigue actuando. Hay que comprobar que no aparece un modo computacional. La alternativa es aplicar la física como tendencia sobre el salto de 2Δt.

### 1.2 Presión en superficie variable

- **Hoy:** `Columna` fija p_s por la altitud (`presion_superficie`).
- **Con el núcleo:** p_s evoluciona, así que `COL.ps`, `ph`, `pm` y `dp` tienen que **actualizarse en cada paso** con la p_s del núcleo.
- **Coste:** recalcular las geometrías de la convección (`_geometria_conveccion` tiene una caché que habrá que invalidar o recalcular) ⚠️.
- **Discrepancia de definiciones a resolver:** `Columna` usa la presión del centro de la capa como la media aritmética de sus seminiveles (`pm = ½(ph_k + ph_{k+1})`), mientras que el núcleo (SB81) usa ln p_k = ln p_{k+1/2} − α_k. En la capa de arriba difieren mucho: σ_k = σ_{1/2}·e^{−ln2}·…
  - 🔶 **Unificar en la definición de SB81**, la del núcleo. Así la temperatura potencial y la saturación se calculan en la misma presión que usa la dinámica.
  - Esto cambia ligeramente la física de la v3.1, y hay que medir cuánto.

### 1.3 Viento en superficie 🔶

- **Hoy:** VIENTO = 5 m/s fijo ("provisional hasta la Fase 6", línea 126).
- **Propuesta:** |v| de la capa más baja, interpolado al centro de la celda, más una ráfaga mínima: √(|v|² + v_r²).
  - La ráfaga mínima es una práctica habitual para que el flujo no se anule con el aire en calma ⚠️ (por citar).
  - Su valor (¿1 m/s?) habrá que justificarlo, porque es un parámetro.
- **Efecto:** los flujos de calor sensible, evaporación y sublimación pasan a depender del viento real. **Esto debería producir los alisios y los desiertos subtropicales** (más evaporación donde sopla más).

### 1.4 Rozamiento con la superficie y capa límite

- **Lo que falta:** M3N no tiene rozamiento sobre el viento, porque hasta ahora no había viento. Held y Suarez usa una relajación de Rayleigh idealizada, que **no** sirve para un modelo realista.
- 🔶 **Propuesta:**
  - **arrastre en superficie** de fórmula aerodinámica, τ = ρ C_D |v| v, con C_D sacado del mismo C_H neutro y del factor de estabilidad de Louis que ya usa M3N (`coeficientes_neutros`, `factor_estabilidad_louis`, que ya da la función para el momento ⚠️ comprobar);
  - **difusión vertical** del momento, del calor y del vapor en la capa límite, con una K que dependa de la estabilidad, implícita en cada columna;
  - versión mínima: el esquema de Frierson et al. (2006) para modelos de radiación gris ⚠️ (por leer en la fuente).
- **Hoy, el calor sensible entra solo en la capa más baja.** Con capas que miden decenas de metros, hace falta esa difusión vertical o se acumulará todo abajo; el ajuste convectivo seco actual lo compensa en parte.

### 1.5 Transporte del vapor: positivo y conservativo 🔶

- **El problema:** el caso 1 de la 6.1 demostró que el esquema centrado da valores negativos con campos poco suaves, y el vapor no lo es.
- **Propuesta:** volúmenes finitos en forma de flujo, con los **mismos flujos de masa del núcleo** (F_u, F_v, W, promediados sobre los 2 pasos de la dinámica para que el vapor se transporte con la misma masa que se movió) y un limitador de flujo (van Leer o similar ⚠️).
  - Esto garantiza que la masa de agua cierra exactamente y que q ≥ 0.
- **Polos:** el filtro polar no puede actuar sobre q, porque rompería la positividad. Alternativas: subpasos zonales o un esquema semilagrangiano en dirección zonal cerca de los polos ⚠️ (por estudiar). Es el punto técnico más delicado de la 6.3.
- **Coherencia:** la ecuación de la temperatura del núcleo usa la masa Δp, y el vapor tiene que usar exactamente la misma, para que q = 1 se conserve como q = 1.

### 1.6 Capa esponja para el techo de M3N

- **Por qué hace falta:** la capa de arriba de M3N va de 0 a 7,35 hPa, con el centro a 3,7 hPa (corregido el 06/10: los documentos decían "0,7 hPa" por un desliz de unidades, 0,735 kPa). La 6.2 (intento 1 de Held y Suarez) demostró que sin disipación en el tope los vientos de las capas altas crecen sin límite.
- 🔶 **Propuesta:** relajación de Rayleigh de u y v hacia 0 (o hacia su media zonal) en las capas con p < ~1 hPa, con un tiempo de amortiguamiento de ~1 día. Hay que fijar el criterio en una fuente: Polvani y Kushner (2002) ⚠️, o la práctica de los modelos con el techo alto.
- **Energía:** la energía cinética que quita se devuelve como calor, para no perder energía ⚠️ (decisión).
- **Alternativa:** reducir el techo del modelo, con otros niveles. Cambiaría la radiación de la v3.0, y no se recomienda.

### 1.7 Difusión horizontal

- La hiperdifusión de la 6.2 se mantiene.
- **Sobre montañas, la difusión de T en superficies σ es incorrecta:** crea calor y frío falsos en las laderas. Hay que difundir T sobre superficies de presión, o corregir con el gradiente de referencia ⚠️ (por diseñar). Con el relieve de la Tierra a 5° el efecto puede notarse.

### 1.8 Energía: el cierre que ya tenemos

- La v3.1 cierra la energía a 10⁻¹² (`energia_total`). Con el núcleo, la energía total pasa a incluir la **energía cinética** y el **Φ_s p_s**.
- **Objetivo:** mantener el diagnóstico de cierre con todas las pérdidas declaradas: la hiperdifusión, el filtro polar, RAW y la esponja, si no devuelve su calor.

## 2. Validación (modo Tierra, decisión ya tomada para las fases anteriores)

Referencias observadas a verificar una a una:
- **viento zonal medio por latitud y altura** (ERA5 o NCEP ⚠️): chorro subtropical de ~30–40 m/s en invierno a ~200 hPa;
- **células de Hadley y Ferrel** (función de corriente meridional: ~10¹¹ kg/s ⚠️);
- **transporte meridional de calor** (Trenberth y Caron 2001 ✅, existencia): ~5–6 PW hacia 35°;
- **humedad relativa en la troposfera libre** (seca en los subtrópicos);
- **temperatura de los polos y hielo marino;**
- **tropopausa** (más alta en el ecuador);
- **precipitación por latitud** (GPCP ⚠️).

## 3. Orden de trabajo propuesto

1. **Cerrar la 6.2:** Held y Suarez en marcha y su análisis.
2. **Núcleo con los niveles de M3N y la capa esponja (1.6),** probado con Held y Suarez en esos niveles.
3. **Transporte positivo de trazadores (1.5),** con su prueba: el caso 1 de la 6.1 ahora sin negativos, y la conservación exacta.
4. **Capa límite y arrastre (1.4),** probados en columna.
5. **Acoplamiento (1.1–1.3, 1.7, 1.8)** con el interruptor I16 y una prueba de cierre de energía y de agua.
6. **Modo Tierra:** simulación de varios años y validación (§2).
7. **P3N** (6.4).

## 4. Decisiones de Carlos

**Tomadas el 06/10/2026:** ver el §6, que sustituye a la lista de este apartado (la de la v1 era: 1.2 presión de SB81; 1.3 viento con ráfaga mínima; 1.4 arrastre y difusión vertical; 1.6 capa esponja; el orden del §3). La ráfaga mínima de 1.3 y la capa esponja de 1.6 cambian respecto a la v1; ver allí por qué.

---

## 5. Avance (05/10/2026): transporte de trazadores (§1.5) — prototipo horizontal hecho

**Código:** `fase6_trazadores.py` y `test_fase6_3.py` (3 pruebas).

### Esquema

- **Forma de flujo:** con la masa de aire de cada celda (m) y la masa que cruza cada cara en el paso (M).
- **Separación por direcciones consistente con la masa:** cada barrido actualiza a la vez m·q y m, y q = m·q/m. Así un campo constante sigue constante aunque cada barrido sea divergente (idea de Easter 1993 ⚠️, cita por verificar). El orden se alterna en cada paso (XY, YX).
- **Reconstrucción lineal** con pendiente limitada (MC), que no crea máximos ni mínimos nuevos.
- **Dirección zonal con Courant arbitrario:** se suman las celdas enteras que cruzan la cara más la fracción de la siguiente (idea de Lin y Rood 1996 ✅, la cita). Esto resuelve las celdas estrechas junto al polo **sin el filtro polar**, que rompería la positividad.

### Prueba con el caso 1 de Williamson

- Campana, 12 días, dt = 3600 s.
- Viento derivado de la función de corriente en las esquinas, así que es exactamente no divergente.

| Eje | Rejilla | l1 | l2 | l∞ | Masa | q mínimo | Courant zonal máx. |
|---|---|---|---|---|---|---|---|
| α = 0 | 5° | 0,691 | 0,502 | 0,463 | 2·10⁻¹⁶ | 0 | 0,25 |
| α = 0 | 2,5° | 0,263 | 0,218 | 0,223 | 4·10⁻¹⁶ | 0 | 0,25 |
| α = π/2 − 0,05 | 5° | 0,732 | 0,531 | 0,553 | 0 | −5·10⁻¹⁰⁹ | 5,7 |
| α = π/2 − 0,05 | 2,5° | 0,254 | 0,216 | 0,250 | 2·10⁻¹⁶ | ~0 | 11,5 |
| α = π/2 | 5° | 0,730 | 0,530 | 0,545 | 2·10⁻¹⁶ | −5·10⁻¹¹⁰ | 5,7 |
| α = π/2 | 2,5° | 0,257 | 0,218 | 0,247 | 0 | ~0 | 11,5 |

- **Frente al esquema centrado** (DISENO_FASE6_1 §4): a 5° pasa de l2 ≈ 1 y mínimos de −450 m a **l2 = 0,50 sin ningún negativo**. Los −10⁻¹⁰⁹ son redondeo en números desnormalizados, es decir, 0.
- **Por encima del polo** (Courant zonal de hasta 11) la precisión es la misma que en el ecuador: el polo no degrada nada.
- **Convergencia** de 5° a 2,5°: ÷2,3. Es el orden ~1,2, típico de los esquemas limitados a resolución gruesa ⚠️ (falta una comparación publicada).

### Pendiente

- Barrido vertical (W, Courant < 1).
- Acoplar con los flujos de masa reales del núcleo, promediados sobre los 2 subpasos.
- Decidir cómo casar la masa de aire del trazador con la p_s del núcleo, que se ve alterada por el filtro RAW y el filtro polar:
  - opción A: el trazador lleva su propia masa de aire y se mide su deriva respecto a p_s;
  - opción B: corrección proporcional ("fijador de masa") ⚠️.

### Transporte 3D acoplado al núcleo (05/10)

- **Barrido vertical:** con W del núcleo, misma reconstrucción limitada; comprueba que |Courant| < 1.
- **`transportar_3d`:** orden simétrico alterno, X Y Z / Z Y X.
- **Flujos de masa:** los del núcleo en el instante n de cada paso leapfrog (F_u, F_v, W; `guardar_flujos`), multiplicados por Δt. El trazador lleva **su propia masa de aire**.

**Prueba: onda baroclínica de JW06 durante 10 días** (dt = 450 s, con hiperdifusión), con dos trazadores: q₁ = 1 y un "vapor" q₂ = 0,02 cos⁴φ σ³.

| Diagnóstico | Resultado |
|---|---|
| q₁ | **1,000000 exacto** en todo el globo |
| Masa de q₂ | conservada a 2·10⁻¹⁶ |
| q₂ mínimo | 3,6·10⁻¹⁵ (sin negativos) |
| Diferencia entre la masa de aire del trazador y la del núcleo | **2,2·10⁻⁴, que no crece con el tiempo** (días 1 a 9) |

Esa diferencia acotada es el efecto del filtro RAW, el filtro polar y la corrección semiimplícita sobre p_s.

**Decisión 🔶 (opción A + conciliación):** en cada paso de la física, q pasa a la masa dinámica conservando la masa de agua, q ← q·m_tr/Δp_dyn, y m_tr ← Δp_dyn.
- **El agua se conserva exacta,** que es lo prioritario para el cierre del ciclo del agua.
- **Un campo constante cambia ≤ 2·10⁻⁴ relativo,** acotado y medido.

---

## 6. Decisiones revisadas con fuentes y aprobadas por Carlos (06/10/2026)

Criterio de Carlos para todo este apartado: "que todo lo que se haga se haga desde la información y la seguridad de que es apropiado, preciso y riguroso al máximo". Las decisiones técnicas de la 6.3 las verifica la IA en las fuentes, elige la opción más rigurosa y se las presenta justificadas para su visto bueno. Las que dependan del worldbuilding o de sus preferencias se le consultan aparte.

### 6.1 Correcciones a la v1 de este documento

1. **§1.2.** Decía que la presión del centro de la capa "difiere mucho" en la capa de arriba entre la media aritmética y SB81. **Es falso.** Con el tope a p = 0, las dos dan p_{3/2}/2 en la capa de arriba. Diferencia relativa calculada con los 20 niveles de M3N: 0 % en la capa 1, −0,68 % en la capa 2 (la máxima), bajando hasta −0,01 % en la capa 20.
2. **§1.4.** Suponía que `factor_estabilidad_louis` "ya da la función para el momento" (⚠️ comprobar). **No la da:** es solo la del calor. Louis, Tiedtke y Geleyn (1982), Tabla 1, verificada en el código de LMDZ (`cdrag_mod.F90`, que la cita) ✅:
   - momento, inestable: f_m = 1 − 2b·Ri / (1 + 3bc·C_N·√(z/z₀·|Ri|));
   - momento, estable: f_m = 1 / (1 + 2b·Ri / √(1 + d·Ri));
   - calor, inestable: f_h = 1 − 3b·Ri / (1 + 3bc·C_N·√(z/z₀·|Ri|));
   - calor, estable: f_h = 1 / (1 + 3b·Ri·√(1 + d·Ri));
   - b = c = d = 5.
   La función del calor de M3N coincide con esta tabla ✅. LMDZ usa √(|Ri|(1 + z/z₀)) en lugar de √(z/z₀·|Ri|): con z = 10 m y z₀ ≤ 0,01 m, la diferencia es ≤ 0,05 % ✅ (cálculo propio).
3. **§1.3.** Suponía que una "ráfaga mínima" es práctica habitual. El esquema publicado más cercano a M3N (radiación gris), Frierson, Held y Zurita-Gotor (2006), usa **ráfaga cero**: "zero gustiness velocity in this model, so surface fluxes are allowed to approach zero" ✅ (tesis de Frierson).
4. **§1.6 y `DISENO_FASE6_2.md`.** "Techo a 0,7 hPa" → la capa de arriba va de 0 a 7,35 hPa y su centro está a 3,7 hPa.
5. **Hallazgo nuevo:** el viento fijo de 5 m/s (`VIENTO`, `fase2b_atmosfera.py`) también entra en el **albedo del océano** (tabla de Cox-Munk) y en la **rugosidad del océano** (Z0M = 2·10⁻⁴ m, "Charnock a 5 m/s"). Con viento real, las tres cosas tienen que usar el viento de cada celda.

### 6.2 Las decisiones (aprobadas)

**1.1 Reparto del paso y acoplamiento** ✅ fuente / 🔶 diseño
- CAM euleriano (leapfrog, como M3N) usa un acoplamiento "process split": la física se calcula sobre el nivel n−1 y se aplica con 2Δt ✅ (descripción de CAM, cap. 3). El paso de la física de M3N (992 s) es exactamente 2 × 496 s.
- **Riesgo detectado:** aplicar la física de golpe cada 2 pasos tocaría solo una de las dos cadenas del leapfrog (pasos pares e impares) y excitaría el modo computacional.
- **Decisión:**
  - la física se calcula una vez cada 992 s;
  - su tendencia de T, u y v entra como forzamiento constante en los 2 pasos de la dinámica (cada cadena la recibe una vez);
  - el vapor, el suelo, el hielo y la superficie se actualizan directamente, con lo que el agua sigue cerrando de forma exacta.
- Aplicar la tendencia repartida y no de golpe reduce el ruido de las ondas de gravedad en CAM-SE ✅ (Gross et al. 2018, revisión del acoplamiento física-dinámica, fig. 5).
- **Prueba obligatoria:** física evaluada en n−1 (como CAM, por defecto) frente a en n; medir la amplitud del modo computacional y el cierre de energía.

**1.2** Unificar la presión del centro de la capa con SB81. Se mide el cambio en la v3.1, que se espera pequeño por el §6.1.1.

**1.3 Viento en superficie**
- El de la capa más baja del núcleo, en el centro de la celda. **Sin ráfaga** (Frierson et al. 2006 ✅), con solo un mínimo numérico para no dividir por cero en Ri.
- Si la validación en modo Tierra muestra un desacoplamiento irreal de las noches sobre tierra (oscilación diaria frente a las observaciones), se añade una ráfaga **con fuente verificada**: medido, no supuesto.
- El albedo y la rugosidad del océano pasan a usar el viento local. La rugosidad, con Charnock: constante ⚠️ por verificar en la fuente antes de programarla.

**1.4 Rozamiento y capa límite**
- **Arrastre en superficie:** τ = ρ·C_N·f_m(Ri)·|v|·v, con la f_m de Louis, Tiedtke y Geleyn (1982) ✅ (§6.1.2) y la misma rugosidad que el calor.
- **Mezcla vertical en la capa límite:** el perfil de K de Frierson et al. (2006) ✅. Capa superficial = 10 % de la capa límite (f = 0,1); tope donde el Ri global llega a 1; el mismo coeficiente para el momento, el calor y el vapor. Implícita en cada columna, sobre u, v, T y q.
- Combinar el arrastre de Louis con la mezcla de Frierson es diseño propio 🔶. Prueba de columna: con estratificación neutra tiene que salir el perfil logarítmico.

**1.5** Confirmado lo hecho (§5): el trazador lleva su propia masa de aire y se concilia conservando exactamente el agua.

**1.6 Capa esponja: no se pone de entrada.** La decide el paso 2 del §3 (§6.4). Si hace falta, tendrá que conservar el momento angular y devolver como calor la energía cinética que quite. Su forma se verificará entonces en la fuente ⚠️ (Shepherd et al., "Sponge layer feedbacks in middle-atmosphere models": la página del artículo no fue accesible el 06/10).

**1.7** Corrección de la difusión de T para que actúe como sobre superficies de presión, como CAM ✅ (`difcor.F90`: tcor ∝ ∇⁴p_s · B_k · ∂T/∂p; en σ pura, B_k = σ).

**1.8 Energía**
- La energía cinética que quitan la difusión del viento y el rozamiento en superficie se devuelve como calor, como hace CAM en `difcor.F90` ✅ ("to conserve total energy").
- El balance incluirá la energía cinética y Φ_s·p_s, y declarará por separado las pérdidas del filtro RAW y del filtro polar (ya no serán exactas a 10⁻¹²).
- Antes de decidir si hace falta un corrector global de energía, **se mide** el error. El umbral para decidirlo, con fuente ⚠️.

### 6.3 Equilibrio, climatología y coste (§1.9 nuevo)

**El problema (hallazgo del 06/10).** El criterio actual de `simular_fase2b` es que ninguna celda cambie su media anual más de 0,015 K, y el clima es el último año. Con el núcleo dinámico hay tiempo meteorológico, el criterio no se cumpliría nunca y un año pasa a ser una muestra, no el clima.

**Decisiones (Carlos, 06/10):**
- **Equilibrio:** medias globales anuales en una ventana de los últimos 5 años, con cuatro condiciones a la vez:
  - balance en el tope de la atmósfera |N| < 0,2 W/m²;
  - tendencia del aire a 2 m < 0,02 K/año;
  - tendencia del área de hielo < 0,1 % del océano por año;
  - tendencia del agua del suelo < 1 % de su capacidad por año.

  Son umbrales de partida 🔶: se fijan en firme midiendo el ruido de la primera simulación en modo Tierra. El modelo de 2 capas (v2.4.3) conserva su criterio actual.
- **Arranque caliente:** el océano, el hielo y el suelo salen del equilibrio del modelo de 2 capas con el mismo mapa y los mismos parámetros. El aire arranca en reposo. Prueba: el mismo equilibrio que desde un arranque frío.
- **Aceleración:** el salto geométrico se aplica solo a lo lento (entalpía del océano, hielo, suelo profundo, agua del suelo), **nunca** al aire ni al viento.
- **Climatología:** se promedian los años necesarios para que el error de la media sea < 0,1 K en la temperatura de cada banda de latitud y < 5 % en su precipitación (σ/√N_eff, con la autocorrelación entre años).
- **Exportación** (claves nuevas de `m3n-clima` v1):
  - medias de N años por día del año;
  - extremos absolutos de los N años;
  - desviación típica entre años;
  - datos horarios como media de N años por día y hora;
  - metadatos: años promediados, error estimado y nivel del modelo.
- **Comparaciones:** como el modelo es caótico, una diferencia entre dos simulaciones solo cuenta si supera unas 2 veces el error de la media. La prueba de "interruptor apagado = bit a bit" sigue igual.
- **Coste:** objetivo ~2 h por simulación de P3N, techo una noche. Todas las optimizaciones posibles, cada una con su prueba de igualdad bit a bit o de diferencia medida, acotada y documentada.
  - Medido el 06/10 en el entorno de la IA: núcleo 43 ms/paso (en el PC ~7), física de la v3.1 con agua 36 ms/paso, vapor 3D 7,5 ms/paso.
  - Estimación en el PC: ~10 min por año de P3N y ~15 min por año terrestre (±50 %).
- **Resolución de 2,5°:** se pensará cuando el planeta esté definido, quizá para una sola simulación muy larga.
- **Jerarquía de modelos** (rápido = v2.4.3 de 2 capas; completo = 6.3 en adelante): aceptada como idea; los Barridos se verán más adelante. `DISENO_H3N.md` no cambia por ahora.

### 6.4 Paso 2 del §3: núcleo con los niveles de M3N (v3.1-pre6)

- `python held_suarez.py --niveles m3n` (Held y Suarez, 1200 días, τ = 0,5 días, dt = 450 s, con los 20 niveles de `sigma_seminiveles`). En el PC de Carlos, ~25 min.
- El registro muestra cada 20 días el máximo de |u| en las capas con el centro por encima de 20 hPa. `analizar_held_suarez.py` compara esas capas entre la 1.ª y la 2.ª mitad del promedio.
- Prueba corta en el entorno de la IA (06/10): 21 días estables; el día 20, chorro de 39,0 m/s en 27,5° S, σ = 0,20; capas altas, |u| máx. 28,3 m/s.
- **Criterios para no necesitar la capa esponja** (🔶, fijados ANTES de ver el resultado):
  1. estable los 1200 días, sin valores no finitos;
  2. capas altas sin deriva: el máximo de |u zonal medio| cambia < 10 % entre la 1.ª y la 2.ª mitad;
  3. troposfera coherente con la corrida de niveles iguales (37,3 m/s en 28,9°): diferencias < 2 m/s y < 2° de latitud. Los niveles son distintos, así que no se espera igualdad exacta.

  Si falla 1 o 2, se diseña la capa esponja (1.6). Si falla solo 3, se analiza antes de decidir.

**Resultado (PC de Carlos, 07/10/2026; 25 min; media de los días 200 a 1200):**

| | Niveles iguales (6.2) | **Niveles de M3N** |
|---|---|---|
| Chorro N | 37,3 m/s en 28,9°, σ = 0,22 | **36,6 m/s en 29,3°, σ = 0,26** (mitades: 37,0 en 29,1° / 36,2 en 29,4°) |
| Chorro S | 37,2 m/s en 29,0° | **36,8 m/s en 29,2°** (mitades: 36,8 en 29,2° / 36,7 en 29,3°) |
| Viento del oeste en superficie | 3,3–3,5 m/s en 33,5° | **4,0 m/s en 34,0–34,1°** |
| Calor por remolinos (máx., σ = 0,85) | 33,1° | 33,9–34,0° |
| Momento por remolinos hacia el polo (máx., σ = 0,25) | 34–35 m²/s² en 24° | 37,4–37,9 m²/s² en 24,3° |
| Capas altas (centro < 20 hPa), \|u zonal medio\| máx. | — | **9,9 m/s** (mitades: 10,3 / 9,7; −6 %) |

- Registro cada 20 días: T entre 189 y 312 K y p_s entre 959 y 1030 hPa durante los 1200 días; el máximo instantáneo de |u| en las capas altas, entre 13 y 22 m/s desde el día ~100, sin tendencia.
- **Criterios:** 1 ✅ (estable); 2 ✅ (−6 %, por debajo del 10 %, y bajando); 3 ✅ (diferencias de 0,4–0,7 m/s y de 0,2–0,4°).
- **Decisión: sin capa esponja.** El viento algo mayor en superficie es coherente con las capas más finas abajo de M3N.
- **Límite de la prueba:** Held y Suarez fuerza una estratosfera isoterma y sin estaciones; no puede mostrar un chorro estratosférico de invierno, que aparecerá con la radiación de M3N. En la validación en modo Tierra (§2) se vigilará el viento de las capas altas con el mismo criterio 2; si se desboca, se diseña la esponja (1.6) con la fuente del §6.5.

### Referencias de este apartado

- CAM, descripción del núcleo euleriano (www2.cesm.ucar.edu/models/atm-cam/docs/description/node9.html) ✅.
- CAM, `difcor.F90` (código de CAM 5.4) ✅.
- Gross, M., et al. (2018), revisión del acoplamiento física-dinámica (arXiv:1605.06480) ✅.
- Frierson, D. M. W., Held, I. M., y Zurita-Gotor, P. (2006), J. Atmos. Sci. 63, 2548–2566 ✅ (detalles verificados en la tesis de Frierson).
- Louis, J.-F., Tiedtke, M., y Geleyn, J.-F. (1982), ECMWF Workshop on Planetary Boundary Layer Parameterization, Tabla 1 ✅ (vía LMDZ, `cdrag_mod.F90`).
- Mendonça, J. M., et al. (2016), THOR (arXiv:1607.05535) ✅: Held y Suarez, chorros hacia ~45° y ~250 hPa.

### 6.5 Paso 4: capa límite en columna (`fase6_capa_limite.py`, 06–07/10/2026)

Módulo autónomo, **sin conectar** (eso es el paso 5, con I16). Pruebas: `test_fase6_capa_limite.py` (7).

**Fuentes verificadas para este módulo** (además de las del §6.2):
- Perfil de K: `diffusivity.F90` de GFDL/Isca (opción `do_simple`, la de Frierson et al. 2006) ✅:
  - z < f h: difusividad de Monin-Obukhov;
  - f h ≤ z < h: K(f h)·(z/f h)·[1 − (z − f h)/(h − f h)]²;
  - z ≥ h: 0;
  - f = 0,1 y Ri_c = 1. h, con el Ri global desde la capa más baja e interpolación lineal.
- φ de Monin-Obukhov: `monin_obukhov.F90` de GFDL/Isca ✅:
  - inestable: φ_m = (1 − 16ζ)^(−1/4), φ_h = (1 − 16ζ)^(−1/2);
  - estable: 1 + ζ(5 + b ζ)/(1 + ζ), con b = 1/rich_crit = 0,5.
- Charnock (para el paso 5): z₀ = 0,11 ν/u* + α u*²/g, con α típico 0,018 ✅ (ECMWF, "Sea surface roughness and drag coefficient as function of neutral wind speed", 2010). α varía con el estado del mar (de 0,01 a 0,04); 0,018 es el valor típico que cita el ECMWF.
- Capa esponja (si hiciera falta): Shepherd, Semeniuk y Koshyk (1996, JGR, doi:10.1029/96JD01994) ✅ (existencia y resumen). Una esponja de relajación se acopla de forma artificial con la dinámica de abajo; la forma concreta, solo si el §6.4 la exige.

**Tres hallazgos de este paso:**
1. **La v3.0 calcula el intercambio con la superficie con coeficientes a 10 m (`Z_REF`), pero con la temperatura de la capa más baja,** cuyo centro está a ~184 m. Con la rugosidad de M3N, C_N(10 m)/C_N(184 m) = 2,02 sobre tierra y 1,61 sobre el mar. Los coeficientes se sobreestiman en ese factor.
   - El módulo usa la altura real de la capa más baja (z_a), como Frierson et al. y los GCM.
   - Al conectarlo (paso 5) cambiarán los flujos de la superficie, y se medirá cuánto.
   - El modelo de 2 capas no se toca.
2. **Difundir s = cp T + g z da un flujo de calor falso en una columna neutra:** con la hidrostática discreta, s no es constante en una columna isentrópica (hasta 77 J/kg entre las 5 capas bajas de M3N). Con K ~ 100 m²/s serían varios W/m².
   - **Solución:** el flujo de calor se calcula con el gradiente de θ (H = −ρ c_p K Π ∂θ/∂z, la forma de la teoría K) y se aplica en forma conservativa de energía (masa·Π como "masa" y D·Π_interfaz como conductancia).
   - Una columna neutra queda exactamente sin flujo, y Σ m c_p T se conserva exactamente.
3. **Calor por rozamiento: devolverlo capa a capa como −ΔKE (como `difcor.F90`) puede enfriar** una capa que recibe momento de otra, porque eso es transporte, no disipación.
   - Demostración propia, para Euler implícito y sumando por partes:
     Σ m (KE₁ − KE₀) = −dt Σ_i D_i |v₁,ᵢ − v₁,ᵢ₊₁|² − dt·c_s·|v₁,ₐ|² − ½ Σ m |v₁ − v₀|².
   - Cada término es ≥ 0 y se deposita donde ocurre:
     - lo de cada interfaz, a medias entre sus dos capas;
     - lo del suelo, en la capa más baja;
     - lo del paso implícito, en su capa.
   - El calentamiento es positivo en cada capa y la energía total cierra a redondeo.

**Pruebas** (`test_fase6_capa_limite.py`, todas superadas):
- Louis para el calor = `factor_estabilidad_louis` de M3N, bit a bit.
- f(0) = 1; las dos funciones de Louis son monótonas; en estable, el momento se frena menos que el calor.
- φ en valores exactos.
- Energía Σ m (c_p T + KE) conservada a < 10⁻¹⁴ (relativo), en columnas inestables, estables y mixtas.
- Momento: solo cambia por la tensión en superficie (10⁻¹⁰).
- Vapor: conservado (10⁻¹³) y ≥ 0.
- Calor por rozamiento ≥ 0 en cada capa, y su suma igual al diagnóstico.
- Columna en reposo e isentrópica: no cambia nada (< 10⁻¹⁰ K).
- **Perfil logarítmico:** columna neutra fina (80 niveles de 2 a 2000 m), estado estacionario con una fuerza uniforme. El salto de viento entre 5 y 60 m coincide con la solución exacta con un error de **0,04 %** (tolerancia del 1 %).

**Comprobación de cordura** (columna con los 20 niveles de M3N, atmósfera estándar, viento de 10 m/s, z₀ = 0,01 m):
- **De día** (suelo +3 K, 200 W/m²): h = 971 m, u* = 0,45 m/s, Ri = −0,19, K_m máx. = 79 m²/s.
- **De noche** (suelo −4 K, −30 W/m²): u* = 0,25 m/s, Ri = +0,25, C_m 3,3 veces menor, K_m máx. = 3 m²/s.
- z_a = 185 m.

**Nota sobre la rejilla de M3N.** La capa más baja tiene el centro a ~185 m. Con h ~ 1 km, f h ~ 100 m queda **por debajo** de ella: en M3N, la ley logarítmica entre z₀ y z_a la lleva la fórmula del arrastre, y las interfaces caen en la parte cuadrática del perfil de K. Por eso la prueba del perfil logarítmico se hace en una columna fina (el módulo vale para cualquier rejilla).

### 6.6 Coste del núcleo: dónde se va el tiempo (06–07/10, entorno de la IA, sin cambiar nada)

Perfil de 600 pasos con los niveles de M3N:

| Parte | Porcentaje del tiempo |
|---|---|
| `tendencias` (lado derecho explícito) | ~37 % |
| FFT (40 transformadas pequeñas por paso: filtro polar, semiimplícito, hiperdifusión) | ~13 % |
| `aplicar_hiperdifusion` | ~13 % |
| `np.roll` (21 llamadas por paso) | ~6 % |
| RAW y resto | el resto |

- Con campos de ~52 000 valores domina el coste de llamar a cada operación de numpy, no la aritmética.
- **Plan de optimización** (cada paso con su prueba de igualdad bit a bit o de diferencia medida y acotada):
  1. fusionar `tendencias` en núcleos de numba;
  2. agrupar las FFT de varios campos en una sola llamada;
  3. sustituir `np.roll` por índices.
- **Estimación ⚠️ (por medir):** ~2 veces más rápido el núcleo.

**Resultado de la optimización (v3.1-pre8, 07/10/2026):**
- **`tendencias` compilada con numba** (`fase6_nucleo_nb.py`):
  - cada expresión reproduce el orden de operaciones de numpy → resultado **idéntico bit a bit**;
  - en el entorno de la IA, de 14,3 a 2,4 ms por llamada (×5,9).
- **Hallazgo:** el `np.log` vectorizado de numpy y el `log` que usa numba pueden diferir en el último bit. Aparecía en cuanto p_s dejaba de ser uniforme (1 ulp en dT, du, dv desde el paso 1). Por eso `log(p_s)` se calcula con numpy y se pasa al kernel.
- **Filtros polares en lote** (`_filtro_varios`): las FFT de cada fila son independientes; de 16 a 8 llamadas por paso, idéntico bit a bit.
- **Comprobaciones:**
  - Held y Suarez con los niveles de M3N: 300 pasos con la v3.1-pre8 = 300 pasos con la v3.1-pre7, **bit a bit** (p_s, T, u, v);
  - pruebas nuevas en `test_fase6_2.py`: tendencias compiladas = numpy (con montañas), filtro en lote = uno a uno, 20 pasos de Held y Suarez con numba = sin numba, todo bit a bit.
- **Paso completo** (entorno de la IA, Held y Suarez con los niveles de M3N): de 33 a 24 ms. **×1,4, no ×2 como se estimó en el §6.6.** Era una estimación sin medir y queda corregida.
  - Lo que queda son bloques de < 1 ms: las partes lineales del semiimplícito, los productos de matrices (BLAS), la resolución por ondas, la hiperdifusión (~3,7 ms, sobre todo productos complejos con BLAS) y la aritmética de RAW.
  - Sustituir las operaciones de BLAS rompería la igualdad bit a bit (su orden interno de suma no es reproducible), y compilar el resto ahorraría ~1–2 ms con bastante código nuevo. **Decisión: parar aquí.**
- **Consecuencia para el coste** (§6.3, estimación ⚠️ por medir con el modelo acoplado): con la física de la v3.1, ~36 ms por paso de física en el entorno de la IA, hoy la parte más cara, el paso acoplado queda en ~20 ms en el PC de Carlos frente a ~25 → ~8 min por año de P3N; con el arranque caliente (~16 años), ~2,1 h. Si hace falta bajar más, lo siguiente sería la física de columna, no el núcleo.
- `USAR_NUMBA = False` (en `fase6_nucleo.py`) vuelve a la versión de numpy; sin numba instalado se usa sola.

### 6.7 Paso 5: acoplamiento física–núcleo (interruptor I16 `nucleo_dinamico`) — DISEÑO (07/10/2026)

Basado en las decisiones aprobadas del §6.2 y del §6.3. Lo que la IA decide aquí por ser técnico va marcado 🔶 y se puede revisar.

#### 6.7.1 Arquitectura
- I16 será un **paso nuevo dentro de `simular_fase2b`** (`paso_acoplado`), como `paso_agua` para la v3.1.
  - Así se reutilizan la superficie, el hielo, el registro horario, la caché, los puntos de control y los diagnósticos.
  - Necesita I10 e I11 (y admite I12–I14).
  - **Con I16 apagado, todo es idéntico bit a bit** (prueba T1).
- **Estado nuevo:** el par leapfrog del núcleo (`ant`, `act`: p_s, T, u, v). T_atm pasa a ser la T del núcleo. La forzante de la física del paso anterior va en el estado y en el punto de control.

#### 6.7.2 Orden de un paso de la física (Δt = 992 s en P3N, 900 s en modo Tierra)
1. **Dinámica:** 2 subpasos de Δt/2 (leapfrog semiimplícito + RAW + hiperdifusión). La tendencia de la física del paso anterior entra como forzamiento constante en los dos subpasos, así que cada cadena del leapfrog la recibe una vez (§6.2, decisión 1.1). Se acumulan los flujos de masa (F_u, F_v, W)·Δt de los dos subpasos.
2. **Vapor:** transporte 3D positivo y conservativo con esos flujos (`fase6_trazadores.py`) y conciliación q ← q·m_tr/Δp_núcleo (§5).
3. **Columna:** `COL.actualizar_ps(p_s del núcleo)` con `presion_capa="sb81"` (decisión 1.2). También se recalcula lo que en la v3.1 se fijaba una sola vez y depende de p_s:
   - la masa de la capa baja;
   - el factor de masa de la luz (`preparar_luz(factor_masa=COL.ps/P0)`) ⚠️: hoy se fija al empezar; habrá que pasarlo en cada paso;
   - las capacidades.
4. **Física de columna**, en el orden de `paso_agua`:
   - radiación, superficie e hielo;
   - calor sensible y evaporación, con el **viento real** de la capa más baja (en los centros) y la **altura real** z_a (§6.5, hallazgo 1);
   - convección y condensación; lluvia y nieve; suelo; hielo;
   - **capa límite** (`fase6_capa_limite.py`) sobre u, v, T y q, con el arrastre y el calor por rozamiento local y positivo.
   - **Se apagan** la difusión horizontal del aire (`fact_atm30`) y la del vapor (`fact_vapor`). **El océano conserva su difusión** hasta la Fase 6b.
5. **Tendencias de la física:**
   - F_T = (T después − T antes)/Δt;
   - F_u y F_v: las de la capa límite, de los centros a las caras (`tendencia_a_caras`);
   - F_ps = 0 🔶: la física no cambia la masa de aire seco. La masa del vapor que llueve o se evapora no se resta de p_s, como en muchos modelos de "masa seca"; hay que medir su orden de magnitud.

#### 6.7.3 Superficie con viento real (§6.2, decisiones 1.3 y 1.4)
- Coeficientes con z_a, la estabilidad de Louis para el calor y para el viento (`fase6_capa_limite.py`), y |v_a| del núcleo, sin ráfaga (solo el mínimo numérico).
- **Albedo del océano (Cox-Munk) con el viento local:** hoy es una tabla en μ para 5 m/s; pasa a ser una tabla en (μ, viento).
- **Rugosidad del océano por Charnock:** z₀ = 0,11 ν/u* + 0,018 u*²/g (§6.5 ✅), iterada con u*.
- El hielo marino y la nieve usan los mismos coeficientes.

#### 6.7.4 Energía (§6.2, decisión 1.8)
- E = Σ (c_p T + K + L_v q) Δp/g + Φ_s p_s/g + las entalpías de la superficie, del hielo y de la nieve.
- **Pérdidas declaradas en el diagnóstico, cada una medida por separado:**
  - filtro RAW;
  - filtro polar;
  - hiperdifusión: su energía cinética se devolverá como calor, con la corrección de la difusión de T a superficies de presión (decisión 1.7); hay que añadirlas al núcleo como opciones, apagadas en Held y Suarez;
  - residuo de la conversión caras-centros (medido el 07/10: ~0,2 % del calor por rozamiento, ~0,004 W/m², §6.7.6).
- **Criterio:** la suma de las pérdidas no explicadas < 0,05 W/m² 🔶. Si se supera, se estudia con fuente un corrector global (decisión 1.8).

#### 6.7.5 Agua
- Cierre exacto, como en la v3.1 (~10⁻¹³). El transporte y la conciliación conservan la masa de agua exactamente (§5).

#### 6.7.6 Piezas ya hechas (v3.1-pre9, 07/10)
- **`Columna.actualizar_ps`** y la opción `presion_capa="sb81"`.
  - Con la opción por defecto, los 18 atributos de la columna y los 28 resultados de una simulación corta de la v3.1 (agua y convección) son **idénticos bit a bit** a los de la v3.1-pre8.
  - Con "sb81", p_m = la del núcleo bit a bit.
- **Bug latente corregido:** la caché de la geometría de la convección húmeda (`fase31_agua._geometria_conveccion`) usaba como clave `id()` de los arrays y tres valores sueltos. Con p_s variable habría devuelto en silencio la geometría de otro paso. Ahora compara el contenido completo; lo comprueba una prueba que cambia una sola celda interior.
- **`fase6_acoplamiento.py`:** `viento_en_centros`, `tendencia_a_caras` (su traspuesta sin pesos) y `calor_rozamiento_exacto` (diagnóstico).
- **Medido tras 30 días de Held y Suarez** con los niveles de M3N (capa límite en los centros, Δt = 992 s):
  - el momento zonal quitado en los centros frente al de las caras difiere un **1,6·10⁻⁴** relativo;
  - el calor por rozamiento: 1,568 W/m² (positivo, en los centros) frente a 1,572 W/m² (la pérdida de K con la definición del núcleo);
  - aplicar la segunda celda a celda enfriaba el 1,2 % de las (celda, capa), hasta 5·10⁻³ K por paso, y escalarla por columnas sale muy mal condicionado (factores de 0,46 a > 10⁵).
  - **Decisión 🔶:** calor local y positivo de la capa límite; el residuo global se declara.
- Pruebas: `test_fase6_acoplamiento.py` (6).

#### 6.7.7 Pruebas del acoplamiento (antes de validar en modo Tierra)
- **T1:** I16 apagado = v3.1-pre8 bit a bit (todas las pruebas existentes).
- **T2:** energía con I16 encendido (2 días): ΔE = TOA neto·Δt + las pérdidas declaradas, con el resto < 0,05 W/m².
- **T3:** agua exacta (10⁻¹²).
- **T4:** atmósfera isoterma en reposo sobre montañas, sin radiación: sigue en reposo.
- **T5:** modo computacional, física evaluada en n−1 frente a n (decisión 1.1): amplitud de la oscilación par-impar.
- **T6:** punto de control: cortar y reanudar da lo mismo bit a bit.

#### 6.7.8 Entregas previstas
| Versión | Contenido |
|---|---|
| pre9 (hecha) | Columna con p_s variable, caché corregida, conversión caras-centros, este diseño |
| pre10 | Opciones del núcleo (energía cinética de la hiperdifusión → calor; corrección de T a superficies de presión); luz con p_s variable; superficie con z_a y viento real (Cox-Munk 2D, Charnock) como funciones probadas, sin conectar |
| pre11 | `paso_acoplado` (I16) con T1–T6 |
| pre12 | Equilibrio y climatología (§6.3), arranque caliente, prueba corta en modo Tierra |
| — | Paso 6: validación en modo Tierra en el PC de Carlos, con los criterios del §2 escritos antes |

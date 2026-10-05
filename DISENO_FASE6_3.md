# Fase 6.3 de M3N: acoplar el núcleo dinámico con la física de M3N (BORRADOR v1)

**Estado: BORRADOR v1 (05/10/2026), sin código.** Requisito: que la 6.2 esté cerrada con Held y Suarez.

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

- **Por qué hace falta:** los 20 niveles de M3N llegan a 0,7 hPa. La 6.2 (intento 1 de Held y Suarez) demostró que sin disipación en el tope los vientos de las capas altas crecen sin límite.
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

- 1.2: unificar la presión del centro de la capa en la definición de SB81.
- 1.3: viento en superficie a partir del núcleo, con una ráfaga mínima.
- 1.4: arrastre de fórmula aerodinámica y difusión vertical en la capa límite.
- 1.6: capa esponja (y si devuelve la energía como calor).
- El orden de trabajo de §3.

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

# BORRADOR — Fase 6 de M3N: circulación atmosférica con un núcleo dinámico PROPIO

> **Nomenclatura (08/10/2026, M3N 3.12.1):** este documento se llamaba `DISENO_FASE6_BORRADOR.md`. Fases y versiones van con la numeración normalizada; la correspondencia con los nombres anteriores está en `NOMENCLATURA.md` y `CHANGELOG.md`.

**Estado: v1 APROBADO (05/10/2026). En desarrollo: 6.1 ✅ (`DISENO_FASE_6.1.md`), 6.2 casi ✅ (`DISENO_FASE_6.2.md`), 6.3 empezada (`DISENO_FASE_6.3.md`).**

**Principio (Carlos, 05/10/2026):** M3N es un modelo propio, construido pieza a pieza y de forma exhaustiva y rigurosa. **No se usan otros modelos**: ni sus núcleos ni su física como sustitutos.

Los artículos y los modelos publicados se usan solo para dos cosas:
- como **fuentes** de las ecuaciones y los métodos numéricos;
- como **pruebas de validación estándar**, con resultados publicados que nuestro código tiene que reproducir.

(La versión v0 de este borrador proponía usar ExoPlaSim o núcleos ajenos; queda descartada.)

Etiquetas: ✅ verificado en la fuente; ⚠️ pendiente de verificar.

---

## 0. Por qué la circulación es ahora lo prioritario

Las validaciones de la v3.0 y la v3.1 (`DISENO_FASE_2.3.md` §13 bis; `DISENO_FASE_5.1.md` §8 bis y §10) dan tres fallos con la misma causa: **todo el transporte horizontal es una difusión con un solo coeficiente**.

1. **Polos demasiado cálidos y sin hielo marino.**
2. **Troposfera libre demasiado húmeda**, con una humedad relativa de ~0,85 en todas partes. En la realidad la seca la **subsidencia**: el aire que baja en los subtrópicos se saturó por última vez muy arriba y muy frío (Pierrehumbert, Brogniez y Roca 2007 ✅). La difusión del vapor no puede producir esa sequedad. Por eso el vapor radiativo (I15) se dispara.
3. **Tropopausa plana**, a ~10 km en todas las latitudes.

## 1. Qué circulación cabe esperar en P3N (cálculo propio)

| Magnitud | P3N / Tierra |
|---|---|
| Número de Rossby térmico gH/(Ω²a²) (mismo Δθ y H) | ~1,17 |
| Anchura de la célula de Hadley (∝ Ro^½, Held y Hou 1980 ✅ que existe; fórmula ⚠️) | ~1,08 |
| β = 2Ω/a | ~1,64 |
| Radio de deformación / radio del planeta | ~1,13 |

- Tres células por hemisferio, como en la Tierra.
- Con una inclinación del eje de 1,7°, casi no hay estaciones: la zona de lluvias ecuatorial apenas se moverá.

## 2. Qué es un núcleo dinámico y qué hay que construir

Hasta ahora M3N calcula, en cada columna de aire, la radiación, la convección, el agua y el calor que intercambia con el suelo. Lo que falta es **el movimiento del aire**: las ecuaciones primitivas, que son
- la conservación de la masa;
- el movimiento horizontal (fuerza de presión, Coriolis, advección);
- la termodinámica;
- el transporte del vapor.

Todas se resuelven en la rejilla de M3N, con las 20 capas en coordenada σ que ya existen (Phillips 1957 ✅).

El "núcleo dinámico" resuelve ese movimiento. La física actual de M3N se le suma en cada paso, y la difusión desaparece de la atmósfera (el océano la mantiene hasta la Fase 6.5).

## 3. Decisiones de diseño (opciones, con fuentes)

### 3.1 Discretización horizontal

**Opción G (puntos de rejilla sobre la rejilla actual de 5°, recomendada):**
- Ventaja: **no hay que cambiar de rejilla**. Todo M3N (mapas, C3N, exportación) usa ya la de 72 × 36.
- Precedente: el GISS Model II, con 8° × 10° de resolución y 9 capas (Hansen et al. 1983 ✅), dice textualmente que "los rasgos principales del clima global pueden simularse de forma realista con una resolución tan gruesa como 1000 km" ✅. Usaba una rejilla B de Arakawa ✅.
- **Rejilla escalonada de Arakawa:**
  - **C:** u en las caras este y oeste de cada celda, v en las caras norte y sur; Arakawa y Lamb 1977 ✅.
  - **B:** u y v juntos en las esquinas.
  - Cuál conviene depende del cociente entre el radio de deformación y el tamaño de la celda ⚠️ (Randall 1994, por verificar). En P3N la celda mide unos 408 km en el ecuador.
- **Problema de los polos:** las celdas se estrechan hacia los polos (a 87,5° miden ~18 km de ancho en P3N), y el paso de tiempo estable se desplomaría.
  - El remedio clásico es un **filtro polar de Fourier**, que amortigua las ondas zonales cortas en las latitudes altas (planetWRF ✅).
  - El GISS usaba un filtro de Shapiro ✅.

**Opción S (espectral):**
- Es el método de los grandes modelos (Bourke 1974 ✅; Hoskins y Simmons 1975 ✅), sin problema de los polos.
- Pero necesita **latitudes gaussianas** (✅ ECMWF) y no la rejilla actual: habría que interpolar entre dos rejillas en cada paso, y C3N y la exportación seguirían en la regular.
- Es más elegante, pero añade una capa de complejidad y de error.

### 3.2 Discretización vertical

- σ, ya en uso (v3.0).
- Esquema de **Simmons y Burridge (1981)** ✅, que conserva la energía y el momento angular (lo dice su título; ecuaciones ⚠️, hay que leerlas). Es el que usa el núcleo espectral del GFDL ✅.

### 3.3 Avance en el tiempo

- **Leapfrog semiimplícito:** las ondas de gravedad (~300 m/s ✅ Lynch) se tratan implícitamente. Así el paso lo limita el viento y no esas ondas (✅ GFDL; Robert, Henderson y Turnbull 1972 ✅).
- **Filtro de tiempo:** Robert–Asselin (Asselin 1972 ✅), o su versión mejorada **RAW** (Williams 2009 ✅: α = 0,53, ν ≈ 0,2), que conserva mejor la media.
- **Objetivo:** que el paso de la dinámica sea el de la física, 992 s, o un divisor suyo.
  - Cálculo propio: en el ecuador, con 408 km y vientos de hasta ~100 m/s, el límite de la advección ronda los 4000 s.
  - Cerca de los polos lo decide el filtro polar.

### 3.4 Transporte del vapor y de otros trazadores

- Esquema conservativo y que nunca dé valores negativos: de volúmenes finitos con limitador, o semilagrangiano ⚠️ (por estudiar).
- La masa de agua tiene que seguir cerrando como ahora (10⁻¹³).

## 4. Plan por etapas, cada una validada antes de pasar a la siguiente

| Etapa | Qué se construye | Prueba que tiene que superar (publicada) |
|---|---|---|
| **6.1** | Aguas someras en la esfera: una sola capa, con su rejilla escalonada, Coriolis y filtro polar | Williamson et al. (1992) ✅: caso 1 (advección de una "campana" por encima del polo, solución exacta), caso 2 (flujo zonal estacionario, solución exacta: u0 = 2πa/12 días, g·h0 = 2,94×10⁴ m²/s², también con el eje girado), caso 5 (flujo sobre una montaña de 2000 m), caso 6 (onda de Rossby-Haurwitz). Conservación de la masa exacta y de la energía casi exacta. |
| **6.2** | Ecuaciones primitivas en seco, con 20 capas σ, Simmons–Burridge y paso semiimplícito | Jablonowski y Williamson (2006) ✅: estado de equilibrio que debe mantenerse 30 días, y onda baroclínica (a 5° solo se puede comparar cualitativamente). **Held y Suarez (1994)** ✅: forzamiento estándar (calentamiento por relajación a una T de equilibrio y rozamiento cerca del suelo); clima esperado, corrientes en chorro de ~30 m/s hacia 45° ✅ (la altura exacta ⚠️). Unos 1000 días de simulación ✅. |
| **6.3** | Acoplamiento con la física de M3N (radiación, convección, agua, superficie, hielo); la difusión atmosférica se apaga | Modo Tierra: vientos medios por latitud, células de Hadley y Ferrel, humedad relativa de la troposfera libre, temperatura polar y hielo marino, tropopausa, transporte de calor (Trenberth y Caron 2001), precipitación (GPCP). Las referencias observadas se verificarán una a una. |
| **6.4** | P3N | Comparación con la v2.4.3 y la v3.1; análisis de dónde quedan los desiertos y los vientos. |
| **6.6** (añadida el 09/10/2026) | Arrastre de ondas de gravedad que conserva el momento (el momento que llegaría al tope se deposita en las capas de arriba), después de la Fase 5.2 | Por decidir en su diseño; motivo y fuentes en `DISENO_FASE_6.3.md` §6.17 (Shaw y Shepherd 2007). Sustituye a la capa esponja. |

Cada etapa, igual que hasta ahora: diseño, decisiones de Carlos, código con interruptores, pruebas automáticas y validación documentada.

## 5. Coste esperado (honesto)

- Es el paso más grande de M3N: **semanas de trabajo**.
- La etapa 6.1 es la más importante para aprender y asentar la base (rejilla, filtro polar, avance en el tiempo). Es corta y se puede validar del todo con soluciones exactas.
- CPU (estimación propia ⚠️): la dinámica a 72 × 36 × 20 debería costar del orden de la física actual por paso, a medir en la 6.1.
- Una simulación de equilibrio de Held y Suarez son ~1000 días: en el PC de Carlos, previsiblemente minutos ⚠️.

## 6. Decisiones de Carlos

1. ¿Puntos de rejilla sobre la rejilla actual (G, recomendada) o espectral (S)?
2. Si es G: ¿rejilla B o C? La IA propone estudiarlo en la 6.1 con el caso 5 de Williamson y la referencia de Randall (1994).
3. ¿Se empieza por la etapa 6.1 (aguas someras) como primer paso concreto?
4. Con 1,7° de inclinación casi no hay estaciones: ¿se revisa la inclinación por motivos de worldbuilding (monzones) antes de construir la circulación, o se deja?

## Referencias

- Williamson, D. L., et al. (1992). J. Comput. Phys. 102, 211–224 ✅ (informe técnico ORNL/TM-11895: osti.gov/servlets/purl/5232139). Valores del caso 5 (u0, h0) ⚠️.
- Held, I. M., y Suarez, M. J. (1994). BAMS 75, 1825–1830 ✅ (constantes vía la documentación de MITgcm y Thatcher y Jablonowski 2016).
- Jablonowski, C., y Williamson, D. L. (2006). QJRMS 132, 2943–2975 ✅.
- Arakawa, A., y Lamb, V. R. (1977). Methods in Computational Physics 17, 173–265 ✅ (cita).
- Hansen, J., et al. (1983). Mon. Wea. Rev. 111, 609–662 ✅ (GISS Model II, 8° × 10°).
- Bourke (1974); Hoskins y Simmons (1975); Robert, Henderson y Turnbull (1972) ✅ (existencia).
- Asselin (1972) ✅; Williams (2009), Mon. Wea. Rev. 137, 2538 ✅.
- Simmons, A. J., y Burridge, D. M. (1981). Mon. Wea. Rev. 109, 758–766 ✅ (cita; ecuaciones ⚠️).
- Phillips, N. A. (1957). J. Meteor. 14, 184–185 ✅.
- Pierrehumbert, Brogniez y Roca (2007) ✅.
- Libros de apoyo: Durran, *Numerical Methods for Fluid Dynamics*, 2ª ed., 2010 ✅; Vallis, *Atmospheric and Oceanic Fluid Dynamics*, 2ª ed., 2017 ✅; Washington y Parkinson, *An Introduction to Three-Dimensional Climate Modeling*, 2ª ed., 2005 ✅.


## 7. Decisiones de Carlos (05/10/2026)

- **Camino B: núcleo dinámico propio.** Esta decisión **sustituye explícitamente** el plan de `DISENO_FASE_5.md` §3.6 ("Fase 6 (decidido)": mantener la D como remolinos más una célula de Hadley parametrizada según Siler, Roe y Armour 2018).
  - La IA no señaló esa sustitución al proponer este borrador; se hizo después, a petición de Carlos.
  - Alternativas que se descartaron: A, la de Siler et al.; y un camino intermedio, A ahora y B después.
- **En la atmósfera la D desaparece.** Las borrascas se resuelven explícitamente. Quedan como parametrizaciones de lo que está por debajo de la escala de la rejilla: la hiperdifusión (que solo disipa en la escala de la rejilla), la capa límite y la convección.
  - Si la validación en modo Tierra muestra que falta transporte hacia los polos, porque los remolinos están mal resueltos a 5°, se estudiará una difusión residual **medida, no supuesta**.
- **Rejilla C, filtro polar desde 60°, RAW con semiimplícito:** confirmadas (`DISENO_FASE_6.1.md` §7).
- **Pendiente:** la inclinación del eje (worldbuilding).

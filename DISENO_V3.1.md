# DISEÑO — M3N v3.1: ciclo del agua (Fase 5a) sobre la atmósfera de N capas

**Estado (05/10/2026): código hecho y probado (conservación exacta). Calibración y validación en curso (§9).** Los interruptores I11–I14 están **apagados por defecto**: el clima exportado de P3N no cambia hasta que Carlos decida encenderlos.

Por decisión de Carlos (05/10/2026, noche), la IA desarrolla la v3.1 con los criterios de siempre y se los explica después. Las decisiones de fondo ya estaban tomadas en `DISENO_FASE5A.md` (§1 y §11). Este documento describe **cómo** se han implementado, qué se ha cambiado respecto al borrador de la 5a por tener ya N capas, y qué queda pendiente.

Etiquetas: (a) física universal, (b) propiedad de P3N, (c) calibración terrestre. Estado: ✅ verificado, ⚠️ pendiente.

---

## 0. Resumen

La v3.1 añade el agua al modelo de N capas de la v3.0:
- **vapor propio en cada capa**: decisión H2, generalizada a N capas;
- **evaporación** del océano y del suelo, y **sublimación** de la nieve y del hielo marino;
- **transporte** del vapor por remolinos, con la misma difusividad que el calor;
- **convección húmeda**: Betts-Miller simplificado de Frierson (2007), traducido del código de Isca;
- **condensación de gran escala**: lluvia o nieve;
- **cubo de agua del suelo** (Manabe 1969);
- **nieve en tierra**, que se funde, tiene tope y tiene albedo;
- **propiedades térmicas del suelo según su agua** (pedido de Carlos);
- **albedo de la nieve y del hielo según el espectro de S3N** (pedido de Carlos, también para el hielo).

Todo **conserva exactamente** la energia (aire c_p·T, vapor L·q, superficie, hielo, nieve) y el agua. Las pruebas lo comprueban: cierre de energía ~10⁻¹², cierres de agua ~10⁻¹³ kg/m².

## 1. Interruptores (todos necesitan I10)

| | Nombre | Contenido |
|---|---|---|
| I11 | `ciclo_agua` | Vapor, evaporación, transporte, condensación, lluvia y nieve, cubo, nieve en tierra |
| I12 | `conveccion_humeda` | Betts-Miller simplificado, y el ajuste seco pasa al adiabático **seco** g/c_p de cada planeta |
| I13 | `suelo_termico_agua` | Conductividad y capacidad del suelo según su agua |
| I14 | `albedo_espectral` | Albedo de la nieve y del hielo con el espectro de la estrella |
| I15 | `vapor_radiativo` | PROTOTIPO (§8 bis): el infrarrojo depende del vapor del modelo |

Con I11–I14 apagados, la v3.0 queda exactamente igual (prueba V0).

**Cambio de numeración:** en `DISENO_FASE5A.md` estaban previstos como I10–I13. I10 lo ocupó la atmósfera de N capas.

## 2. Saturación (a) ✅

- e_s sobre agua y sobre hielo con la fórmula de Ambaum (2020). Constantes de MetPy, como en `DISENO_FASE5A.md` §3.1.
- q_s = ε·e_s / (p − (1 − ε)·e_s).
- dq_s/dT exacta: d ln e_s/dT = L(T)/(R_v T²).
- En la convección, una **tabla** cada 0,005 K (como `escomp` de Isca), con un error < 3×10⁻⁷.
- **Calor latente constante en la energía** (decisión 1, opción A): L_v = 2,50084×10⁶ J/kg; L_f = 3,34×10⁵ J/kg.
  - L_f es el **mismo valor que ya usaba el hielo marino** (y que usa Isca), para que la nieve y el hielo cierren la energía entre sí. MetPy da 3,337×10⁵ (diferencia de 0,09 %).
- La saturación sobre agua se usa en el océano, el suelo y toda la atmósfera; sobre hielo, en el hielo marino y en la nieve.
  - **Limitación:** la condensación en la atmósfera usa siempre la saturación sobre agua, como en Isca.

## 3. Evaporación, sublimación y rocío (a + c) — `DISENO_FASE5A.md` §3.3

E = ρ·C_H·U·β·(q_s(T_sup) − q_aire), con C_E = C_H y la misma corrección de estabilidad de Louis que el calor sensible.
- **Implícita en el vapor de la capa baja** (DISENO_FASE5A §3.3, Riesgo 4) y explícita en la temperatura de la superficie.
  - Margen de estabilidad comprobado: con el suelo más fino (seco, ~4×10⁴ J/m²/K) y ~20 W/m²/K de derivada del calor latente, el factor explícito es < 0,5.
- **Tierra:** fracción sin nieve, con β del cubo y vapor sobre agua; fracción con nieve, con β = 1 y vapor sobre hielo (sublimación, L_s).
- **Hielo marino:** sublimación con L_s, dentro del Newton de su superficie. La masa sublimada sale del hielo.
- **Límites:** no se evapora más agua del cubo ni más nieve de las que hay.
- E < 0 (rocío o escarcha) va al cubo o a la nieve.
- **Viento: 5 m/s también en el modo Tierra.** La decisión 7 de Carlos pide el viento medio real de la Tierra, pero **no se ha encontrado una fuente verificada** (los valores rondan 6–7,5 m/s según el producto). Por la regla de no usar en el código valores ⚠️, se mantiene 5 m/s.
  - Consecuencia: en el modo Tierra la evaporación sale algo por debajo de la real.
  - **PENDIENTE:** fuente para el viento oceánico medio.

## 4. Transporte del vapor — cambio respecto al borrador C′

El borrador (DISENO_FASE5A §3.6, opción C′) difundía la **energía estática húmeda cerca de la superficie**, h = c_p·T + g·z + L·q_aire, con un solo D.

Con N capas y vapor propio en cada una se probó primero exactamente eso: masa de vapor transportada = (D/c_p)·∇q_aire, repartida en la columna según el perfil de vapor. **Resultó inestable.** Cuando la capa baja se seca un paso (condensación o mezcla) y el vapor de la columna sigue arriba, la "capacidad" efectiva W/q_aire se dispara. Cada paso se mueven cantidades enormes de vapor entre vecinas, y se bombea vapor hacia las capas altas: humedades de 0,7 en la capa más alta en 40 pasos.

**Solución adoptada:**
- Cada capa del transporte (σ ≥ 0,25, las mismas que reciben el calor seco) **difunde su propia humedad con la misma difusividad de remolinos que el calor**, κ = D/(c_p·M_tr), siendo M_tr la masa de esas capas.
- Es implícita, monótona (nunca da vapor negativo) y conserva el agua capa a capa: en coordenada σ, la fracción de masa de cada capa es la misma en todas las celdas.
- La matriz es fija y se factoriza una sola vez.

**Consecuencias:**
- **Mismo κ para el calor y el vapor** (un solo parámetro libre, D, como pedía la decisión 2).
- El flujo total de vapor es (D/c_p)·∇q̄, con q̄ la humedad media en masa de la troposfera, en vez de la de superficie. Como el vapor se concentra abajo, ∇q̄ es menor que ∇q_aire, y **el transporte latente por unidad de D es menor que en la C′ pura**.
- El reparto latente/seco que salga es **validación independiente**. En la Tierra es ~50/50 en latitudes medias (Hwang y Frierson 2010, citando a Trenberth y Stepaniak 2003, ✅ secundaria).

## 5. Convección húmeda (I12) — Betts-Miller simplificado (Frierson 2007)

Traducción fiel de `qe_moist_convection.F90` de Isca (código fuente ✅):
- CAPE de una parcela desde la capa baja, con temperatura virtual;
- perfil de referencia: el adiabático húmedo de la parcela, con RH_ref = 0,8;
- relajación con τ = 7200 s;
- convección profunda: llueve, y la entalpía húmeda de la columna se conserva exactamente;
- convección somera: no llueve; conserva la entalpía y el agua por separado.

Diferencias con Isca:
1. El LCL se resuelve por Newton directo; Isca lo tabula con el mismo Newton.
2. Saturación de Ambaum tabulada, en lugar de la tabla de Isca.
3. **Tope de la convección en p = 0,1·p_superficie**: unos 16 km en la Tierra, alrededor de la tropopausa tropical.
   - Isca, con muchas capas, nunca llega arriba: la flotabilidad se acaba antes de la tropopausa.
   - Con 20 capas y la capa más alta en p → 0, la parcela podía seguir flotando hasta arriba y meter en esa capa minúscula temperaturas y humedades absurdas.
4. La capa más alta, con p_arriba = 0, usa el doble del tramo entre su centro y su base en el ln p de la CAPE. Isca obtendría infinito.

Otros detalles:
- **RH_ref (c) es incierto:** 0,6 en el control de Frierson (2007b) y 0,8 por defecto en Isca. **Prueba de sensibilidad obligatoria.**
- **Ajuste seco:** con I12, al adiabático **seco** g/c_p de cada planeta (9,05 K/km en P3N), en vez de los 6,5 K/km provisionales. Así el gradiente vertical sale de la física, como pedía la decisión 6.
  - El ajuste seco ahora **también mezcla el vapor** en los tramos que mezcla.
  - Sin esa mezcla, el vapor evaporado se quedaba en la capa baja, de ~430 kg/m².

## 6. Condensación de gran escala y precipitación

- Isca `lscale_cond` con `hc = 1` y sin reevaporación (✅ código; Manabe et al. 1965 ✅). Un paso de Newton con el calentamiento latente.
- **RH crítica = 1** (y no < 1 como preveía `DISENO_FASE5A.md` §3.7). La RH < 1 era para dos capas gruesas; con 20 capas, la resolución vertical es la de un modelo de capas finas (Isca, Frierson 2006).
  - **Sensibilidad pendiente:** 0,9.
  - Los modelos con nubes usan RH críticas de 0,7–0,9 para la **fracción de nubes** (CAM3 ✅, Met Office ✅), no para la condensación.
- **Lluvia o nieve:** rampa lineal con el aire a 2 m, toda nieve a ≤ 0 °C y toda lluvia a ≥ 2 °C (CLM5 ✅). Jennings et al. (2018) ✅ dan el umbral del 50 % observado en tierra en 1,0 °C de media, el centro de la rampa.
- **Energía de la nieve:** congelarse suelta L_f en la capa baja del aire (como el CLM5 ✅).
  - Sobre el océano libre, la nieve se funde con calor del océano.
  - Sobre el hielo marino, se suma a su masa.
  - En tierra, se acumula.

## 7. Suelo y nieve en tierra

- **Cubo de Manabe:** W_fc = 150 kg/m², β = min(1, W/(0,75·W_fc)) (Vallis et al. 2018, Isca, ✅ secundaria). Lo que sobra es escorrentía, que vuelve al océano al instante.
- **Nieve:** reserva S (kg/m² de agua).
  - Se funde con la energía del suelo por encima de 0 °C.
  - El agua de la fusión va al cubo.
- **Tope:** S_max = 1000 kg/m² (el del CLM4.5 ✅, vía `DISENO_FASE5A.md` §3.10). Lo que sobra va al océano como descarga glaciar y se funde a costa del océano (agua y energía conservadas).
- **Albedo gradual:** fracción cubierta S/(S + 10 kg/m²). La constante de 10 kg/m² es (c) ⚠️, sin fuente: unos 4 cm de nieve de 250 kg/m³ cubren media celda.
- **Fundiéndose:** el albedo baja linealmente en el último grado antes de 0 °C (CICE ✅).
- **Limitaciones:**
  - sin aislamiento térmico de la nieve (ni en tierra ni sobre el hielo marino; PENDIENTE desde la 5a);
  - sin envejecimiento de la nieve.

## 8. Suelo según su agua (I13) y albedo según la estrella (I14)

**Suelo:** fórmulas del CLM5 ✅ (Farouki 1981; de Vries 1963; número de Kersten de Johansen), con un suelo mineral medio (50 % arena, 50 % arcilla, porosidad 0,43).

| | Conductividad (W/m/K) | Capacidad (J/m³/K) | Inercia térmica |
|---|---|---|---|
| Seco | 0,22 | 1,29×10⁶ | ~530 |
| Saturado | 2,15 | 3,09×10⁶ | ~2580 |
| Antes (v2.4–v3.0, "roca densa") | 2,5 | 2,5×10⁶ | 2500 |

- La saturación se mide con el cubo, Sr = W/W_fc. Es una aproximación (c): el cubo representa la capa activa del suelo.
- Las propiedades se actualizan **una vez al año** con el agua media del año anterior, y en el equilibrio son las de la humedad climatológica.
- Al cambiar la capacidad, **la temperatura del suelo no cambia**: el agua entra y sale a la temperatura del suelo, y su calor sensible no se contabiliza (opción A).
  - Ocurre solo entre años, nunca en el año final donde se mide el cierre.
  - La primera versión conservaba C·(T − 0 °C). En suelos muy fríos, como la Antártida, eso daba saltos de decenas de grados e impedía converger: se detectó en la simulación de la Tierra.
- **Limitación:** no capta la variación estacional de la humedad del suelo.
- **Respuesta a la indicación de Carlos** ("tener en cuenta el agua para disminuirla"): el agua **sube** la inercia; un suelo seco es el que cambia de temperatura deprisa. El valor de base, el del suelo seco, sí es mucho menor que el de antes, como Carlos intuía.

**Albedo:**
- Por bandas, de CICE/Icepack ✅: nieve fría 0,98 en el visible y 0,70 en el infrarrojo cercano; fundiéndose, −0,10/−0,15; hielo desnudo, 0,78/0,36.
- El reparto entre bandas sale de la ley de Planck con la temperatura de la estrella (S3N 5420 K: 44,2 % de la luz por debajo de 0,7 µm; el Sol, 48,8 %).
- **Resultados:**
  - nieve fría 0,824 en P3N, frente a 0,837 con el Sol;
  - nieve fundiéndose 0,696, frente a 0,711;
  - el hielo marino observado en la Tierra (0,65, SHEBA) se multiplica por 0,966, y queda en 0,628.
- **Comprobación independiente:** Shields et al. (2013) ✅, con espectros completos, obtienen una reducción **mayor** (nieve 0,796 con una estrella G y 0,748 con una K). El cálculo por dos bandas es **conservador**. Limitación declarada.
- Con I14 apagado se usan los valores calculados con el Sol (la nieve también existe sin I14).

## 8 bis. Vapor radiativo (I15, PROTOTIPO, 05/10/2026) — pendiente de decisión de Carlos

**Motivo.** La validación de la v3.0 (`DISENO_V3.0.md` §13 bis): polos demasiado cálidos y sin hielo. Con un espesor óptico fijo, el aire polar, frío y seco, es tan opaco al infrarrojo como el tropical.

**Esquema:**
- d(τ) = (A + B·q)·dp/P0, la forma de Byrne y O'Gorman (2013), verificada en el código de Isca (esquema `byrne`, que con su modelo usa A = 0,8678 y B = 1997,9).
- Se recalcula en cada paso con el vapor de cada capa.
- **Calibración para M3N**, con el mismo método que el τ fijo de la v3.0 (`calibrar_tau_vapor.py`):
  - temperatura de la Atmósfera Estándar 1976 con el suelo a 289 K;
  - perfil de vapor con la forma de la atmósfera estándar de EE. UU. de la AFGL (Anderson et al. 1986), escalado a 24,9 kg/m² (Trenberth y Smith 2005 ✅);
  - objetivo: DLR = 314 y OLR = 267 (Wild et al. 2019).
  - Resultado: **A = 0,5656, B = 773,1**. Espesor total de la referencia: 2,45, el mismo que la v3.0.

**Columnas tipo de la AFGL** (temperatura y vapor de cada perfil):

| Columna | Suelo | Agua precipitable | τ fijo (v3.0): OLR / DLR | τ con vapor: OLR / DLR |
|---|---|---|---|---|
| Tropical | 299,7 K | 38,9 | 312 / 371 | 308 / 411 |
| Latitud media, verano | 294,2 K | 27,8 | 300 / 350 | 304 / 363 |
| Estándar de EE. UU. | 288,2 K | 13,6 | 265 / 314 | 284 / 258 |
| Subártico, invierno | 257,2 K | 4,1 | 203 / 220 | 214 / 124 |

**Lectura:**
- En lo alto de la atmósfera, el cambio en los polos es moderado: +11 W/m² de enfriamiento.
- **En la superficie, el cambio es enorme:** el infrarrojo que llega al suelo en el invierno subártico baja de 220 a 124 W/m². Es lo que en la realidad deja que la superficie polar se enfríe muchísimo en la noche polar y que se forme hielo, el arranque de la retroalimentación hielo-albedo.
- El diagnóstico de §13 bis se matiza: el exceso de calor polar viene sobre todo de la superficie sin hielo, que absorbe ~140 W/m² de luz (en la Tierra real, mucho menos). Lo que impide el hielo es el infrarrojo hacia el suelo, demasiado alto con el τ fijo.

**Limitaciones:**
1. Con solo dos parámetros y una sola columna de referencia, A y B no quedan bien separados. El DLR de las columnas secas probablemente sale **demasiado bajo**: en la realidad, el CO₂ y el continuo del vapor emiten mucho cerca del suelo incluso con aire seco.
   - No se ha podido contrastar con cálculos línea a línea de referencia: la página de Springer que los daba quedó bloqueada por límite de peticiones, y RRTMG no se instala con pip en este entorno.
   - **Pendiente:** un tercer dato verificado (columna seca) para separar A y B.
2. La absorción de luz solar por el vapor sigue fija (Wild 2019).

**Prueba en el modelo completo (modo Tierra, 05/10/2026):**
- Con I15 la Tierra simulada **se calienta sin freno**: 29,6 °C de media en la superficie el primer año y 31,6 °C el segundo, con los polos a +20 °C. El agua precipitable llega a 71 kg/m², frente a 24,9 observados.
- Sin I15, el mismo modelo da 20–22 °C con 38 kg/m².
- Se paró la simulación: no aportaba más.

**Causa:** la humedad relativa de la troposfera libre sale ~0,85 en todas partes. La convección relaja hacia RH_ref = 0,8, y no hay nada que seque el aire. En la realidad lo seca la **subsidencia** de la circulación (célula de Hadley, Fase 6), y la RH de la troposfera libre tropical ronda 0,3–0,5. Con tanto vapor y el vapor radiativo, la retroalimentación del vapor se dispara.

**Conclusión:** el vapor radiativo **no es utilizable** mientras la humedad del modelo no sea realista. Necesita, como mínimo, la sensibilidad a RH_ref (en curso, RH_ref = 0,6) y, de verdad, la circulación de la Fase 6.

**Estado:** programado y probado (conserva energía y agua), y **apagado**. No se recomienda encenderlo todavía. La decisión de adelantarlo es de Carlos.

## 9. Validación

### 9.1 Pruebas exactas (`test_v31.py`)

- Conservación en la convección húmeda y en la condensación.
- Columna seca y estable: la convección no hace nada.
- Saturación frente a valores de referencia.
- Ley de Planck frente a la integral numérica.
- Suelo seco y saturado.
- Rampa de la nieve.
- Interruptores apagados por defecto.
- Simulación corta con todo encendido: cierre de energía < 10⁻⁹, cierres de agua < 10⁻⁹ kg/m², sin recortes.

### 9.2 Modo Tierra y P3N

Ver §10, que se rellena con los resultados.

## 10. Resultados

Estado al entregar la v3.1-pre2 (05/10/2026, mediodía). Todas las simulaciones están hechas en el entorno de la IA (~30 min por año terrestre) y **ninguna ha llegado aún al equilibrio**:

| Modo Tierra, superficie | Año 1 | Año 2 | Año 3 |
|---|---|---|---|
| v3.1, RH_ref = 0,8: media global | 21,9 °C | 20,2 °C | 19,5 °C |
| — polos N / S | 11,5 / −1,2 | 10,0 / −3,5 | 9,2 / −4,6 |
| — agua precipitable | 38,1 kg/m² | — | 31,3 kg/m² |
| v3.1, RH_ref = 0,6: media global | 21,7 °C | | |
| — agua precipitable | 31,8 kg/m² | | |
| v3.1 + vapor radiativo (I15) | 29,6 °C | 31,6 °C | (parada) |

Lectura provisional:
- El modelo se va enfriando hacia ~19 °C, **sin hielo marino** en ningún hemisferio. El problema de los polos de la v3.0 sigue con el agua.
- La humedad relativa de la troposfera libre (~0,83 con RH_ref = 0,8) es demasiado alta, porque falta la subsidencia de la Fase 6. Con RH_ref = 0,6 baja a ~0,58 y el agua precipitable se acerca a la observada.
- **Pendiente:** la calibración de D con el vapor y la validación completa (`validar_v30.py --v31`), mejor en el PC de Carlos.

**Aviso:** la validación de la v3.0 (`DISENO_V3.0.md` §13 bis) ya mostró que, con un espesor óptico fijo, los polos salen demasiado cálidos y sin hielo. La v3.1 no lo arregla, y el vapor radiativo, con la humedad actual, lo empeora (§8 bis).

**Punto de control (05/10/2026):** `simular_fase2b(..., archivo_estado=ruta)` guarda el estado completo al acabar cada año y, si se relanza con el mismo archivo y los mismos parámetros, continúa donde se quedó. El resultado es el mismo bit a bit que sin interrupción (`test_v31.py`).

## 11. Pendiente

- Recalibrar D (y D_oc) con el vapor en el modo Tierra frente a Trenberth y Caron (2001).
- Sensibilidad: RH_ref (0,6 / 0,8), RH crítica (0,9 / 1), paso de tiempo.
- Viento medio oceánico verificado para el modo Tierra.
- Aislamiento por la nieve (tierra y hielo marino).
- Opción C del calor latente (entalpía completa del agua).
- Exportación a H3N de la precipitación, la nieve y la humedad del suelo (claves nuevas en `m3n-clima`).
- Lluvia tropical **no válida** hasta la Fase 6 (sin célula de Hadley), como estaba declarado.

### 10 bis. Simulaciones largas del modo Tierra con la v3.1 completa (I11–I14 encendidos, I15 apagado), 05/10

Se hicieron en el entorno de la IA, que se reinició varias veces; se retomaron siempre desde el último año guardado.

| Año | T sup. global, RH_ref = 0,8 | T sup. global, RH_ref = 0,6 |
|---|---|---|
| 1 | 21,90 °C | 21,71 °C |
| 2 | 20,15 | 19,89 |
| 3 | 19,48 | 19,22 |
| 4 | 19,21 | 18,96 |
| 5 | 19,10 | (pausada) |
| 9 | 19,00 | — |
| 10 | 18,99 | — |

- **RH_ref = 0,8: estabilizada en ~19,0 °C a partir del año 9**, con los polos a +8,5 °C (norte) y −6,2 °C (sur). La corrida con RH_ref = 0,6 se perdió en el año 4 por un reinicio; no se retoma.
- **Polos:** el norte en torno a +9 °C y el sur entre −5 y −6 °C. **Sin hielo marino** en ninguno.
- **Sensibilidad a RH_ref:** de 0,8 a 0,6 solo cambia ~0,25 K. La humedad excesiva **no** se corrige con ese parámetro.
- **Coherente con el diagnóstico del §8 bis:** falta la circulación, que es la Fase 6 (camino B, decidido el 05/10).
- Las corridas no se continuarán: la v3.1 queda apagada y la Fase 6 la sustituirá como base de la validación.

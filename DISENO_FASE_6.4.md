# DISENO_FASE_6.4.md — Fase 6.4: P3N con el núcleo dinámico (clima oficial 4.0.0)

**Estado:** en curso. La 6.3 se cerró el 11/10/2026 (`DISENO_FASE_6.3.md` §6.25). Lo decidido está en el §1 (antes de la 6.4) y en el §7 (decisiones del 11/10). **La simulación definitiva de P3N está hecha y analizada (12/10/2026, §8)**; queda investigar el desplazamiento al polo de los vientos del oeste (§8.4), las sensibilidades del océano (§4.2), la exportación a H3N y la aprobación de la 4.0.0.
**Convención:** ✅ = verificado en la fuente citada; ⚠️ = por verificar; 🔶 = propuesta de la IA pendiente de Carlos; (c) = valor calibrado en la Tierra.

## 0. Qué es la Fase 6.4 y qué entrega

- **Qué es:** llevar a P3N el modelo validado en modo Tierra en la 6.3: núcleo dinámico (I16) con toda la física de la 6.3 y la configuración de `clima_dinamico.py`.
- **Qué entrega:**
  1. la climatología de 30 años de P3N (mapa activo) con el núcleo;
  2. su análisis frente al clima oficial actual (2.4.3, modelo de 2 capas) y frente a la v3.1 sin núcleo;
  3. si cumple el criterio que se apruebe antes (§4.7), **la 4.0.0**: el clima oficial que M3N exporta para H3N pasa a ser el del núcleo. Es provisional hasta tener nubes (decisión del 09/10/2026).
- **Cuándo empieza:** cuando la 6.3 esté cerrada con su criterio (`DISENO_FASE_6.3.md` §6.22). Este diseño se prepara antes para que las decisiones estén tomadas a tiempo.

## 1. Lo que ya está decidido (no se reabre aquí)

| Decisión | Fecha | Dónde |
|---|---|---|
| El clima oficial pasa al núcleo con la 4.0.0, tras la 6.4, sin esperar a las nubes (provisional) | 09/10/2026 | `DISENO_FASE_6.3.md` §6.17 |
| Climatología de 30 años (normal de la OMM) y prueba de deriva a posteriori | 09/10/2026 | §6.17 |
| La estrella no se recalibra hasta tener nubes (Fase 5.2) | 05/10/2026 | `DISENO_FASE_2.3.md` §13 |
| Rejilla de 5°; la de 2,5° se decide con el planeta definido | 06/10/2026 | `DISENO_FASE_6.3.md` §6.3 |
| Mapa de trabajo: pizarra1 | 06/10/2026 | briefing del 06/10 |
| Día solar de 19,84 h (71 424 s); hora de P3N = 1/24 del día | 04/10/2026 | `parametros.py` |
| 1 bar al nivel del mar, hasta el barrido final (Fase 9) | 01/10/2026 | `fase30_multicapa.P_NIVEL_MAR` |
| Principio: no asumir valores de la Tierra para P3N | 04/10/2026 | `DISENO_FASE_6.3.md` §6.16 |
| Modelo propio; nunca ajustar parámetros para esconder un error | permanente | — |

## 2. P3N tal como lo ve el modelo (leído del código de la 3.18.1 y del mapa activo, 11/10/2026)

| Magnitud | P3N | P3N / Tierra | Origen |
|---|---|---|---|
| Radio | 4676,7 km | 0,734 | `parametros.P3N_CIRCUNFERENCIA` |
| Gravedad en superficie | 9,053 m/s² | 0,923 | de la masa y el radio |
| Día solar | 71 424 s (19,84 h) | 0,827 | decisión del 04/10 |
| Velocidad de giro sideral Ω | 8,824·10⁻⁵ rad/s | 1,210 | `parametros.OMEGA_SIDERAL` |
| Año | 23 337 580 s = 326,75 días de P3N = 270,1 días terrestres | 0,739 | `ORBITA_PERIODO` |
| Excentricidad | 0,046 | 2,75 | `ORBITA_EXCENTRICIDAD` |
| Inclinación del eje | 1,7° | 0,073 | `INCLINACION_AXIAL` (§4.1) |
| Estrella S3N | 0,874 M☉, 5420 K, 0,627 L☉ | — | `S3N_MASAS_SOLARES` (calibrada el 02/10 para el modelo de 2 capas) |
| Irradiancia a la distancia media (a = 0,782 UA) | 1395,5 W/m² | 1,025 | cálculo de `parametros` |
| Irradiancia en el periastro / el apoastro | 1533 / 1275 W/m² (cociente 1,202) | — | (1 ± e)⁻² |
| Paso de la física / del núcleo | 992 s (72 por día) / 496 s | 1,10 | `temperatura.PASO_TIEMPO` |
| Pasos por año | 23 526 | 0,67 | — |
| Celda en el ecuador / a 87,5° | 408 km / 17,8 km de ancho | 0,734 | rejilla de 5° |
| Mapa pizarra1: tierra | 39,6 % del área (norte 40,8 %, sur 38,4 %) | — | puente de C3N (04/10) |
| Mapa pizarra1: altitud | **0 m en todas las celdas** | — | puente de C3N (§4.6) |

## 3. Qué circulación cabe esperar

- **Escalas** (cálculo propio con la tabla del §2; ya en `DISENO_FASE_6.md` §1):
  - número de Rossby térmico gH/(Ω²a²), con el mismo contraste de temperatura y la misma altura: **×1,17**;
  - anchura de la célula de Hadley ∝ (número de Rossby térmico)^½ (Held y Hou 1980, *J. Atmos. Sci.* 37, 515 ✅ como proporcionalidad, citada así en Guendelman y Kaspi 2019, *Geophys. Res. Lett.*, doi:10.1029/2018GL080752 ✅): **×1,08**, es decir, ~2–3° más ancha en latitud que en el modo Tierra con la misma física;
  - β = 2Ω/a: ×1,65; radio de deformación / radio del planeta: ×1,13.
  - Lectura: tres células por hemisferio, como en la Tierra. Es una expectativa de orden de magnitud, no una predicción: la teoría de Held y Hou es para una atmósfera sin remolinos.
- **Estaciones con 1,7° de inclinación** (cálculo propio, luz media diaria en el tope, fórmula estándar de la insolación diaria):
  - a 45°, la inclinación hace oscilar la luz entre 299 y 329 W/m² (30 W/m²);
  - la excentricidad (0,046) la hace oscilar entre el periastro y el apoastro en ~58 W/m² a esa latitud, y en los **dos hemisferios a la vez**.
  - Por tanto, con 1,7° el ciclo anual de P3N lo manda la **distancia a la estrella**, no la inclinación: un "verano" global cerca del periastro y un "invierno" global cerca del apoastro, con la banda de lluvias casi quieta en el ecuador.
  - Linsenmeier, Pascale y Lucarini (*Planet. Space Sci.*, doi:10.1016/j.pss.2014.11.003; arXiv:1401.5323 ✅): con poca inclinación el ecuador recibe más luz que los polos y el contraste de temperatura entre ellos es marcado; las desviaciones estacionales grandes son propias de las inclinaciones altas.
- **Más cálido que el clima oficial actual (importante):** la estrella se calibró el 02/10 para que el modelo de **2 capas** diera ~15–16 °C en el aire a 2 m. El modelo con núcleo no tiene nubes (Fase 5.2): en modo Tierra sale +5 K por encima de la Tierra real (§6.14 de la 6.3). Lo esperable es que P3N salga también **bastante más cálido** que los 14,5 °C del clima oficial de hoy (2.4.3, pizarra1). No se corrige: es la decisión del 05/10 (la estrella espera a las nubes). Pero **la 4.0.0 cambiará el clima que se ve en H3N**, y conviene saberlo antes de aprobarla.

## 4. Cuestiones de diseño

### 4.1 La inclinación del eje (worldbuilding: Carlos y Ozan)

- **Hoy:** 1,7°, un valor provisional. Ya se anotó que conviene fijarla antes de la 6.4.
- **Por qué importa ahora:** cambia el clima entero (estaciones, banda de lluvias, hielo, monzones) y, por tanto, la simulación larga de la 6.4 y la 4.0.0. Con otro valor habría que repetirla (~6,5 h).
- **Qué cambia según el valor** (luz media diaria en el tope a la distancia media, W/m², cálculo propio; la del equinoccio no depende de la inclinación):

| Inclinación | Ecuador, solsticio | 45°: verano / invierno | 80°: verano / invierno |
|---|---|---|---|
| 1,7° (hoy) | 444 | 329 / 299 | 99 / 58 |
| 10° | 437 | 400 / 228 | 239 / 0 (noche polar) |
| 23,44° (la Tierra) | 408 | 512 / 119 | 547 / 0 |
| 35° | 364 | 606 / 40 | 788 / 0 |

  (Equinoccio: 444 en el ecuador, 314 a 45°, 77 a 80°. Con la excentricidad, cada valor sube hasta un ~10 % en el periastro y baja otro tanto en el apoastro.)
- **M3N admite cualquier valor** (`parametros.INCLINACION_AXIAL`); el calendario (día 1 = solsticio de invierno del norte) se recalcula solo.
- **Decisión de worldbuilding, no técnica:** la IA no propone un valor. Lo técnico es el **cuándo**: antes de la simulación definitiva de la 6.4. Si no está decidida para entonces, hay dos caminos: (a) esperar; (b) hacer la 6.4 con 1,7° y repetirla cuando se decida (otra noche).

### 4.2 El océano: capa de mezcla y difusión (el único transporte del océano hasta la 6.5)

**Lo que hay (leído en el código y los documentos):**
- **Capa de mezcla fija de ~48,7 m.** No viene de una fuente: es la capacidad que daba una fórmula heredada de la Fase 0 (inercia·√(día/π)) con un día de 24 h, y se fijó en ese valor el 05/10 al quitarle la dependencia del día (`DISENO_FASE_2.3.md` §14).
- **Difusión del océano D_oc = 0,27 (c)** (`D_OCEANO_V30`): calibrada en la Fase 2.3 en modo Tierra frente a Trenberth y Caron (2001), **junto con una atmósfera difusiva** (D_atm = 1,55). Con I16 la atmósfera ya no difunde, pero esta D_oc sigue activa. No se escala con la rotación.
- **Lo que dice el modo Tierra con I16** (3.17.0): el océano lleva la parte correcta del transporte a 35° N (0,19 frente a 0,22), demasiada a 35° S (0,34 frente a 0,08) y mucha a 50° S (1,45 PW): una difusión no tiene la barrera de la corriente circumpolar antártica. Es un problema de **forma**, no de un coeficiente.
- **El océano real** (de Boyer Montégut et al. 2004, *J. Geophys. Res.* 109, C12003, doi:10.1029/2004JC002378 ✅): la capa de mezcla va de menos de 20 m en verano a 150 m en invierno en subtrópicos y latitudes medias (~70 m de media anual en latitudes medias), más de 300 m en el invierno del océano Austral y 500–740 m en las zonas de convección profunda. Es una propiedad **dinámica** (viento, estratificación, estaciones), no una constante del agua.

**El problema de principio:** las dos cosas son propiedades de un océano que M3N todavía no tiene (6.5), y dependen de su viento, su rotación y la forma de sus cuencas. El principio del 04/10 dice que no se asuman valores de la Tierra; pero sin modelo de océano, cualquier valor es una suposición.

**Opciones:**
- **(A) Mantener 48,7 m y D_oc = 0,27 como valores provisionales declarados** (c), los mismos con los que se validó el modelo en modo Tierra, y **medir cuánto dependen de ellos** los resultados de P3N con pruebas de sensibilidad (abajo). Coherente con la 6.3, honesto con la incertidumbre.
- **(B) Recalibrar D_oc en modo Tierra con I16** frente a Trenberth y Caron. Sería más coherente con la atmósfera de ahora, pero sigue siendo un valor terrestre, y un solo coeficiente no arregla la forma del transporte (punto anterior). Coste: varias simulaciones del modo Tierra.
- **(C) Escalar con P3N** (por ejemplo, con la rotación): no hay una ley verificada para una difusión de océano equivalente; sería inventar.

**Decisión técnica (11/10/2026; Carlos la delega en el criterio técnico): (A).**
- **Cuatro climatologías exploratorias de P3N** (12 años cada una, WMO-No. 1203 §5.2.3, como ya permite `--anos`), cambiando un solo valor cada vez, a la mitad y al doble: D_oc × 0,5 y × 2; capa de mezcla 24,4 m y 97,4 m.
- El factor 2 en los dos sentidos es simétrico y cubre el rango observado de la capa de mezcla fuera de las zonas de convección profunda (de 20 a 150 m en subtrópicos y latitudes medias; de Boyer Montégut et al. 2004).
- Coste: ~3 h cada una en el PC de Carlos (~12 h, dos noches o en paralelo si el PC lo permite).
- El resultado se escribe en la exportación de la 4.0.0 y en este documento como incertidumbre del clima oficial ("tal región, ±x K según el océano").
- **Hecho (11/10/2026):** `clima_dinamico.py --anos 12 --oceano-difusion F` y `--oceano-capa-mezcla F` (F = 0,5 o 2).
  - La capa de mezcla se cambia con la inercia del agua, de la que sale su capacidad (`fase2_inercia_multicapa.py`); también la usa el modelo de partida del arranque.
  - D_oc pasa al modelo con núcleo con `fase6_clima.simular_clima(..., opciones_nucleo={"d_oceano": ...})`. El modelo de partida conserva la suya, porque solo es el punto de partida y el equilibrio lo decide el modelo con núcleo.
  - Cada prueba va a su carpeta (`p3n_dos_capas_exploratoria12_oceanoD2`, `..._capa0.5`) con `sensibilidad.json`. No se puede continuar una simulación con otro océano: una carpeta con simulación y sin ese archivo cuenta como factores 1 y 1.
  - Los factores van en `resultado.pkl`, en el resumen y en la exportación (`climatologia.sensibilidad`). Esa exportación no puede ser la que lee H3N (`--activo`).
  - Con los factores en 1, el resultado es idéntico bit a bit al de antes (comprobado con `--dias 1`).
  - 2 pruebas (`test_sensibilidad_oceano.py`).
  - Órdenes, una por prueba: `python clima_dinamico.py --anos 12 --oceano-difusion 0.5`, `... --oceano-difusion 2`, `... --oceano-capa-mezcla 0.5` y `... --oceano-capa-mezcla 2`. Después, `python analizar_p3n.py --comparar <las cuatro carpetas>`.

### 4.3 El arranque de P3N

- **Lo que hay:** la cadena arranca del equilibrio del modelo de 2 capas (§6.10.2 de la 6.3). En el mapa de pruebas antiguo (continente rectangular de 600 m) ese equilibrio tenía el ecuador a ~45 °C, ~14 K más que la v3.1. No se pudo comprobar si I16 arranca bien desde ahí; existe la alternativa `--arranque v31`.
- **Prueba en el entorno de la IA (11/10/2026):** P3N con pizarra1, la 3.18.1 y la configuración real (órbita completa, arranque de 2 capas, `clima_dinamico.py` sin `--dias`). Resultado:
  - **Modelo de partida:** equilibrio en **23 años** (cambio máximo de una celda 0,0116 K < 0,015), como el clima oficial de pizarra1 del 05/10 (23 años). Del año 9 al 18 el cambio máximo se quedó entre 0,9 y 4,7 K en celdas de mar a 67,5° N y 67,5° S, es decir, en el **borde del hielo marino** (comprobado con `DEPURAR`: la celda del cambio máximo y los saltos del hielo grueso de los años 6, 11 y 17). Del año 19 al 23 bajó solo a 0,19, 0,037, 0,024, 0,018 y 0,012 K. No es un fallo: es el tiempo que tarda en asentarse la extensión del hielo.
  - **Núcleo con P3N, por primera vez:** arranca sin inestabilidades y sigue estable durante el 72 % del primer año (~236 días de P3N). Ahí se paró a propósito, para no gastar más: el clima de verdad lo da la simulación en el PC de Carlos.
  - **Pruebas cortas (`--dias 1` y `--dias 3`):** con P3N se rompen en las primeras horas del núcleo (Courant vertical de los trazadores 1,3 en una columna a 7,5° S, por un ascenso de ~3·10¹⁴ kg por subpaso en la troposfera media). Con la Tierra, `--dias 1` pasa (`--dias 3` también se rompe). Esas pruebas parten de un estado irreal: un modelo de partida de 2 "años" sin equilibrar, calentado con la luz de un solo día repetido (el día 0 de P3N está junto al periastro, 1524 W/m²). Ya estaba documentado en la 6.3 que esos estados están al borde de la inestabilidad (§6.21). La prueba que vale es la de la configuración real. Consecuencia práctica: **las pruebas automáticas de P3N con I16 no pueden usar `--dias 1`**; usarán el modo Tierra, como hasta ahora, o un estado de P3N ya equilibrado.

### 4.4 La composición de la atmósfera (supuesto declarado)

- El infrarrojo gris (τ₀, f) y la luz absorbida están calibrados con la **columna media de la Tierra** (Wild et al. 2019; desde la 3.18.0, con su masa observada). En P3N se aplican escalados con su masa de aire (1 bar al nivel del mar).
- Eso supone **los mismos gases que la Tierra** (CO₂, vapor, oxígeno y ozono en las mismas proporciones). Está decidido que P3N tiene oxígeno (04/10); la cantidad de CO₂ no.
- **Propuesta 🔶:** declararlo como supuesto del clima oficial hasta el barrido final (Fase 9), sin tocar nada.

### 4.5 El "1 bar" de P3N: ¿total o de aire seco?

- Pregunta abierta desde la 3.18.0 (§6.21). Hoy es de aire **seco**; el vapor añadiría 2–3 hPa al nivel del mar.
- Cambia muy poco el clima. Conviene decidirlo antes de la 4.0.0 para que el número oficial signifique lo que se quiere.
- **Decidido (11/10/2026): aire seco (§7.4).**

### 4.6 El mapa: pizarra1 no tiene relieve

- En el puente que lee M3N (`mapa_activo_m3n.json`, del 04/10), **todas las celdas están a 0 m**. pizarra1 es un mapa de trabajo de tierra y mar sin montañas.
- Consecuencias: (1) la 4.0.0 sería el clima de un P3N **llano**: sin lluvia orográfica ni sombras de lluvia, y con los desiertos donde los pongan la circulación y la forma de los continentes; (2) la estabilidad del núcleo sobre relieve (los problemas antárticos de la 6.3) no se pone a prueba en P3N.
- **Decidido (11/10/2026):** falta pintar el relieve en C3N; la 6.4 se hace con pizarra1 llano (§7.5).

### 4.7 Cómo se juzga P3N sin observaciones (criterio para la 4.0.0; aprobado por Carlos el 11/10/2026)

*(Rehecho el 11/10 tras la objeción de Carlos: la primera versión exigía tres células por hemisferio y la banda de lluvias cerca del ecuador. Era un error de planteamiento: eso son **resultados** que dependen de la rotación, el tamaño del planeta y la inclinación, y no se pueden exigir. En simulaciones idealizadas, con más rotación las células de Hadley se encogen, aparecen más chorros y células de Ferrel, y crece el contraste de temperatura entre el ecuador y los polos (Kaspi y Showman 2015, *Astrophys. J.* 804, 60 ✅). P3N gira 1,21 veces más rápido que la Tierra pero es más pequeño, y la teoría de Held y Hou da para él células algo más anchas, no más estrechas (§3): lo previsible es una estructura parecida a la terrestre, pero es una previsión, no un requisito.)*

P3N no tiene observaciones. Lo que se puede exigir es que **el modelo funcione igual de bien que en el modo Tierra** (lo validado en la 6.3 es el código, no el planeta) y que **lo que salga se entienda**. Propuesta, fijada antes de la simulación:
- **A. Bloqueantes, los mismos de la 6.3** (`DISENO_FASE_6.3.md` §6.22):
  1. toda la simulación con una sola versión del código;
  2. sin inestabilidades;
  3. equilibrio alcanzado;
  4. energía y agua cerradas cada año (energía < 10⁻⁷ en términos relativos, agua < 10⁻¹¹ kg/m²);
  5. el error numérico de la dinámica, devuelto por el corrector y sin tendencia;
  6. sin deriva en los 30 años;
  7. el análisis y la exportación a H3N terminan sin errores.
- **B. Comprensión de lo que salga** (sin resultados exigidos):
  - el análisis describe la circulación (cuántas células y dónde, chorros, alisios, banda de lluvias, desiertos, hielo, transporte de calor) y la compara con lo que da la teoría para la rotación y el tamaño de P3N (§3) y con el modo Tierra de la misma versión;
  - cada diferencia grande con lo esperado se explica con fuente, o se investiga antes de la 4.0.0. Nada se rechaza por ser distinto de la Tierra.
- **C. Comparación documentada** con el clima oficial actual (2.4.3) y con la v3.1 sin núcleo: temperatura, hielo, lluvia y vientos, con las diferencias grandes explicadas (la mayor se sabe ya: sin nubes, más calor).
- **D. Sensibilidad al océano** (§4.2): cuánto cambia el clima con los valores del océano a la mitad y al doble. Se escribe como incertidumbre del clima oficial; no aprueba ni suspende.
- **Carlos aprueba la 4.0.0.**

### 4.8 La 4.0.0: qué cambia en la exportación y en H3N

- **Mecánica (ya programada en la 3.16.0):** `python exportar_clima.py --dinamico --activo` escribe la climatología del núcleo en el archivo oficial (`clima_activo_m3n.json`) en vez de en `clima_dinamico_m3n.json`.
- **La 4.0.0** es el paquete que convierte eso en lo normal: la exportación por defecto pasa a ser la del núcleo, `CHANGELOG` y documentos lo dicen, y el clima de 2 capas sigue disponible.
- **Consecuencia de flujo de trabajo (importante) 🔶:** hoy, tras editar el mapa en C3N, `exportar_clima.py` tarda ~6 min (modelo de 2 capas). Con el núcleo, un clima nuevo son **~6,5 h** (una noche). Propuesta: mantener el modelo de 2 capas como **vista rápida** en H3N (marcada así) mientras llega el clima del núcleo, y que el indicador de estado diga cuál se está viendo. Decisión de Carlos (11/10): **idea anotada, no prioritaria**; se diseñará aparte, en H3N (ver §7.7).
- **H3N:** con la 4.0.0, el clima oficial ya trae lluvia, viento y humedad. Hay que revisar el selector de fuente de la v0.2.2: pasaría a "oficial (núcleo)" frente a "vista rápida (2 capas)".

### 4.9 Coste

- Con la velocidad medida en el PC de Carlos en el modo Tierra (~40 pasos/s), un año de P3N (23 526 pasos) son ~10 min. Con ~10 años de equilibrio y 30 de climatología: **~6,5 h**, más el arranque (minutos). La validación se adapta (§5).
- Las cuatro sensibilidades del océano (§4.2): climatologías exploratorias de 12 años, ~3 h cada una (~12 h en total).

### 4.10 Comprobaciones técnicas propias de P3N

- **Paso de tiempo:** el subpaso del núcleo (496 s) frente a la celda ecuatorial (408 km) da un número de Courant de la advección ~1,5 veces el del modo Tierra (450 s y 556 km) para los mismos vientos. El paso semiimplícito trata las ondas de gravedad, y cerca de los polos manda el filtro polar. Lo comprueba la prueba del §4.3.
- **Hiperdifusión:** τ = 0,5 "días" son 0,5 × 86 400 s (`fase6_nucleo.preparar_hiperdifusion`), con la celda ecuatorial del propio planeta: está bien definida en P3N (amortigua la onda más corta en el mismo tiempo físico).
- **Coeficientes del modelo de 2 capas escalados con la rotación** (`FACTOR_ROTACION_D` = 0,683): con I16 la atmósfera no difunde y no se usan, salvo en el arranque (modelo de 2 capas).

## 5. Plan de trabajo propuesto 🔶

1. Decisiones de Carlos (§6), con Ozan para la inclinación.
2. Prueba de arranque de P3N (§4.3). **Hecha (11/10/2026).**
3. Una herramienta de análisis de P3N: `validar_i16.py` sin las referencias terrestres, más lo que pide §4.7 (células, vientos, banda de lluvias, transporte, hielo, desiertos), y la comparación con la 2.4.3 y la v3.1. Con pruebas.
   - **Hecho (11/10/2026): `analizar_p3n.py`** (§5.1), con 8 pruebas (`test_analizar_p3n.py`). La comparación con la v3.1 sin núcleo no está: la v3.1 no se exporta y habría que simularla aparte (~1 h en el PC); se decide al llegar a ese punto.
4. Los cambios que pidan las decisiones (inclinación, "1 bar", relieve), cada uno con su versión.
5. La simulación definitiva de P3N (PC de Carlos, ~6,5 h), y las sensibilidades del océano si se aprueban. **La simulación, hecha y analizada (12/10/2026, §8)**; las sensibilidades, aprobadas (§7.2) y pendientes.
6. Aplicar el criterio del §4.7 → Carlos aprueba → paquete de la **4.0.0** y versión de H3N.

### 5.1 La herramienta de análisis: `analizar_p3n.py` (11/10/2026)

Uso, después de `python clima_dinamico.py` con el mapa activo: `python analizar_p3n.py` (opciones `--referencia-tierra`, `--clima-oficial`, `--comparar`, `--sin-ano-extra`). Escribe `outputs/analisis_p3n/informe.txt` y `campos.npz`. No cambia nada del modelo. Reutiliza las funciones de `validar_i16.py` y de `herramientas/criterio_cierre_6_3.py`, para que P3N se mida igual que la Tierra.
- **A:** los bloqueantes de la 6.3, leídos del `resumen.txt` con las mismas funciones del criterio de cierre (A7, la exportación y H3N, a mano).
- **B, sin suponer la forma del clima** (lo que pide el §4.7 corregido):
  - **todas** las células de la función de corriente, como tramos de signo a σ 0,5 con su extremo y sus bordes (TropD, Adam et al. 2018), en lugar de las dos células de `validar_i16.py`, que se buscan en ventanas pensadas para la Tierra (20° S–40° N y 40° S–20° N);
  - **todos** los núcleos de chorro (criterio de §6.16 de la 6.3, del oeste y del este);
  - **todas** las bandas de viento de la capa baja y todos los máximos y mínimos de la lluvia zonal;
  - lo que es menos del 5 % de lo mayor se marca "débil", pero no se oculta.
- **Estaciones de P3N:** como las manda la distancia a la estrella (§3), el año se mide desde el periastro: el periastro y el apoastro son el máximo y el mínimo de la luz media diaria en el tope, el ciclo anual va en doceavos de año orbital y hay dos temporadas de un cuarto de año centradas en el periastro y en el apoastro. Para eso `validar_i16.ano_extra` acepta otro calendario; por defecto sigue usando los meses del modo Tierra, y su prueba pasa igual.
- **B':** la anchura de las células tropicales del modo Tierra de la misma versión (`--referencia-tierra`, el `campos.npz` de `validar_i16.py`) multiplicada por el factor de Held y Hou para P3N (1,082, solo g, Ω y a; §3). Es una orientación, no un requisito.
- **C:** comparación con el clima oficial de 2 capas (el json que lee H3N), solo si es del mismo mapa (huella).
- **D:** diferencias con otras simulaciones de P3N del mismo mapa (las sensibilidades del océano).
- **Probado:**
  - las piezas, con casos construidos a mano (células alternas, una célula débil, varios chorros y uno del este, bandas de viento, máximos en los bordes, ceros exactos, el calendario con la órbita real: periastro y apoastro a medio año y luz en la proporción ((1 + e)/(1 − e))²);
  - el conjunto, con los datos de una prueba corta.
  - Falta probarla con una simulación de P3N de verdad (§4.3).

## 6. Decisiones de Carlos

*(Todas tomadas el 11/10/2026: ver §7. La lista queda como registro de lo que se preguntó.)*

1. **Inclinación del eje** (§4.1): el valor (con Ozan) y si se espera a tenerlo o se hace la 6.4 con 1,7°.
2. **Océano** (§4.2): (A) valores provisionales declarados más sensibilidad (recomendada), (B) recalibrar D_oc en modo Tierra o (C) otra. ¿Con sensibilidad de la capa de mezcla también?
3. **Composición de la atmósfera** (§4.4): ¿se declara como supuesto (Tierra) hasta la Fase 9?
4. **"1 bar"** (§4.5): ¿total o de aire seco?
5. **Relieve de pizarra1** (§4.6): ¿llano a propósito o falta pintarlo?
6. **Criterio para la 4.0.0** (§4.7): ¿de acuerdo, o cambios?
7. **Vista rápida de 2 capas en H3N** tras la 4.0.0 (§4.8): ¿sí o no?

## 7. Decisiones tomadas (11/10/2026)

1. **Inclinación del eje:** se queda en 1,7° durante la 6.4. Ahora importa más que el modelo sea preciso que el worldbuilding; si cambia, será al final (decisión de Carlos).
2. **Océano:** opción (A) con las cuatro sensibilidades del §4.2 (Carlos delega en el criterio técnico).
3. **Composición del aire:** como en la Tierra mientras no haga falta otra cosa (Carlos). Comprobado que hoy no hay ninguna incoherencia que obligue a decidirla: el aire seco (R, c_p) y las calibraciones de la radiación ya son los de la Tierra, y el oxígeno de P3N (04/10) es coherente con ellos. Volverá a importar en la Fase 5.2: la radiación por bandas y el ozono dependen de la composición.
4. **"1 bar" = aire seco** (Carlos delega en el criterio técnico). Motivo: la masa de aire seco es una propiedad fija del planeta (no aparece ni desaparece) y es lo que conserva el núcleo; el vapor es un resultado del clima que cambia con él. Fijar la presión **total** en 1 bar obligaría a ajustar el aire seco al vapor que saliera, de forma circular. Es lo mismo que se hizo en el modo Tierra con la masa observada (3.18.0). La presión total al nivel del mar se informa como diagnóstico (unos hPa por encima de 1 bar).
5. **Relieve:** falta pintarlo en C3N. La 6.4 se hace con pizarra1 llano. Cuando se pinte: una prueba corta del núcleo con el relieve nuevo y una simulación nueva del clima (una noche); no hace falta versión nueva de M3N (cambia el mapa, no el modelo).
6. **Criterio para la 4.0.0:** rehecho (§4.7) sin resultados exigidos (la estructura depende de la rotación y el tamaño; Kaspi y Showman 2015). **Aprobado por Carlos (11/10/2026).**
7. **Vista rápida en H3N** (§4.8): a Carlos le parece bien la idea, pero **queda anotada como idea, no prioritaria y fuera del alcance de la 6.4**: antes hay que pensar bien cómo funcionaría e integrarla en H3N. La 4.0.0 no depende de ella (el selector de fuente de la v0.2.2 sigue funcionando).

## 8. La simulación definitiva de P3N (12/10/2026)

### 8.1 Qué se simuló y cuánto costó (PC de Carlos)

- `python clima_dinamico.py --carpeta outputs/clima_dinamico/p3n_3.19.0` con la 3.19.0 (código `8154dd9b484f222cba1e671b`), mapa pizarra1 (huella `851bc12fc081346f`), arranque de 2 capas. Empezó el 11/10 a las 21:30.
- **Coste: 7,85 h** (estimado: ~6,5 h). La diferencia es el equilibrio: **15 años** en lugar de los ~10 supuestos. Dentro del techo de una noche (06/10). Velocidad del año extra del análisis: 36 pasos/s.
- Análisis: `python analizar_p3n.py --carpeta outputs/clima_dinamico/p3n_3.19.0 --referencia-tierra outputs/validacion_i16_3.18.1` (12/10, 05:37).

### 8.2 Criterio §4.7, parte A: los bloqueantes

| | Criterio | Resultado |
|---|---|---|
| A1 | un solo código | ✅ un tramo, 3.19.0 |
| A2 | estable, cadena completa | ✅ 30 de 30 años |
| A3 | equilibrio | ✅ en 15 años (N +0,132 W/m², T2m +0,0032 K/año, hielo −0,00016, agua del suelo +0,0037) |
| A4 | energía < 1·10⁻⁷ y agua < 1·10⁻¹¹ kg/m² cada año | ✅ máximos 4,5·10⁻⁸ y 3,2·10⁻¹³ |
| A5 | error de la dinámica controlado | ✅ media −0,047 W/m² (lo devuelve el corrector), sin tendencia (p 0,81) |
| A6 | sin deriva | ✅ −0,0002 ± 0,0011 K/año (p 0,87); ninguna banda; precipitación sin deriva (p 0,36) |
| A7 | análisis y exportación | análisis ✅; exportación y lectura en H3N, pendientes |

**El error de la dinámica es ~20 veces menor que en el modo Tierra** (−0,047 frente a −1,04 W/m²). Hipótesis: pizarra1 no tiene relieve, y en el modo Tierra buena parte de ese error viene de la fuerza de presión sobre las montañas (los problemas antárticos de §6.11–§6.12 de la 6.3). Si es así, crecerá al pintar el relieve (§7.5).

### 8.3 Parte B: el clima de P3N (medido; la lectura, marcada como tal)

**Lo que se entiende:**
- **Balance:** luz que llega 349,2 W/m², igual a S₀/4·(1 − e²)^−½ = 349,25 ✅; albedo 0,144 (sin nubes); neto en el tope −0,04 W/m². Precipitación = evaporación = 3,483 mm/día (modo Tierra 3,47).
- **Temperatura:** aire a 2 m 21,57 °C (norte 21,06, sur 22,08); +7,1 K sobre el clima oficial de 2 capas (14,45 °C). Era lo esperado (§3): sin nubes y con la estrella calibrada para el modelo de 2 capas; no se corrige hasta la 5.2 (decisión del 05/10). El ecuador queda igual que en el modelo de 2 capas y los polos, 23–38 K más templados: el transporte por el viento, lo que la Fase 6 tenía que aportar (como en el modo Tierra, §6.15 de la 6.3).
- **Agua precipitable 44,0 kg/m²**: la del modo Tierra (36,6) escalada con Clausius-Clapeyron (~7 %/K) por ~2,5 K más da ~43. Coherente.
- **Estaciones:** periastro el día 309 (1533,3 W/m²), apoastro el día 146 (1275,5). La temporada del periastro es 1,32 K más cálida que la del apoastro, en los dos hemisferios a la vez, como preveía el §3.
- **Banda de lluvias en el ecuador** (7,25 mm/día a 2,5° N; 7,28 a 2,5° S en la temporada del periastro y 7,70 a 2,5° N en la del apoastro); mínimos subtropicales a ±22,5° y máximos de latitudes medias a ±42,5°.
- **Coherencia de Donohoe et al. (2013):** la atmósfera cruza el ecuador hacia el sur (−0,22 PW) y la lluvia queda ligeramente al norte.
- **Transporte hacia los polos:** máximo total 2,81 PW a 40° N y 2,72 a 40° S (atmósfera 2,43 y 2,23; latente 1,47 y 1,39). La Tierra lleva ~5,5 PW; escalado por el área (0,734² = 0,54), ~3,0 PW: el mismo orden.
- **Células de Hadley:** +4,71 y −4,59·10¹⁰ kg/s, con bordes a +27,8° y −27,7°. La orientación de Held y Hou (modo Tierra × 1,082) daba +31,4° y −29,8°: 2–4° más estrechas, en el sentido de lo que hace más rotación (Kaspi y Showman 2015). No es un requisito (§4.7).
- **Hielo marino:** norte 4,23 / 3,28 y sur 3,72 / 2,07 millones de km² (máximo / mínimo), con poco ciclo anual (1,7°).
- **Tropopausa:** ~15,8 km en los trópicos y ~10,5–11 km en los polos; la tropical, demasiado caliente (−46,5 °C): radiación gris (limitación 5).

**Lo que no era lo esperado (§3 preveía tres células por hemisferio, como en la Tierra):**
1. **Dos células por hemisferio, sin célula polar:** la indirecta (tipo Ferrel) va de ~28° hasta el polo (−1,12 a 60° N; +0,92 a 55° S).
2. **Viento del oeste en superficie de ~31° al polo, con el máximo a 82,5° N (10,4 m/s) y 72,5° S (9,5 m/s).**
3. **Vórtice estratosférico de 91,2 y 85,0 m/s** a 62,5°, 3,7 hPa (modo Tierra 77 y 74; ERA-Interim 16,5 y 33,7).
4. **Viento del oeste en la alta troposfera ecuatorial** (+12,8 a +15,6 m/s a σ ~0,25).
5. **Tierra más cálida que el mar** (23,85 frente a 20,07 °C). Depende de dónde está la tierra en pizarra1; pendiente de mirarlo con el mapa.

**El control con el modo Tierra** (`herramientas/control_tierra_p3n.py`; detalle en `DISENO_FASE_6.3.md` §6.26) aplicó las mismas funciones a la climatología del modo Tierra de la 3.18.1:
- **1 es del código:** el modo Tierra tampoco forma célula polar.
- **4 es del código, más intenso en P3N:** el modo Tierra ya tiene +2 a +5 m/s (ERA-Interim, −3 a +0,3). La explicación por la célula de Hadley que cruza el ecuador queda descartada (en jun–ago del modo Tierra esa célula es fuerte y el sesgo es mayor).
- **2 es propio de P3N en su magnitud:** el máximo del oeste en superficie está a 52,5° N y 57,5° S en el modo Tierra (ERA-Interim 47,5° N y 52,5° S) y a 82,5° N y 72,5° S en P3N. Con el núcleo seco y los niveles de M3N (Held y Suarez con la geometría de la Tierra, 3.4.0) estaba a ~34°. Es la diferencia grande que el criterio §4.7 B obliga a explicar o investigar antes de la 4.0.0.
- **3:** limitación declarada (sin arrastre de ondas de gravedad, 6.6; radiación gris, 5.2), más intensa en P3N. Con 1,7° de inclinación no hay noche polar ni invierno: el vórtice es de todo el año en los dos hemisferios (lectura, no demostración).

### 8.4 Investigación del punto 2 (decisión de Carlos, 12/10/2026, 06:17: "A")

**Pregunta:** ¿qué parte del desplazamiento al polo de los vientos del oeste es geometría de P3N (radio, gravedad, rotación) y qué parte es la física (vórtice estratosférico, radiación)?

**Experimento (3.20.0):** Held y Suarez (núcleo seco, sin la física de M3N) con los niveles de M3N y la geometría de P3N, frente a la misma prueba con la geometría de la Tierra **con el mismo código**:
- `python held_suarez.py --niveles m3n --planeta p3n` y `python held_suarez.py --niveles m3n` (~25 min cada una en el PC; se pueden lanzar a la vez, cada una usa un hilo), después `python analizar_held_suarez.py`, que compara las latitudes de las dos.
- El forzamiento de HS94 es el mismo en las dos (mismas temperaturas de equilibrio, mismo contraste de 60 K entre el ecuador y el polo y mismos tiempos de relajación en segundos); solo cambian a, Ω y g. Sin ciclo diario, la duración del día no interviene.
- La corrida de la Tierra se repite con el código actual, para que las dos usen exactamente el mismo (la de la 3.4.0 se hizo con una versión anterior); la antigua se aparta a otra carpeta antes.
- Probado (entorno de la IA): 40 días de P3N estables (chorro de 42–45 m/s a 32,5° S, presión en superficie 969–1019 hPa); 3 pruebas automáticas (`test_held_suarez_planeta.py`).

**Cómo se leerá (fijado antes de ver el resultado) 🔶:** la latitud del máximo del oeste en superficie, P3N menos Tierra, con la misma configuración:
- **< 5°** (una celda de la rejilla, el límite de lo que se puede distinguir): la geometría de P3N no lo explica; la causa está en la física (sospechoso principal, el vórtice: Polvani y Kushner 2002) y queda ligada a la 6.6 y la 5.2 como limitación declarada.
- **≥ 5° pero mucho menos que los 25–30° que se separan P3N y la Tierra con el modelo completo:** parte geometría y parte física; se escriben las dos con su peso.
- **Comparable a esos 25–30°:** es sobre todo geometría del planeta, un resultado de P3N y no un fallo, más el sesgo hacia el polo que el código ya tiene en el modo Tierra (~+5°).

**Resultado (PC de Carlos, 12/10/2026, 07:07; 1200 días cada una, media de 1000; `analizar_held_suarez.py`):**

| Máximo del oeste en superficie | Tierra (código de la 3.20.0) | P3N | P3N − Tierra |
|---|---|---|---|
| Norte | 34,0° | 34,5° | +0,5° |
| Sur | 34,1° | 36,2° | +2,1° |
| Chorro (norte / sur) | 29,3° / 29,2° | 30,4° / 31,7° | +1,1° / +2,5° |

- **Lectura con el criterio fijado antes: < 5°.** La geometría de P3N (radio, gravedad y rotación) **no explica** el desplazamiento al polo de los vientos del oeste del modelo completo (de ~35° a 72,5–82,5°). La causa está en la física. Sospechoso principal: el vórtice estratosférico sin arrastre de ondas de gravedad (Polvani y Kushner 2002), más fuerte en P3N que en el modo Tierra (91 frente a 77 m/s) y de todo el año, sin estaciones. Queda como **limitación declarada**, ligada a la 6.6 y a la 5.2. No está demostrado que sea el vórtice; sí que no es la geometría.
- **Reproducibilidad:** la corrida de la Tierra con el código de hoy da exactamente lo mismo que la de la 3.4.0 (chorro 36,6 m/s a 29,3°; oeste de superficie 4,0 m/s a 34,0°).
- La diferencia entre el norte y el sur de P3N (30,4° y 31,7°), con un forzamiento simétrico, es ruido de muestreo: las dos mitades del promedio se separan ~1°.
- En seco, las capas altas (< 20 hPa) de P3N llegan a 16,6 m/s frente a 9,9 en la Tierra, y crecen un 14 % entre las dos mitades: también en seco, la estratosfera de P3N tiende a acelerarse más.
- Con esto, el punto 2 del §8.3 queda **investigado** como pide el §4.7 B.

### 8.5 La climatología oficial se rehace con las medias zonales por nivel (decisión de Carlos, 12/10/2026, 07:40: camino "a")

- **Por qué:** la vista Circulación de H3N (el corte latitud × altura con las células, los chorros y la temperatura; `AUDITORIA_H3N_EXPLORAR.md` de H3N, P2.1) necesita esas medias, y la climatología del 11–12/10 no las guardó.
- **Cómo (3.21.0):** `fase6_zonal.py` las acumula por día del año en los 30 años registrados, solo leyendo el núcleo (idéntico bit a bit, `test_fase6_zonal.py`). `clima_dinamico.py --continuar-desde outputs/clima_dinamico/p3n_3.19.0` hace una climatología **nueva** de 30 años desde el estado final de aquella, que está en equilibrio y sin deriva, sin repetir el arranque ni los 15 años de equilibrio. El origen queda en `origen.json`, en el resumen y en la exportación.
- **Coste:** ~5 h (30 años de P3N a ~37 pasos/s), después de las sensibilidades.
- **Consecuencia:** el clima oficial de la 4.0.0 será el de esta climatología nueva, otros 30 años del mismo estado de equilibrio. Las cifras de §8.2–§8.3 cambiarán dentro del ruido entre muestras de 30 años: el error de la media es ≤ 0,25 K por banda (§8.1). El análisis se repetirá sobre ella.

### 8.6 Parte D: sensibilidad al océano (PC de Carlos, 12/10/2026; `analizar_p3n.py --comparar`, 3.21.0)

Cuatro climatologías exploratorias de 12 años (3.20.0, ~6,1–6,5 h cada una, en paralelo; equilibrio en 13–15 años; todas cumplen el criterio de equilibrio), cada una menos la climatología de 30 años de la 3.19.0:

| Prueba | Global | Lluvia | Entre 50° S y 50° N (filas cada 10°) | Al polo de ~60° |
|---|---|---|---|---|
| Capa de mezcla × 0,5 | −0,03 K | −0,005 mm/día | de −0,16 a +0,05 K | de −0,24 a 0 K |
| Capa de mezcla × 2 | +0,01 K | +0,008 | de −0,06 a +0,05 K | de −0,39 a +0,29 K |
| D_oc × 0,5 | −0,50 K | +0,031 (+0,9 %) | de −0,32 a +0,43 K | hasta −8,06 K (82,5° S), −4,88 (72,5° S), −2,80 (77,5° y 67,5° N) |
| D_oc × 2 | +0,61 K | −0,060 (−1,7 %) | de −0,72 a +0,35 K | hasta +11,62 K (82,5° S), +9,68 (77,5° N), +7,11 (87,5° N) |

**Ruido.** El error de la media por banda es de 0,17–0,46 K en las de 12 años y ≤ 0,25 K en la de 30 (resúmenes, σ/√N_eff): diferencias por debajo de ~0,5 K en una fila no se distinguen de la variabilidad entre climatologías. Desde la 3.22.0, D da el error de cada diferencia.

**Lectura (medida):**
- **La capa de mezcla no cambia el clima medio**: todo dentro del ruido. Lo que sí debería cambiar es la **amplitud de las estaciones** (una capa más fina se calienta y se enfría más deprisa); las medias anuales no la miden. 🔶 Pendiente: medirla con la sección D de la 3.22.0 antes de escribir la incertidumbre definitiva.
- **La difusión del océano importa poco entre 50° S y 50° N (< 0,8 K) y mucho cerca de los polos (de varios K a más de 10 K).** El signo tiene sentido físico: más difusión lleva más calor del ecuador a los polos (el ecuador se enfría 0,7 K a 2,5° S; los polos se calientan).
- **Asimetría:** doblar la difusión calienta los polos más de lo que enfría reducirla a la mitad. 🔶 Hipótesis, no comprobada: la realimentación del hielo marino (menos hielo, mar más oscuro, más luz absorbida). Se comprueba con el hielo de la sección D de la 3.22.0.

**Consecuencia para la 4.0.0 (incertidumbre declarada):** el clima de P3N es robusto al océano entre los 50° de cada hemisferio (menos de 1 K); al polo de ~60°, la incertidumbre es de **varios grados, hasta más de 10 K**, hasta que M3N tenga corrientes (Fase 6.5). Coincide con lo sabido del modo Tierra (§4.2): una difusión no tiene la forma del transporte real del océano. La exportación y H3N lo dirán.

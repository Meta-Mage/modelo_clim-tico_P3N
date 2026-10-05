# DISEÑO — M3N v3.0: atmósfera de varias capas (en seco)

**Estado (04/10/2026, noche): código hecho, PENDIENTE DE CALIBRAR Y VALIDAR.** `fase30_multicapa.py` + interruptor I10 (apagado por defecto), modo Tierra y scripts de calibración (§13). Hay pruebas con un modelo de una columna (§9).
Por decisión de Carlos (04/10/2026), este diseño lo desarrolla la IA con los criterios de siempre:
- física verificable;
- nada terrestre sin justificar, con las etiquetas (a) física universal, (b) propiedad de P3N y (c) calibración terrestre;
- comparar con la Tierra;
- no esconder errores ajustando parámetros.

Luego se le explica el resultado. Lo que sea decisión de Carlos (worldbuilding, coste) se le pregunta aparte.

Antecedentes: `DISENO_FASE5A.md` (aviso inicial). La v3.0 se hace **antes** de la Fase 5a (que pasa a ser la v3.1), para construir el ciclo del agua una sola vez, sobre la estructura vertical definitiva.

---

## 0. Resumen

La atmósfera pasa de 2 capas (capa límite 1000–900 hPa + "troposfera" 900–0 hPa) a **N capas en coordenada sigma**, de momento sin vapor de agua.

Con N capas, el perfil vertical de temperatura, la tropopausa y la estratosfera **salen de la física** (radiación + convección) en lugar de reconstruirse con supuestos.

Además, la presión sobre las montañas pasa a ser la real **en toda la física**: es la opción C de la presión, que estaba pendiente.

---

## 1. Coordenada vertical: sigma (σ = p / p_s)

- Cada columna se divide en N capas, que son fracciones de **su propia** presión en superficie. Es la coordenada estándar de los modelos atmosféricos y sigue el relieve sin que las montañas "atraviesen" capas.
- **Presión en superficie de cada celda (b):** equilibrio hidrostático con la gravedad de P3N, p_s = p₀·exp(−g·z / (R·T̄)), con p₀ = 1 bar (decisión de Carlos del 01/10) y T̄ la temperatura media del aire entre el nivel del mar y la superficie.
  - **Consecuencia:** sobre una montaña hay menos aire encima. Eso significa menos capacidad calorífica, menos efecto invernadero y menos luz absorbida.
  - **Resuelve la opción C** que `DISENO_FASE5A.md` dejaba pendiente.
- **Reparto de los niveles:** más finos cerca del suelo (capa límite) y espaciado regular en altura aproximada por arriba (zona de la tropopausa). Receta de la prueba: z_k = 40 km·(k/N)^1,6, con σ = exp(−z/7,5 km). Es solo una receta de espaciado; no impone ningún dato físico.
- **Número de capas N:** se decide por **convergencia** (§9), no copiándolo de otro modelo. Provisional: 20–30.

## 2. Infrarrojo: transferencia "gris" de dos flujos

- **Ecuaciones (a):** flujos ascendente y descendente de infrarrojo, con la absorción de cada capa dada por su espesor óptico Δτ.
  - Fuente de Planck **lineal en τ** dentro de cada capa, la forma estándar de los esquemas de dos flujos. Con ella el resultado converge con menos capas.
  - Comprobado: en una atmósfera isoterma coincide exactamente con la fuente constante.
- **Perfil del espesor óptico** (forma de Frierson, Held y Zurita-Gotor 2006; ✅ verificada en el código de Isca, `two_stream_gray_rad.F90`):
  τ(p) = τ₀·[ f·(p/p₀) + (1 − f)·(p/p₀)⁴ ]
  - El término lineal representa los gases bien mezclados (CO₂…).
  - El término de grado 4 representa el vapor de agua, que se concentra abajo con una escala de altura ~¼ de la del aire.
  - **Presión absoluta** (p/p₀, no σ): sobre las montañas queda menos absorbente encima, como en la realidad.
- **Calibración (c), mismo método que en la Fase 2b:**
  - Con la Atmósfera Estándar US 1976 y el suelo a 289 K, τ₀ y f deben reproducir los valores **sin nubes** de Wild et al. (2019): DLR = 314 W/m² y OLR = 267 W/m².
  - Resultado, convergido con el número de capas: **τ₀ = 2,452, f = 0,240** (N = 20: 2,460 / 0,236; N = 160: 2,452 / 0,240).
- **Rechazado — dependencia con la latitud de Isca/Frierson** (τ = 6 en el ecuador y 1,5 en el polo): imita el vapor de la Tierra. P3N no tiene base para copiarla. El vapor del propio P3N entrará en la Fase 5b.
- **Provisional hasta la 5b:**
  - El exponente 4 es (c). Habrá prueba de sensibilidad con 3 y 5.
  - En la 5b, el término del vapor se sustituirá por el vapor calculado.

## 3. Luz del sol

- **Totales sin cambios:** transmisión media con τ_sw = 0,18 y 73 % de lo atenuado absorbido por la atmósfera (Wild et al. 2019, como en la v2.4.2). La parte que llega al suelo **no cambia**.
- **Reparto vertical de lo absorbido por la atmósfera:** con peso ∝ Δ[(p/p₀)⁴], la forma de absorción del vapor de Isca (`solar_exponent = 4` ✅). Antes iba entero a la "troposfera".
- **Sobre montañas:** la masa de aire que atraviesa la luz se escala con p_s/p₀ (física, a).

## 4. Convección

- **Ajuste convectivo seco** que conserva la entalpía (Manabe y Strickler 1964). Las capas que se enfrían con la altura más deprisa que un gradiente crítico se mezclan hasta quedar justo en él.
- **Gradiente crítico 6,5 K/km (c), PROVISIONAL.**
  - Es el mismo valor que ya usa la v2.4.2 a través de la Atmósfera Estándar, así que la v3.0 es comparable con ella.
  - En la **v3.1 (Fase 5a)** se sustituye por el **ajuste convectivo húmedo**, con el gradiente adiabático húmedo calculado para P3N (decisión de Carlos: "debemos calcularlo").

## 5. Superficie y capa límite

- **Calor sensible** (Louis 1982, implícito, como hoy) hacia la **capa más baja**.
- **De día** (aire inestable), el ajuste convectivo lo reparte hacia arriba: la altura de la capa límite **sale sola**.
- **De noche** (estable), se queda en la capa baja: las **inversiones nocturnas salen solas**.
- **Aire a 2 m:** perfil logarítmico entre la superficie y el centro de la capa más baja (a ~100–300 m, frente a ~410 m de la capa límite actual). Es más preciso que hoy.
- **Océano, hielo marino y columna de suelo: sin cambios.** Su interfaz pasa a ser la capa más baja.

## 6. Transporte horizontal (decisión 2 de la 5a, opción C′, en seco)

- **Flujo difusivo de la energía estática seca cerca de la superficie**, s = c_p·T + g·z, con **un solo coeficiente**.
- **Por qué incluir g·z (física, a):** una masa de aire que sube por una ladera sin intercambiar calor **conserva s**, así que s no genera un transporte falso entre llanuras y montañas. La temperatura sola sí lo generaría.
  - Esto resuelve el pendiente "si h incluye g·z" de `DISENO_FASE5A.md` §3.6.
  - En la v3.1, s pasa a ser la energía estática húmeda (se añade L·q).
- **Reparto vertical de lo que llega:** proporcional a la masa de las capas troposféricas. Implementado como las capas con σ ≥ 0,25 (`SIGMA_TOPE_TRANSPORTE`), una aproximación fija de la troposfera, en vez de la tropopausa diagnosticada en cada paso: así el operador es lineal y la energía se conserva exactamente (ΔT igual en esas capas = ΔX). **Provisional**, con prueba de sensibilidad frente a un reparto más concentrado abajo. Lo de fondo es que el transporte por remolinos es máximo hacia 850 hPa (Peixoto y Oort 1992; ⚠️ verificar).
- **Calibración (c):** D en modo Tierra (siempre 24 h) frente al transporte total de Trenberth y Caron (2001).
  - Cruce de control: el D obtenido, expresado como difusividad de remolinos, debe salir del orden de 1,06×10⁶ m²/s (Hwang y Frierson 2010).
  - Es previsible que el transporte siga "saturando" en seco; se resuelve en la v3.1 con el calor latente.
- **Rotación:** se implementa ya el escalado D ∝ 1/Ω² (Williams y Kasting 1997) con la velocidad de giro **sideral** de cada planeta, ya que D se recalibra en esta versión.
  - Hecho en la v2.4.3 (`parametros.FACTOR_ROTACION_D`), también para el modelo de dos capas. Con el día de 19,84 h, el factor es **0,684**. En el modo Tierra vale 1.
- **Océano (D_oc):** sin cambios de método; se recalibra.

## 7. Numérica y coste

- Radiación y ajuste convectivo **vectorizados** sobre las 2592 celdas.
- **Coste medido** en el entorno de la IA, que es ~4,5× más lento que el PC de Carlos: para N = 20, ~2,2 ms por paso la radiación y ~2,3 ms el ajuste, frente a ~2,7 ms por paso del modelo actual.
- **Estimación en el PC de Carlos:** de ~6 min a **~15 min por simulación completa**.
- **Margen:** llamar a la radiación cada 2 pasos (práctica habitual en los modelos climáticos), solo si una prueba demuestra que no cambia el resultado.
- La aceleración de la convergencia (salto geométrico y salto del hielo) se extiende a las N capas.
- La caché incluye los módulos nuevos en `MODULOS_FISICA`.

## 8. Interruptor y compatibilidad

- **I10 `atmosfera_multicapa`.** Apagado, el modelo debe reproducir la **v2.4.3 bit a bit** (prueba V0: ✅ superada el 04/10/2026, 1 año en el mapa de la pizarra, los 16 resultados idénticos).
- **I10 está APAGADO por defecto** hasta que pase la calibración y la validación (§10, §13): las exportaciones del clima de P3N siguen siendo las de la v2.4.3.
- **Exportación m3n-clima:** mismas claves. Se añaden como claves nuevas, sin subir la versión del formato, la presión en superficie por celda y, opcionalmente, la altura de la tropopausa.

## 9. Pruebas con una columna (04/10/2026, modo Tierra, media global)

Equilibrio radiativo-convectivo con la luz absorbida de Wild et al. sin nubes (214 W/m² en el suelo y 73 W/m² en la atmósfera), τ₀ y f calibrados y gradiente crítico de 6,5 K/km.

| N | Suelo (K) | Tropopausa | Estratosfera a 25 km (K) | Tiempo |
|---|---|---|---|---|
| 10 | 296,35 | 11,2 km, 229 K | 226,3 | 1,8 s |
| 15 | 297,09 | 10,6 km, 233 K | 225,4 | 4,3 s |
| 20 | 297,44 | ~8,6 km, 236 K | 225,1 | 10 s |
| 30 | 297,74 | ~8,8 km, 235 K | 225,1 | 27 s |
| 40 | 297,85 | ~9,0 km, 236 K | 225,1 | 61 s |

**Lectura:**
- **El suelo converge.** Las diferencias son de 0,3 K entre N = 20 y 30, y de 0,1 K entre 30 y 40. El N final se decidirá con el modelo completo, con el criterio "doblar N cambia la media global < 0,2 K y las medias por latitud < 0,3 K".
- **La estratosfera** sale a **225,1 K**, y la teoría de la atmósfera gris predice la temperatura de piel T_e/2^¼ = **224,3 K** (con T_e = 266,7 K para OLR = 287 W/m²). La física del esquema es correcta, y además **la tropopausa y la temperatura de piel salen solas**: la decisión 4 de la 5a (fórmula de la temperatura de piel) **deja de ser necesaria**.
- **El suelo sale a ~297 K** (la Tierra real, ~288 K) porque una Tierra **sin nubes** en equilibrio emite 287 W/m², no 267. Es lo mismo que ya se vio en la Fase 2b (§ "Calibración del infrarrojo"). No es un fallo.
- El diagnóstico de la tropopausa (criterio de 2 K/km) es ruidoso con capas gruesas. En el modelo completo se usará la definición de la OMM bien implementada.

## 10. Validación del modelo completo

1. **V0:** con I10 apagado, v2.4.3 bit a bit. ✅
2. **Energía:** desequilibrio < 10⁻³ (hoy ~10⁻⁵).
3. **Modo Tierra:**
   - Flujos (luz 247/73/53; DLR ≈ 314 y OLR ≈ 267 solo como comprobación de la calibración).
   - Transporte (Trenberth y Caron 2001).
   - Tropopausa: debe salir más alta en el ecuador que en los polos. Valores típicos de la Tierra: ~16–17 km en los trópicos y ~8–10 km en los polos (⚠️ verificar la fuente).
   - Estratosfera.
   - Hielo marino frente a la v2.4.2 y a las observaciones.
4. **P3N con prueba1 y con el mapa de la pizarra:** comparación por latitudes con la v2.4.2. Cada diferencia importante, explicada.
5. **Convergencia en N** (§9) y **paso de la radiación**.
6. **Sensibilidad:** exponente del vapor (3/4/5) y reparto vertical del transporte.

## 11. Efecto sobre las decisiones de la Fase 5a

- **Decisión 3 (humedad H2):** se generaliza a **humedad propia en cada una de las N capas**. Desaparece la integración sobre un perfil reconstruido, y la humedad relativa crítica de la condensación se puede acercar a la de los modelos de capas finas.
- **Decisión 4 (tropopausa por temperatura de piel):** **ya no hace falta**. La tropopausa sale de la radiación y la convección (§9).
- **Decisión 2 (C′):** el reparto vertical se resuelve aquí (§6).
- **Presión local (opción C):** resuelta (§1).

## 12. Limitaciones declaradas

1. **Radiación gris:** no distingue longitudes de onda. Revisión en la 5b / Fase 8.
2. **Sin ozono:** P3N **tiene oxígeno** (Carlos, 04/10/2026), así que tendrá capa de ozono, que absorbe ultravioleta y calienta la estratosfera. El esquema gris no lo incluye: la estratosfera del modelo sale más fría de lo que sería, y la tropopausa, peor marcada. Efecto en superficie pequeño. Se añade con la radiación por bandas (5b / Fase 8).
3. **Absorción del vapor con forma fija terrestre** (exponente 4) hasta la 5b.
4. **Gradiente crítico de 6,5 K/km impuesto** hasta la v3.1.
5. **Sin nubes.**
6. **Transporte horizontal difusivo** (remolinos) hasta la Fase 6.

## 13. Modo Tierra y herramientas de calibración (04/10/2026)

**Modo Tierra.** Con la variable de entorno `M3N_MODO=tierra`, `parametros.py` sustituye P3N por la Tierra:
- valores de la Earth Fact Sheet de la NASA ✅ (verificada el 05/10/2026);
- irradiancia 1361,0 W/m² a 1 UA (coincide con Kopp y Lean 2011, 1360,8 ± 0,5 ✅);
- año sideral 365,256 días, excentricidad 0,0167;
- el día 1 a las 0 h es el solsticio de diciembre (anomalía media 347,476°, calculada con la longitud del perihelio, 102,94719° heliocéntrica);
- oblicuidad 23,44°, día de 86 400 s, radio 6371 km ✅; g = 9,80665 m/s² (estándar, CODATA 2018) ✅;
- paso de tiempo 900 s.

`modo_tierra.py` da el mapa: tierra/agua con `global-land-mask` (20 × 20 puntos por celda, mayoría), **28,2 % de tierra** (real 29,2 %). Altitud: Antártida 2300 m, Groenlandia 2000 m, resto 0 m, igual que la calibración de la v2.2c.

**Limitación:** faltan las demás montañas. En la v3.0 la presión en superficie depende de la altitud, así que esto hace la Tierra simulada un poco más cálida sobre los continentes. Mejora pendiente: altitud media real por celda (ETOPO).

**Diagnóstico nuevo: transporte implícito.** Con I10 encendido, la simulación devuelve `flujos`, las medias anuales por celda de:
- la convergencia del transporte atmosférico (exacta, sale del propio paso implícito);
- la convergencia del transporte oceánico;
- el balance en lo alto de la atmósfera.

`fase30_multicapa.transporte_meridional()` las integra desde el polo norte y da los PW que cruzan cada latitud.

**`calibrar_v30.py`** (en el PC de Carlos, `M3N_MODO=tierra`). Es una rejilla D_atm ∈ {0,40; 0,55; 0,70; 0,90} × D_oc ∈ {0,12; 0,20}:
- sin hielo marino, igual que en la v2.2c (el hielo se valida después, con D fijado);
- simulaciones en paralelo;
- compara con Trenberth y Caron (2001, resumen ✅): transporte atmosférico máximo **5,0 ± 0,14 PW a 43° N**, parecido hacia 40° S; a 35° la atmósfera lleva el 78 % (N) y el 92 % (S) del total.

Con `--convergencia-n`, el mismo script repite la simulación con N = 10, 20, 30 y 40.

**Resultado de la calibración (05/10/2026, en el PC de Carlos).** Se hicieron 18 simulaciones en tres tandas, todas convergidas en 8–12 años:
- tanda 1: D_atm 0,40–0,90 × D_oc 0,12 / 0,20;
- tanda 2: D_atm 1,0–1,3 × D_oc 0,16 / 0,22;
- tanda 3: D_atm 1,35 / 1,5 × D_oc 0,22 / 0,28.

Selección:

| D aire | D océano | Pico atm. N / S (PW) | % océano a 35° N / S | Aire a 2 m: global / ecuador / polos N, S |
|---|---|---|---|---|
| 0,55 | 0,12 | 3,61 / 2,97 | 11,8 / 19,6 | 19,9 / 31,0 / −5,5, −24,1 |
| 0,90 | 0,20 | 4,55 / 3,63 | 12,3 / 21,2 | 20,0 / 28,7 / 0,3, −18,4 |
| 1,30 | 0,22 | 5,32 / 4,22 | 10,3 / 18,3 | 20,0 / 27,3 / 4,0, −15,2 |
| 1,50 | 0,28 | 5,53 / 4,30 | 11,4 / 20,4 | 20,0 / 26,6 / 5,6, −13,7 |
| Observado | | 5,0 a 43° N / parecido a 40° S | 22 / 8 | |

**Lectura:**
- **El transporte ya no satura** como en la v2.4: sigue creciendo con D, cada vez algo menos.
- **La temperatura global no depende de D:** 19,9–20,0 °C en las 18 simulaciones. El exceso de calor no viene del transporte (§13, más abajo).
- **El sur transporta siempre ~20 % menos que el norte.** En la realidad son parecidos. Causas probables, por cuantificar con la validación:
  - el océano del sur se lleva más parte;
  - la meseta antártica, por el término g·z de la energía estática;
  - el sur es casi todo mar, con gradientes más suaves.
  En la realidad pesan además las borrascas más intensas del océano austral y el calor latente, que una difusión con un solo D no puede representar.
- **El reparto entre océano y atmósfera tiene la asimetría al revés** que la observada: el modelo da más océano en el sur. El océano solo difunde; sin corrientes reales (Fase 6b) no se puede corregir.

**Criterio** (acordado con Carlos): no favorecer a ningún hemisferio.
- La media de los dos hemisferios del pico atmosférico = 5,0 PW.
- La media de la parte del océano a 35° = 15 %.

**Ajuste:** una superficie cuadrática sobre las 14 simulaciones con D_atm ≥ 1,0 (residuo máximo 0,003 PW y 0,02 %), resuelta para los dos objetivos. Da **D_atm = 1,55** y **D_oc = 0,27**.
- Queda 0,05 por encima del D_atm más alto simulado, una extrapolación mínima.
- **La simulación de validación con estos valores sirve de comprobación directa.**
- Valores previstos: pico 5,62 N / 4,38 S; océano 10,7 % N / 19,3 % S.

**Cruce con Hwang y Frierson (2010):**
- Expresado como difusividad de remolinos sobre la masa que recibe el transporte (σ ≥ 0,25), D_atm = 1,55 equivale a K ≈ 8×10⁶ m²/s.
- En forma de D sobre la energía de toda la columna, 1,55 frente a 0,27 W/m²/K: unas **6 veces el de Hwang y Frierson** (1,06×10⁶ m²/s, aplicado a la energía estática HÚMEDA).
- Explicación parcial:
  - el gradiente meridional de la energía húmeda es aproximadamente el doble que el de la seca;
  - sin nubes, la diferencia de temperatura entre el ecuador y los polos del modelo es unas 2 veces menor que la real: ~31 K frente a ~60 K. Para el valor real se toman el ecuador a ~27 °C y los polos a −18 °C (N) y −49 °C (S), como en `DISENO_FASE3.md` §6.1.
- Juntas explican un factor ~4. **El resto (~1,5) queda sin explicar**: se revisa en la v3.1, cuando el vapor entre en el transporte y D baje.

**Para P3N:** 1,55 × 0,684 = **1,06**. D_oc = 0,27 no se escala.

**`validar_v30.py`** (05/10/2026). Va después de la calibración, con el D elegido y el hielo encendido.
- **En modo Tierra** compara con la realidad:
  - balance global, frente a Wild et al. 2019;
  - aire a 2 m por bandas;
  - perfil vertical por bandas: tropopausa con la definición de la OMM, buscada entre 550 y 75 hPa como Reichler, Dameris y Sausen (2003) ✅, temperatura a 25 km y gradiente medio entre 0 y 6 km;
  - hielo marino, frente a NSIDC;
  - transporte.
- **En P3N** compara la v3.0 con la v2.4.3 en el mapa activo.
- **Pendiente:** una climatología observada del aire a 2 m por latitudes, con fuente verificada, para la comparación por bandas en modo Tierra.

Para eso la simulación devuelve también, en `flujos`, las medias anuales de:
- DLR y OLR por celda;
- luz absorbida en el suelo y en el aire;
- perfil de temperatura y presiones.

**Primera prueba (1 año sin converger, D_atm = 0,55, D_oc = 0,12):**
- pico atmosférico 3,45 PW (40° N) y 2,88 PW (40° S);
- océano a 35°: 12 % (N) y 20 % (S);
- aire a 2 m: global 20,0 °C, ecuador 30,7 °C.

El esquema funciona. Los valores definitivos saldrán de la rejilla convergida.

**Coste medido** (entorno de la IA, ~4,5× más lento que el PC de Carlos): 262 s por año terrestre (35 065 pasos), unos 7,5 ms por paso. Se reparte así:
- infrarrojo ~1,3 ms;
- ajuste ~0,8 ms;
- resoluciones implícitas ~0,5 ms;
- el resto, operaciones de numpy dentro del paso.

Llamar a la radiación cada 2 pasos ahorraría solo un ~15 %, así que **no se hace** (no compensa el riesgo).

**La v3.0 sale más cálida que la v2.4.3 (pendiente de entender del todo; NO se ajusta nada para esconderlo):**
- **P3N**, mapa de la pizarra, 1 año sin converger: 19,8 °C frente a 14,4 °C.
- **Modo Tierra:** 20,0 °C frente a 14,5 °C de la v2.2c.

**Explicación más probable (a):**
- Las dos versiones simulan una Tierra **sin nubes**, que absorbe ~287 W/m² y en equilibrio tiene que emitir lo mismo, 20 W/m² más que los 267 observados con cielo despejado.
- En la v3.0, las capas que emiten al espacio están atadas al suelo por el gradiente vertical (6,5 K/km). Para emitir 20 W/m² más se calienta toda la columna: unos 6–8 K, lo esperado para la respuesta de Planck de unos 3 W/m²/K.
- En la v2.4.3, la "troposfera" es un solo bloque que puede calentarse por su cuenta, con un umbral de ajuste de 40,5 K. Emite el exceso calentando ese bloque, sin calentar el suelo.
- La v3.0 es la físicamente más realista. La "Tierra sin nubes a 14,5 °C" de la v2.2c era, en parte, un artefacto de las dos capas.

**Consecuencia:**
- La masa de S3N (0,874) se calibró en la v2.2c para un aire a 2 m de 15,94 °C con el modelo de dos capas.
- Con la v3.0, el mismo Sol da un P3N unos 5 °C más cálido mientras no haya nubes.
- **Decisión de Carlos, cuando la v3.0 esté validada:** recalibrar la estrella ya, o esperar a las nubes (Fase 5b), que enfriarán.

## 13 bis. Validación de la v3.0 (05/10/2026, en el entorno de la IA, con hielo)

**Modo Tierra** (D_atm = 1,55, D_oc = 0,27; 10 años hasta el equilibrio; desequilibrio 1,6×10⁻⁵):

| | Modelo | Tierra sin nubes (Wild et al. 2019) |
|---|---|---|
| Luz absorbida por el suelo | 225,2 | 214 |
| Luz absorbida por el aire | 71,4 | 73 |
| Luz reflejada | 43,7 | 53 |
| DLR | 334,0 | 314 (atmósfera real, con nubes) |
| OLR | 296,6 | 267 (atmósfera real, con nubes) |

**Lo que sale bien:**
- **El transporte confirma la calibración:** 5,62 PW en el norte y 4,37 PW en el sur, frente a 5,62 y 4,38 previstos. La pequeña extrapolación era correcta.
- **Gradiente medio entre 0 y 6 km:** 6,4–6,5 K/km entre 60° N y 60° S, el impuesto.
- **Estratosfera a 25 km:** −42 a −51 °C.

**Lo que sale MAL** (no se ajusta nada para esconderlo):
1. **Polos demasiado cálidos y ningún hielo marino.**
   - Aire a 2 m entre 90° y 60° N: +7 °C (real, unos −15 °C).
   - En el sur, el océano a +13 °C.
   - Hielo marino: 0, frente a 15,6/18,5 millones de km² de máximo.
   - El aire a 2 m medio es 20,0 °C.
2. **Tropopausa demasiado baja y casi plana:** unos 10,3 km en todas las latitudes (real, ~16–17 km en el trópico y ~8–9 en los polos).
3. **OLR casi plano:** 271 W/m² en los polos y 314 en el ecuador (real, ~180 y ~290).
   - Los polos reciben por transporte lo que en la Tierra (balance en lo alto, −130 W/m², ≈ real).
   - Pero tienen que emitir mucho más, porque absorben más luz al no haber hielo.
   - Y porque **el espesor óptico gris es el mismo en todas partes**: en la Tierra, el aire polar, frío y seco, es mucho más transparente al infrarrojo que el tropical.

**P3N** (mapa de la pizarra, D = 1,55 × 0,683 = 1,06):

| | v3.0 | v2.4.3 (con la corrección de la capa de mezcla) |
|---|---|---|
| Aire a 2 m global | 23,0 °C | 14,5 °C |
| 90°–60° N | +5,0 °C | −24,6 °C |
| 10° N – 10° S | 32,0 °C | 31,1 °C |
| Hielo marino | 0 | 6 % N / 11 % S |

El calentamiento es casi todo polar: +30 K en los polos y +1–3 K en el trópico.

**Diagnóstico:**
- La causa principal es estructural: **el vapor de agua no es radiativamente activo**. El espesor óptico τ es fijo y está calibrado con la media global, que dominan los trópicos húmedos. Eso da a los polos un efecto invernadero demasiado fuerte y a los trópicos uno algo débil.
- Lo empeoran tres cosas:
  - la falta de nubes;
  - un D grande, necesario para transportar los 5 PW solo con calor seco;
  - la retroalimentación del hielo: sin frío no se forma hielo, y sin hielo los polos absorben más luz.
- En el modelo de dos capas, los polos eran fríos por un motivo equivocado: el bloque alto podía calentarse sin calentar el suelo.

**Matiz (análisis con columnas, `DISENO_V3.1.md` §8 bis):**
- Lo que más cambia con un τ que depende del vapor no es lo que sale por arriba, sino **el infrarrojo que llega al suelo polar**: en el invierno subártico, 220 W/m² con τ fijo frente a 124 con vapor.
- Con τ fijo, la superficie polar no puede enfriarse lo bastante para formar hielo, y sin hielo absorbe mucha más luz.

**Consecuencia:** la v3.0 tal cual **no debe encenderse para P3N**.

**La solución física:** que el espesor óptico dependa del vapor que calcula el propio modelo (vapor radiativamente activo, previsto para la 5b). La v3.1 ya tiene ese vapor. **Decisión de Carlos:** adelantar esa parte de la 5b.

Se descartó usar el τ dependiente de la latitud de Frierson/Isca porque es una forma terrestre (c) impuesta.

## 14. Revisión del código y correcciones (05/10/2026)

Una revisión independiente del código, hecha por un agente aparte con experimentos propios, **no encontró errores de energía, signos ni índices en la física**:
- el cierre de energía paso a paso queda por debajo de 10⁻⁷ W/m² con 2 y con N capas, con y sin hielo;
- el infrarrojo coincide con la solución analítica;
- el ajuste convectivo con numba coincide con la forma cerrada hasta 2×10⁻¹¹ K.

Corregido:
1. **La capa de mezcla del océano dependía de la duración del día** (`fase2_inercia_multicapa.py`).
   - Su capacidad se calculaba como inercia·√(día/π), una fórmula heredada de la capa diurna de la Fase 0.
   - Con 24 h daba los ~48 m documentados; desde la v2.4.3 (19,84 h) daba **43,3 m**, un 9 % menos de capacidad que nadie había decidido.
   - Ahora es fija: la capacidad con 24 h, ~48,7 m de agua de mar.
   - **Cambia algo el clima de P3N exportado con la v2.4.3:** un poco más de inercia en el océano, con estaciones algo más suaves en el mar. En modo Tierra no cambia nada, así que la calibración de D sigue valiendo.
2. `radiacion.py` (y `atmosfera.py`) no estaban en la huella de la caché. Ahora la caché se invalida también si cambian.
3. **La caché escribe de forma atómica** (archivo temporal + renombrar): dos procesos ya no pueden dejar un archivo a medias.
4. **El informe y la exportación** muestran los D que usa de verdad la simulación (`D_usados`, ya escalados con la rotación), no los de la Tierra.
5. **Convergencia:** `simular_fase2b` devuelve un indicador explícito, `convergido`. Antes, el aviso de "no convergió" saltaba también si convergía justo en el año 50.
6. **Modo Tierra:** la caja de Groenlandia ya no pone 2000 m en Baffin ni en Ellesmere.
7. **`validar_v30.py`:**
   - el gradiente de 0 a 6 km leía de menos: daba 6,29 en un perfil exacto de 6,5 K/km;
   - el filtro de celdas bajas estaba descrito de tres formas distintas.
8. Código sin uso eliminado y comentarios desfasados (año de 270 días) corregidos.

**Optimización:** infrarrojo con funciones compiladas (numba), unas 2 veces más rápido. Mismo resultado salvo redondeo: diferencia máxima 2×10⁻¹³ W/m².

## Referencias

✅ = verificado en esta sesión.
- Frierson, D. M. W., Held, I. M., y Zurita-Gotor, P. (2006). A gray-radiation aquaplanet moist GCM. Part I. *J. Atmos. Sci.*, 63. Forma τ(p) verificada en el código de Isca ✅.
- Isca (ExeClim): `two_stream_gray_rad.F90` (τ_eq = 6, τ_pole = 1,5, `linear_tau` = 0,1, exponentes 4) ✅.
- Hwang, Y.-T., y Frierson, D. M. W. (2010). *GRL* ✅.
- Manabe, S., y Strickler, R. F. (1964). Thermal equilibrium of the atmosphere with a convective adjustment. *J. Atmos. Sci.*, 21. ⚠️ cita por verificar.
- Pierrehumbert, R. T. (2010). *Principles of Planetary Climate*, cap. 4 (temperatura de piel). ⚠️ cita por verificar; la relación se ha comprobado numéricamente (§9).
- Peixoto, J. P., y Oort, A. H. (1992). *Physics of Climate*. ⚠️ por verificar.
- Trenberth, K. E., y Caron, J. M. (2001); Wild, M. et al. (2019); Williams, D. M., y Kasting, J. F. (1997); Louis, J.-F. (1982): ya usadas en fases anteriores.
- U.S. Standard Atmosphere (1976).

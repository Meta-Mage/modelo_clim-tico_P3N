# DISEÑO — FASE 5 de M3N (borrador de la 5.1, antes «5a»): ciclo del agua (v2.5 prevista; se hizo como v3.1)

> **Nomenclatura (08/10/2026, M3N 3.12.1):** este documento se llamaba `DISENO_FASE5A.md`. Fases y versiones van con la numeración normalizada; la correspondencia con los nombres anteriores está en `NOMENCLATURA.md` y `CHANGELOG.md`.

**Estado: BORRADOR v0 (04/10/2026), para revisión de Carlos. No hay código escrito.**

> **EN PAUSA — decisión de Carlos (04/10/2026, opción 2): antes de la 5.1 se construye la ATMÓSFERA DE VARIAS CAPAS** (adelanto de la antigua Fase 8, hoy Fase 2.3), en seco, como fase propia, y la 5.1 se hará después **encima** de ella, para no construir el ciclo del agua dos veces. Estado de las decisiones de este documento tras ese cambio:
> - Decisión 1 (calor latente constante, opción A; C pendiente): **sigue vigente**.
> - Decisión 2 (transporte C′, energía estática húmeda cerca de la superficie): **sigue vigente**; el reparto vertical del calor que llega se resolverá con las capas.
> - Decisión 3 (humedad H2) y decisión 4 (tropopausa por temperatura de piel, opción B, aceptada el 04/10): **a revisar** en el marco de varias capas; previsiblemente se simplifican (el perfil y la tropopausa saldrían de la física).
> - Decisiones 5 (albedo espectral de la nieve) y 6 (v2.4.2, ya hecha): sin cambios.
> - Numeración de versiones: **pendiente de confirmar por Carlos** (propuesta: la atmósfera de varias capas es la reescritura mayor del núcleo prevista como v3.0, así que pasaría a ser v3.0 y la 5.1 la v3.1).
Base: M3N v2.4.1. Este documento sigue el método de las fases anteriores: diseño → decisiones de Carlos → implementación con interruptores → validación → exportación.

Cada parámetro lleva una **etiqueta de origen** y un **estado de verificación**:

| Etiqueta | Significado | Regla |
|---|---|---|
| **(a)** | Física universal (propiedades del agua, del aire, leyes físicas) | Vale para P3N tal cual |
| **(b)** | Propiedad de P3N (gravedad, radio, rotación, órbita, estrella, presión) | Se calcula para P3N; nunca se copia de la Tierra |
| **(c)** | Calibración o valor empírico obtenido en la Tierra | Solo se usa si depende de cosas que P3N comparte con la Tierra; se justifica caso a caso |

| Estado | Significado |
|---|---|
| ✅ verificado | Comprobado en la fuente primaria o en el código fuente de un modelo publicado (se indica cuál) |
| ⚠️ pendiente | Valor orientativo; NO se usará en código hasta verificarlo |

---

## 0. Resumen

La Fase 5.1 añade el **ciclo del agua sin nubes radiativas**: evaporación, vapor de agua en el aire y su transporte, condensación y precipitación (lluvia y nieve), agua del suelo, nieve en tierra y **convección húmeda** (el gradiente vertical de temperatura pasa a calcularse para P3N en lugar de suponerse 6,5 °C/km).

Lo que la 5.1 **no** hace (y por qué):
- **Nubes y efecto radiativo del vapor** → Fase 5.2. Las emisividades actuales ya contienen el vapor terrestre; tocarlas ahora lo contaría dos veces.
- **Célula de Hadley y vientos** → Fase 6. Sin ella, la lluvia **tropical** de la 5.1 NO es válida (ver §3.6).
- **Lluvia orográfica direccional** (barlovento/sotavento) → Fase 6.
- **Presión variable con la altitud en toda la física** (opción C) → revisión aparte, PENDIENTE.

---

## 1. Decisiones de Carlos (04/10/2026)

1. La Fase 5 se parte en **5.1** (antes «5a», ciclo del agua) y **5.2** (antes «5b», nubes y radiación).
2. **No se fija** si P3N es más húmedo o más seco que la Tierra: sale lo que dé la física; se revisará después.
3. La **lluvia orográfica direccional espera a la Fase 6**.
4. La 5.1 usa **difusión pura** para el vapor, con el código preparado para separar "remolinos" y "circulación media" en la Fase 6. La lluvia tropical de la 5.1 se declara no válida.
5. **Presión local solo para la saturación** (opción B). La opción C (presión local en toda la física) queda **PENDIENTE de forma explícita**.
6. El gradiente vertical de temperatura **DEBE calcularse para P3N**, no suponerse → **ajuste convectivo húmedo**.
7. **Viento**: en modo Tierra, el viento medio real de la Tierra; en P3N, 5 m/s provisional hasta la Fase 6 (opción B).
8. **Sublimación** sobre el hielo marino con su calor latente; la nieve que cae sobre el hielo se suma a su masa; el **aislamiento por la nieve** queda como mejora pendiente.
9. Riesgos 1–7 aceptados con las soluciones propuestas (§7).
10. Principio general: **no asumir valores ni ecuaciones pensados para la Tierra** sin justificarlo (etiquetas a/b/c).

Pendiente de Carlos (con Ozan): **periodo de rotación definitivo** (§9).

---

## 2. Prerrequisito: v2.4.2 — quitar el "día de 24 h" escondido

**Estado: HECHO en la v2.4.2 (04/10/2026)**: ángulo horario general, comprobación de que el día es múltiplo del paso, registro horario y exportación con el número real de horas. Verificado: datos de la órbita idénticos bit a bit con 24 h; 21 pruebas superadas; simulación de prueba de un año con un día de 18 h (360 días, 18 registros horarios por día).

Hallazgo de la revisión del 04/10: el código supone días de 24 h en sitios donde no debería.

| Sitio | Problema |
|---|---|
| `geometria.py`, `angulo_horario()` | `(hora - 12) * PI / 12`: el sol da una vuelta cada 24 h aunque `ROTACION_PERIODO` cambie. Con un día de 30 h daría 1,25 vueltas por día, sin error. |
| `temperatura.py` y llamadas | `hora = (t % ROTACION_PERIODO) / 3600`: la "hora" es de 3600 s; hay que decidir qué es una hora en P3N |
| `fase2b_atmosfera.py`, `nieve_permanente_posible` | "meses" de 270/12 = 22,5 días: depende del número de días del año |
| Registro horario y exportación | Supone 24 registros por día |
| H3N | Textos y ejes con "270 días" |

**Propuesta:** una versión **v2.4.2** antes de la 5.1 que exprese todo en función de `ROTACION_PERIODO` (fracción del día en lugar de "hora"). **Prueba de aceptación:** con `ROTACION_PERIODO = 86400` el resultado debe ser **idéntico bit a bit** al de la v2.4.1.

**Sutileza a decidir:** el código trata `ROTACION_PERIODO` como el **día solar** (de mediodía a mediodía). El periodo de rotación sideral es algo más corto: 1/P_sol = 1/P_sid − 1/P_orb (giro en el mismo sentido que la órbita). Con 270 días por año la diferencia es ~0,37 %. Hay que dejar escrito cuál de los dos se define.

---

## 3. Física

### 3.1 Saturación del vapor (a) ✅

Presión de vapor de saturación sobre **agua líquida** y sobre **hielo**: solución de **Ambaum (2020)**, ecs. 13 y 17, que integra la ecuación de Clausius-Clapeyron con un calor latente que varía con la temperatura (ley de Kirchhoff):

- e_s,agua(T) = e₀ · (T₀/T)^((c_pl − c_pv)/R_v) · exp[(L_v0/T₀ − L_v(T)/T) / R_v], con L_v(T) = L_v0 − (c_pl − c_pv)(T − T₀)
- e_s,hielo(T): igual con L_s0, c_pi y L_s(T) = L_s0 − (c_pi − c_pv)(T − T₀)

Verificado contra la implementación de **MetPy** (`metpy.calc.thermo`, funciones `_saturation_vapor_pressure_liquid/solid`). Comprobación numérica propia (04/10) frente a la de Alduchov y Eskridge (1996, coeficientes 610,94 Pa / 17,625 / 243,04 °C, ✅): diferencias < 0,3 % entre −40 y +40 °C.

**Por qué no la Clausius-Clapeyron con L constante:** sobrestima e_s un **2,7 % a 30 °C** y un **4,7 % a 40 °C** (cálculo propio, 04/10). Justo en el océano tropical, donde más se evapora: daría ~5 % de evaporación de más.

Humedad específica de saturación: q_s = ε·e_s / (p − (1 − ε)·e_s), con ε = M_agua/M_aire = 0,62198 (a).

Sobre agua y sobre tierra sin nieve se usa e_s,agua; sobre hielo marino y nieve, e_s,hielo.

### 3.2 Hallazgo: coherencia energética del calor latente — DECIDIDO (opción A)

Si el calor latente varía con la temperatura (como en 3.1), el agua que se evapora a 30 °C y se condensa a −10 °C se lleva una cantidad de energía y devuelve otra distinta. La diferencia, (c_pl − c_pv)·ΔT ≈ 2360 × 40 ≈ 9×10⁴ J/kg (~4 % de L), en la realidad la compensa el calor sensible del propio vapor y de las gotas. M3N no lleva la cuenta de eso: su aire es "seco" a efectos de capacidad calorífica.

Con ~80 W/m² de calor latente global, una contabilidad inconsistente fugaría del orden de **2–3 W/m²**, frente a un balance actual de 10⁻⁵. Inaceptable.

**Decisión de Carlos (04/10/2026): opción A.** Calor latente **constante** para el **balance de energía** (L_v0 = 2,50084×10⁶, L_f = 3,337×10⁵, L_s0 = L_v0 + L_f J/kg; valores de MetPy, §5) y fórmula de Ambaum **solo para la saturación**. La energía cierra exactamente.

Práctica documentada ✅ (verificado en el código fuente, 04/10/2026):
- Isca (`src/shared/constants/constants.F90`): `HLV = 2.500e6`, `HLF = 3.34e5`, `HLS = HLV + HLF`, constantes; los usan su flujo de superficie y su convección.
- CESM (`CESM_share/src/shr_const_mod.F90`): `SHR_CONST_LATVAP = 2.501e6`, `SHR_CONST_LATICE = 3.337e5`, `LATSUB = LATICE + LATVAP`, constantes.

Error que se acepta, cuantificado (cálculo propio con la ley de Kirchhoff, L(T) = L_v0 − (c_pl − c_pv)(T − T₀)): cada kilo evaporado a 30 °C se lleva un **2,8 %** más de energía de la real; cada kilo condensado a −10 °C suelta un **0,9 %** menos. Es menor que otras incertidumbres de la 5.1 (viento fijo, RH_ref).

**PENDIENTE PARA EL FUTURO (opción C, decisión de Carlos):** sustituir por una contabilidad completa de la entalpía del agua (calor latente variable con la temperatura **y** calor sensible del vapor y de la precipitación), que es la solución exacta. Se abordará cuando el resto del modelo sea lo bastante preciso para que este 1–3 % importe. Opción B (calor latente variable sin más) descartada: crea energía (~2–3 W/m²).

### 3.3 Evaporación, sublimación y rocío (a + c)

Fórmula bulk, la misma estructura que el calor sensible actual:

E = ρ · C_E · U · β · (q_s(T_sup, p_local) − q_aire)

- **C_E = C_H** (mismo coeficiente y misma corrección de estabilidad de Louis que el calor sensible), con z₀q = z₀h. Es la hipótesis de "longitudes de rugosidad iguales para temperatura y humedad" de los modelos idealizados de Frierson et al. (✅ citado en Frierson 2007b; ✅ también en el código de Isca). El valor neutro sobre el océano del código actual, C_HN = 1,1×10⁻³ (Large y Pond 1982; Smith 1988), es (c): depende de la física del intercambio aire-agua, no del planeta.
- **β** (disponibilidad de agua): 1 en el océano; en tierra, el "cubo" (§3.5).
- **E < 0 permitido**: rocío o escarcha.
- **Sobre el hielo marino**: sublimación con L_s y e_s,hielo; se resta del espesor del hielo y entra en el balance de su superficie (que ya se resuelve por Newton en `paso()`).
- **Resolución implícita** (Riesgo 4): q_s(T^{n+1}) ≈ q_s(T^n) + (dq_s/dT)·(T^{n+1} − T^n), resuelto a la vez que el calor sensible. Es un sistema lineal pequeño por celda.
- **Límite**: la evaporación de un paso nunca vacía el cubo por debajo de 0.

### 3.4 Presión local para la saturación (b) — opción B

p_local = p_s · exp(−g·z / (R_d·T̄)), con g = 9,053 m/s² de P3N y T̄ la temperatura media de la capa de aire inferior en esa celda. **Solo** entra en q_s. El resto del modelo sigue con 1 bar. **Limitación declarada: la opción C queda PENDIENTE** (§10).

A 5000 m, con T̄ ≈ 255 K: p ≈ 540 hPa. Sin esta corrección, la capacidad de vapor en las cumbres sería casi la mitad de la real.

### 3.5 Agua del suelo: el "cubo" de Manabe (1969) (c)

- Capacidad W_fc = **150 mm** (Manabe 1969). ✅ Valor por defecto en el código de Isca (`max_bucket_depth_land = 0.15`, "default from Manabe 1969").
- β = min(1, W / (0,75·W_fc)). ✅ Así está implementado en Isca (`surface_flux.F90`). ⚠️ No lo he verificado en el artículo original de Manabe; la fuente secundaria del ECMWF (2002) lo describe solo como "proporcional al contenido".
- Lo que pasa de W_fc es escorrentía y **vuelve al océano al instante** (sin ríos; simplificación declarada).
- Es (c): describe suelos y vegetación terrestres. Para P3N es **provisional** hasta los biomas (Fase 7). Limitación conocida del cubo: sobrestima la evaporación del suelo desnudo (ECMWF 2002 ✅).

### 3.6 Transporte de calor y vapor — DECIDIDO (opción C′, Carlos 04/10/2026)

**Historia de la decisión.** El borrador v0 proponía atar el vapor a D_atm con una sola difusividad K. Al revisarlo (04/10) se vio que en M3N eso es incorrecto: el calor se difunde en la TR (900–0 hPa, 90 % de la masa, gradiente meridional pequeño) y el vapor vive abajo; el D_atm = 2,4 equivale a K ≈ 1,06×10⁷ m²/s en modo Tierra (D·a²/C_TR con g terrestre), unas 10 veces el K de Hwang y Frierson (2010, 1,06×10⁶ m²/s ✅), y ligar el vapor a ese K exageraría su transporte. Se consideraron: (A) K única ligada a la TR — descartada; (B) dos coeficientes calibrados con dos observaciones (transporte total y latente); (C) difundir calor y vapor juntos en la capa límite — descartada tal cual (metería todo el transporte en el primer km); (C′) — **elegida**.

**Opción C′ (marco de Hwang y Frierson 2010 ✅ y Siler, Roe y Armour 2018 ✅):**
1. El flujo meridional y zonal de energía de la atmósfera es difusivo en la **energía estática húmeda cerca de la superficie**, h = c_p·T + L_v0·q (+ g·z, ver abajo), con **un solo coeficiente**. En el código de Bonan et al. (✅) la difusión actúa sobre h/c_p ("temperatura equivalente") con D·p_s/(a²·g)·c_p.
2. La **convergencia de energía** que llega a cada celda se reparte: la parte de vapor va a la humedad (§3.7) y la parte de calor sensible se reparte entre la CL y la TR. **Regla de reparto: PENDIENTE** de justificar con fuentes (candidatas: proporcional a la masa de cada capa, o según el perfil vertical observado del transporte por remolinos).
3. **Un solo parámetro libre** (D_atm), recalibrado en modo Tierra frente al transporte TOTAL de Trenberth y Caron (2001). El reparto latente/sensible que salga es **validación independiente**. Comprobación cruzada: el D resultante, expresado como K para la Tierra, debe ser del orden del de Hwang y Frierson (1,06×10⁶ m²/s); si sale muy distinto, hay que entender por qué.
4. El océano (D_oc) sigue con su difusión propia en la capa de mezcla, recalibrada.
5. Traslado a P3N: D en unidades de M3N (W/m²/K, laplaciano sobre la esfera unidad) no depende del radio (Williams y Kasting 1997), así que K ∝ a²; y D ∝ 1/Ω² con la rotación de P3N.

**Por qué C′ y no B:** marco publicado y validado; un parámetro en vez de dos; calor y vapor coherentes por construcción; y la parametrización de la célula de Hadley de la Fase 6 (Siler et al. 2018) está **construida sobre este mismo marco**.

**Coste aceptado:** cambia dónde actúa la difusión de la Fase 2.2 (interruptor I6, `difusion_reubicada`) → hay que **revalidar la 2.2 y el hielo** (Fase 3) en modo Tierra y en P3N. La resolución numérica debe ser implícita y estable, con la humedad dentro.

**Física que lo apoya (⚠️ verificar en la fuente antes del código):** en latitudes medias, el transporte de calor por remolinos es máximo en la troposfera baja-media (~850 hPa), con un segundo máximo cerca de la tropopausa; el de vapor, casi todo por debajo de ~700 hPa (Peixoto y Oort 1992). El gradiente cerca de la superficie es el que "ve" ese transporte.

**Pendiente de decidir:** si h incluye g·z (energía potencial). Sobre montañas, la superficie está más alta; Hwang y Frierson usan h en superficie sin topografía (mundo acuático o Tierra suavizada). A decidir con la opción de presión (§3.4) y con las fuentes.

**Fase 6 (decidido):** el código separa ya el término de remolinos, para poder multiplicarlo por un peso w(φ) cuando entre la célula de Hadley. Forma de Siler et al. (2018) ✅ (código de Bonan et al.): w = 1 − exp(−(sen φ / σ)²), con σ = 0,3 para la Tierra. **σ para P3N deberá deducirse en la Fase 6**: depende de Ω·a (§9).

> **SUSTITUIDO (05/10/2026, decisión de Carlos):** la Fase 6 no usará esta célula de Hadley parametrizada, sino un núcleo dinámico propio (camino B; ver `DISENO_FASE_6.md` §7).

### 3.7 Estructura vertical de la humedad — DECIDIDO (opción H2, Carlos 04/10/2026)

**Opción H2:** el modelo guarda la **humedad específica de cada capa** (q_CL y q_TR), igual que guarda la temperatura de cada una.
- La evaporación entra en la CL; la convección (§3.8) sube vapor a la TR; la difusión C′ (§3.6) usa q_CL.
- **Saturación de la TR:** NUNCA con la temperatura media de la capa (la saturación es muy no lineal). Se integra q_s a lo largo del **perfil de temperatura reconstruido** de la TR (adiabático húmedo desde la CL hasta la tropopausa, §3.8, e isotermo por encima), ponderando por masa. Lo mismo, con su propio tramo, en la CL.
- **Condensación de gran escala:** a partir de una **humedad relativa crítica** de capa (< 1), porque una capa entera rara vez está saturada de media aunque partes de ella sí lo estén. Parámetro (c): ⚠️ valor y fuente por verificar; prueba de sensibilidad obligatoria.
- Depende de la decisión sobre la tropopausa (§3.8).

**Por qué H2** (y no H1, la recomendación del borrador v0): coherencia con las dos capas de temperatura de M3N; la humedad cerca de la superficie que necesita la C′ es propia del modelo y no sale de una forma impuesta; menos datos terrestres de tipo (c); y prepara la Fase 5.2, porque el efecto invernadero del vapor depende sobre todo de la humedad de la troposfera alta, que así calcula el modelo.

Descartadas: H1 (agua de la columna con forma vertical fija de Manabe y Wetherald 1967, de tipo c) y H3 (solo q_CL con la relación de Smith 1966, más c).

**Coste aceptado:** más complejidad y pruebas de estabilidad numérica específicas.

### 3.8 Convección húmeda y gradiente vertical (decisión 6) (a + b)

- **Gradiente adiabático húmedo en coordenadas de presión** ✅ (Bakhshaii y Stull 2013; implementación de MetPy `moist_lapse`):
  dT/dp = (1/p) · (R_d·T + L_v·r_s) / (c_pd + L_v²·r_s·ε / (R_d·T²))
  **La gravedad no aparece** en coordenadas de presión: vale para P3N sin adaptar. Solo al pasar a K/km entra la g de P3N (seco: g/c_p = 9,0 K/km).
- **Ajuste convectivo húmedo** (Manabe, Smagorinsky y Strickler 1965): la convección actúa cuando la diferencia CL−TR supera la que correspondería a un perfil adiabático húmedo **calculado en esa celda y en ese momento**. Sustituye al umbral fijo de 40,5 K, que salía del perfil de la Atmósfera Estándar terrestre.
- **Dónde acaba el perfil** (la TR incluye la estratosfera): **propuesta de la IA**, la tropopausa se sitúa donde el adiabático húmedo alcanza la **temperatura de piel** de una atmósfera gris, T_piel = 2^(−1/4)·T_e, con T_e = (OLR/σ)^(1/4) del propio modelo (Pierrehumbert 2010, cap. 4). Comprobación: en la Tierra, T_e ≈ 255 K da T_piel ≈ 214 K, cerca de la tropopausa observada. Sale del propio P3N, sin datos terrestres. ⚠️ Pendiente de verificar la referencia exacta y de confirmar. Limitación: ignora el calentamiento por ozono, que en P3N dependería de si hay oxígeno (worldbuilding).
- **Precipitación convectiva**: relajación tipo Betts-Miller simplificado (Frierson 2007): la humedad se relaja hacia una humedad relativa de referencia RH_ref con un tiempo τ, y el calor latente se suelta en la TR.
  - τ = 7200 s ✅ (valor por defecto en Isca `qe_moist_convection`; control en Frierson 2007b).
  - RH_ref: **0,8** por defecto en Isca ✅, **0,6** de control en Frierson (2007b) ✅, con experimentos de 0,8 a 0,95. **No hay un valor único**: es (c) y es **incierto** → prueba de sensibilidad obligatoria y documentada, NO ajuste para cuadrar resultados.
- **Condensación de gran escala**: Isca usa por defecto humedad relativa = 1 (`hc = 1.0` en `lscale_cond` ✅), pero Isca tiene muchas capas finas; con las dos capas gruesas de M3N hace falta una **humedad relativa crítica de capa < 1** (§3.7, ⚠️ pendiente de fuente). El calor se suelta en la capa donde se condensa. Sin reevaporación en la 5.1.
- **Validación clave:** en modo Tierra, el gradiente medio que salga debe acercarse a 6,5 K/km. Si sale, la física funciona; no se impone.

### 3.9 Lluvia o nieve (a + c)

Fracción de lluvia lineal con la temperatura del aire a 2 m entre T_toda_nieve y T_toda_lluvia, acotada entre 0 y 1. ✅ Forma verificada en el código del CLM (`repartition_rain_snow_one_col`). ⚠️ Los valores de los umbrales **no los he podido verificar** (en el CLM vienen de un fichero de parámetros). Físicamente dependen de si el copo se funde al atravesar el aire templado cercano al suelo, que es física (a) con el perfil de P3N. Pendiente.

### 3.10 Nieve en tierra (Riesgos 1 y 3)

- Reserva de nieve S (kg/m² de agua equivalente). Se funde con la energía sobrante cuando la superficie llega a 0 °C, con L_f, igual que el hielo marino.
- **Tope** S_max: lo que lo supere vuelve al océano como "descarga glaciar", **conservando agua y energía** (fundirlo cuesta L_f al océano). Referencia ✅: el CLM5 usa 10 000 mm por defecto y el CLM4.5 usaba 1000 mm (`h2osno_max`). El valor para M3N se justificará con el criterio de que más espesor ya no cambie el clima. Pendiente.
- El espesor por encima del tope **no entra en la convergencia**; sí la cobertura y la temperatura.
- **Albedo gradual** según la fracción cubierta (nunca un escalón; lección del intento perdido de la Fase 3).
- ⚠️ **Hallazgo P3N-específico — el albedo de la nieve depende de la estrella (a + b).** La nieve refleja muchísimo en el visible y bastante menos en el infrarrojo cercano. S3N (~5420 K) emite proporcionalmente más infrarrojo cercano que el Sol (~5772 K), así que **el albedo "de banda ancha" de la nieve en P3N es menor que en la Tierra** (Shields et al. 2013, ⚠️ pendiente de verificar). Lo mismo afecta al **albedo del hielo marino (0,65), que hoy es un valor terrestre**. Propuesta: calcular el albedo de la nieve y del hielo integrando su albedo espectral sobre el espectro de S3N. Necesita datos espectrales de una fuente que hay que verificar. Puede ir en la 5.1 para la nieve, y como corrección pendiente para el hielo marino.

### 3.11 Hielo marino

- Sublimación (§3.3).
- La nieve que cae sobre el hielo se suma a su masa (agua y energía conservadas).
- **PENDIENTE**: el aislamiento térmico de la nieve sobre el hielo. Motivación: el hielo de M3N sale demasiado grueso (12–23 m en P3N; 9 m frente a 2–3 m observados en modo Tierra).

---

## 4. Interruptores (con todos apagados, v2.4.1 bit a bit — prueba V0)

| | Nombre provisional | Contenido | Requiere |
|---|---|---|---|
| I10 | `vapor` | Evaporación, vapor, transporte, condensación y precipitación | I2, I3 |
| I11 | `suelo_y_nieve` | Cubo y nieve en tierra | I10 |
| I12 | `conveccion_humeda` | Ajuste convectivo húmedo (sustituye a I5) | I10, I5 |
| I13 | `presion_local_saturacion` | Opción B (§3.4) | I10 |

---

## 5. Parámetros

| Parámetro | Valor | Etiqueta | Fuente | Estado |
|---|---|---|---|---|
| L_v0 (a T₀ = 273,16 K) | 2,50084×10⁶ J/kg | a | MetPy | ✅ |
| L_f | 3,337×10⁵ J/kg | a | MetPy | ✅ |
| L_s0 | L_v0 + L_f | a | MetPy | ✅ |
| R_v | 461,52 J/kg/K (R/M_agua, M = 18,015268 g/mol) | a | MetPy | ✅ |
| c_pv | 1860 J/kg/K | a | MetPy | ✅ |
| c_pl | 4219,4 J/kg/K | a | MetPy | ✅ |
| c_pi | 2090 J/kg/K | a | MetPy | ✅ |
| e₀ | 611,2 Pa | a | MetPy | ✅ |
| ε | 0,62198 | a | M_agua/M_aire | ✅ |
| g | 9,053 m/s² | b | `parametros.py` | ✅ |
| p_s | 1 bar | b | Decisión de Carlos, 01/10 | ✅ |
| Viento (P3N) | 5 m/s provisional | b | Hasta la Fase 6 | ✅ decidido |
| Viento (modo Tierra) | Media real sobre el océano | c | Kent et al. (2013) u otra fuente | ⚠️ valor por verificar |
| C_HN océano | 1,1×10⁻³ | c | Large y Pond 1982; Smith 1988 (ya en v2.4.1) | ✅ en código |
| W_fc | 150 mm | c | Manabe 1969, vía Isca | ✅ |
| Umbral β | 0,75·W_fc | c | Isca | ✅ código / ⚠️ original |
| τ convección | 7200 s | c | Isca; Frierson 2007b | ✅ |
| RH_ref convección | 0,6–0,8 (incierto) | c | Isca; Frierson 2007b | ✅ el rango / sensibilidad obligatoria |
| RH crítica de condensación de gran escala (por capa) | < 1, por determinar | c | Isca usa 1,0 con capas finas (✅); para capas gruesas, ⚠️ fuente pendiente | ⚠️ |
| Umbrales lluvia/nieve | — | a + c | — | ⚠️ |
| S_max | — | — | Criterio físico (§3.10) | ⚠️ |
| Albedo de la nieve | — | a + b | Integración espectral con S3N | ⚠️ |
| D_atm (ahora sobre la energía estática húmeda, §3.6), D_oc | Se recalibran | c | Modo Tierra frente a Trenberth y Caron 2001; cruce con Hwang y Frierson 2010 | Se recalcula |

---

## 6. Numérica

- Evaporación y calor sensible implícitos y acoplados (§3.3).
- Variables nuevas: vapor (o W), agua del suelo, nieve, todas ≥ 0. **Un recorte a cero crearía agua: debe registrarse y salir en el balance.**
- **Salto de aceleración (Riesgo 2):** el vapor NO se extrapola (se reajusta en días); el agua del suelo SÍ, acotada a [0, W_fc]; la nieve solo por debajo del tope. El año del salto no cuenta para el balance de agua.
- **Convergencia:** se añade el agua del suelo (cambio anual en mm, con umbral a justificar); no el espesor de nieve sobre el tope.
- **Caché:** los módulos nuevos se añaden a `MODULOS_FISICA`.
- **Paso de tiempo:** prueba de sensibilidad con 450 s.

---

## 7. Validación

**V0** Con todos los interruptores apagados, v2.4.1 bit a bit.
**V1** Balance de agua anual: evaporación = precipitación ± Δalmacenamiento ± descarga glaciar (tolerancia a fijar, del orden del balance de energía).
**V2** Balance de energía < 10⁻³ (actual: ~10⁻⁵), con el calor latente incluido.
**V3** Paso de tiempo a la mitad: resultados sin cambios apreciables.
**V4** Salto de aceleración activado frente a desactivado: el mismo clima.
**V5 Modo Tierra** (con 24 h, viento real de la Tierra):
- Calor latente global ≈ 82 W/m² y sensible ≈ 21 W/m² (Wild et al. 2019; ya citado en DISENO_FASE2B).
- Precipitación ≈ 2,7–2,8 mm/día (GPCP; ⚠️ cifra exacta por verificar).
- Tiempo de residencia del vapor 8,5–8,9 días ✅ (van der Ent y Tuinenburg 2017).
- Agua precipitable global ⚠️ (~25 kg/m², por verificar).
- Transporte atmosférico máximo ~5,0 PW (Trenberth y Caron 2001), **ahora alcanzable** si el diagnóstico de la saturación es correcto. Reparto latente/sensible como validación independiente.
- Gradiente vertical medio ≈ 6,5 K/km, que debe **salir**, no imponerse.
- Lluvia: válida en latitudes medias y altas; en los trópicos se espera que falle (sin Hadley), y se documenta.

**V6** Sensibilidad a RH_ref (0,6 / 0,7 / 0,8) documentada.
**V7** P3N con prueba1: resultados, comparación con la v2.4.1 y limitaciones.

**Regla:** ningún parámetro se ajusta para esconder un error que pertenece a otra fase.

---

## 8. Exportación (m3n-clima v1, claves nuevas sin subir de versión)

Precipitación (total y nieve), evaporación, humedad del aire a 2 m o punto de rocío, agua del suelo y nieve, diarias. **Advertencia en los metadatos:** "precipitación tropical no válida hasta la Fase 6". H3N: vista "Precipitación y nubes" y climograma de Walter y Lieth (este último con la misma advertencia).

---

## 9. Dependencia de la rotación — DECIDIDA (Carlos y Ozan, 04/10/2026)

**Día solar de 19,84 h** (71 424 s), aplicado en la v2.4.3: hora de P3N = 1/24 del día (2976 s), paso de 992 s, año de 326,75 días y día sideral de 19 h 46 min 46 s. Equivale en circulación (Ω·a) a un día terrestre de ~27,0 h. D de la atmósfera escalado por (Ω_Tierra/Ω_P3N)² = 0,683. Pendiente: prueba de sensibilidad del paso (`prueba_paso.py`).

Notas originales:

- Con la v2.4.2 el código acepta cualquier día que cumpla dos condiciones: **múltiplo exacto del paso de 900 s** y, para tener registro horario (y exportación a H3N), **número entero de horas de 3600 s**. Dentro del intervalo de Carlos y Ozan (16,5–18,7 h) eso deja **17 h o 18 h**. Otros valores (p. ej. 17,6 h = 63 360 s, que no es múltiplo de 900 s) exigirían cambiar el paso de tiempo o redefinir la "hora" de P3N.
- **Fijarlo antes de la calibración de la 5.1** evita calibrar dos veces.
- D se calibra en modo Tierra (24 h) y se traslada a P3N con D ∝ 1/Ω² (Williams y Kasting 1997).
- Forma de la circulación (Fase 6): depende de Ω·a; día equivalente en la Tierra = día de P3N × 1,362 (a_Tierra/a_P3N = 6371/4676,7).
- Rango de validez razonable de M3N: ~16–48 h.

---

## 10. Limitaciones que se declaran desde ya

1. Lluvia tropical no válida hasta la Fase 6 (sin célula de Hadley).
2. Sin lluvia orográfica direccional hasta la Fase 6.
3. **Presión local solo en la saturación; opción C PENDIENTE.**
4. Calor latente constante en el balance de energía (§3.2). **PENDIENTE: pulir hacia la opción C** (entalpía completa del agua).
5. Emisividades fijas con vapor terrestre "de fábrica" hasta la 5.2.
6. Escorrentía instantánea al océano (sin ríos).
7. Cubo de suelo terrestre provisional hasta los biomas.
8. Sin aislamiento por la nieve sobre el hielo marino (PENDIENTE).
9. Viento fijo de 5 m/s en P3N hasta la Fase 6.
10. Tropopausa por temperatura de piel gris, sin ozono.

---

## 11. Decisiones de la IA pendientes de confirmar por Carlos

1. ~~Calor latente constante en el balance de energía; Ambaum solo para la saturación (§3.2).~~ ✅ **Confirmado por Carlos (04/10): opción A**, con la C pendiente.
2. ~~Vapor ligado a D_atm con una única K (§3.6).~~ ✅ **Decidido por Carlos (04/10): opción C′** (difusión de la energía estática húmeda cerca de la superficie, un solo coeficiente). Pendientes: regla de reparto vertical del calor que llega; si h incluye g·z.
3. ~~Estructura vertical de la humedad: H1 recomendada (§3.7).~~ ✅ **Decidido por Carlos (04/10): opción H2** (humedad propia en cada capa). Pendiente: humedad relativa crítica de la condensación de gran escala.
4. Tropopausa por temperatura de piel (§3.8).
5. Albedo espectral de la nieve con el espectro de S3N (§3.10), y si se corrige también el del hielo marino.
6. Versión v2.4.2 previa (§2).

---

## Referencias

Verificadas (✅) en este borrador:
- Ambaum, M. H. P. (2020). Accurate, simple equation for saturated vapour pressure over water and ice. *Q. J. R. Meteorol. Soc.* — vía la implementación de MetPy.
- Alduchov, O. A., y Eskridge, R. E. (1996). Improved Magnus form approximation of saturation vapor pressure. *J. Appl. Meteor.*, 35, 601–609.
- Bakhshaii, A., y Stull, R. (2013) — gradiente adiabático húmedo, vía MetPy.
- Bonan, D. B., et al. — código del modelo de balance de energía húmedo (github.com/dbonan/energy-balance-models).
- Frierson, D. M. W. (2007). The dynamics of idealized convection schemes and their effect on the zonally averaged tropical circulation. *J. Atmos. Sci.*, 64, 1959–1976.
- Frierson, D. M. W. (2007b). Convectively coupled Kelvin waves in an idealized moist general circulation model. *J. Atmos. Sci.*
- Hwang, Y.-T., y Frierson, D. M. W. (2010). Increasing atmospheric poleward energy transport with global warming. *GRL*, 37.
- Isca (ExeClim/Isca): `qe_moist_convection.F90`, `lscale_cond.F90`, `surface_flux.F90`, `idealized_moist_phys.F90`.
- CTSM/CLM5 (ESCOMP/CTSM): `atm2lndMod.F90`, `namelist_defaults_ctsm.xml`.
- MetPy (Unidata): `constants/default.py`, `calc/thermo.py`.
- Siler, N., Roe, G. H., y Armour, K. C. (2018). Insights into the zonal-mean response of the hydrologic cycle to global warming from a diffusive energy balance model. *J. Climate*, 31, 7481–7493.
- van der Ent, R. J., y Tuinenburg, O. A. (2017). The residence time of water in the atmosphere revisited. *HESS*, 21, 779–790.
- ECMWF (2002). Review of parametrization schemes for land surface processes.

Citadas pero ⚠️ pendientes de verificar en su fuente original antes de usarlas en código: Manabe (1969, la regla del 75 %), Manabe y Wetherald (1967), Manabe, Smagorinsky y Strickler (1965), Pierrehumbert (2010, temperatura de piel), Shields et al. (2013), Smith (1966), Kent et al. (2013), GPCP, agua precipitable global.

Ya usadas en fases anteriores: Large y Pond (1982), Louis (1982), Trenberth y Caron (2001), Wild et al. (2019), Williams y Kasting (1997).

# DISENO_FASE_5.2.md — Fase 5.2: nubes y radiación por bandas (y, después, la estrella)

**Estado:** BORRADOR DE INVESTIGACIÓN (12/10/2026). Recoge qué tiene que resolver la fase, qué hay hoy en M3N, las opciones con sus fuentes y las decisiones que tocan a Carlos. **Decidido el 12/10/2026 (Carlos, 07:38): la opción (B) del §3.1, por recomendación técnica, y MT_CKD como dato declarado (§2).** Lo demás del §5 sigue abierto.
**Convención:** ✅ = comprobado en la fuente (se dice qué se comprobó); ✅t = comprobados solo el título, la revista y los autores; ⚠️ = por comprobar; 🔶 = propuesta de la IA pendiente de Carlos; (a) física universal, (b) propiedad de P3N, (c) calibración terrestre.

---

## 0. Qué es la Fase 5.2 y por qué es la más importante de las que quedan

- **Qué añade:** (1) la radiación **por bandas** de longitud de onda, en lugar de la radiación gris de hoy; (2) el **ozono** y su calentamiento de la estratosfera; (3) las **nubes**: dónde se forman, cuánta agua tienen y cómo reflejan la luz y atrapan el infrarrojo. Después, **la estrella se recalibra** (decisión del 05/10/2026, reafirmada el 08/10).
- **Por qué es la más importante:** casi todos los sesgos medidos en las fases 6.3 y 6.4 apuntan aquí:

| Sesgo medido | Valor | Dónde | Lo que lo explica (lectura ya escrita) |
|---|---|---|---|
| Modo Tierra demasiado cálido | +5,16 K (19,13 °C frente a 13,97, CRU) | 3.18.1, `DISENO_FASE_6.3.md` §6.25 | sin nubes |
| P3N más cálido que el clima oficial | +7,1 K (21,57 frente a 14,45 °C) | `DISENO_FASE_6.4.md` §8.3 | sin nubes; estrella calibrada para el modelo de 2 capas |
| Albedo planetario | 0,139 (Tierra) y 0,144 (P3N), frente a 0,291 observado | §6.15 y §8.3 | sin nubes |
| Tropopausa tropical baja y caliente | 13,3 km y −44 °C (obs. 16,5 km y −80 °C) | §6.15 | radiación gris, sin ozono (`DISENO_FASE_2.3.md` §12) |
| Estratosfera de verano con el signo cambiado | +65 y +71 m/s (obs. −7 y −6) | §6.25 punto 5 | sin ozono |
| Vórtice estratosférico | 2–5 veces el real; en P3N 91 m/s | §6.25, §8.3 | radiación gris y sin ozono (5.2) + sin ondas de gravedad (6.6) |
| Banda de lluvias al sur, hielo antártico escaso | — | §6.25 puntos 3 y 7 | sin nubes del océano Austral (Hwang y Frierson 2013) |
| Ciclo del agua +25 % | 3,47 mm/día (obs. 2,67) | §6.25 | +5 K por Clausius-Clapeyron |

## 1. Qué hay hoy en M3N (leído en el código y los documentos, 3.20.0)

- **Infrarrojo:** dos flujos **gris**, τ(p) = τ₀·[f·(p/p₀) + (1 − f)·(p/p₀)⁴] (Frierson et al. 2006), con τ₀ y f (c) calibrados para dar los flujos **sin nubes** de Wild et al. (2019) en la columna media terrestre (3.18.1: τ₀ = 2,4450, f = 0,2326). El vapor del modelo **no** entra en el infrarrojo (I15, prototipo, no recomendado).
- **Luz:** transmisión media con τ_sw = 0,18 y el 73 % de lo atenuado absorbido por el aire (Wild et al. 2019), repartido en la vertical con peso ∝ Δ[(p/p₀)⁴]. Sin dependencia de la longitud de onda, salvo el albedo de la nieve y el hielo (I14), que usa dos bandas de Planck a 5420 K (44,2 % de la luz por debajo de 0,7 µm; el Sol, 48,8 %; `DISENO_FASE_5.1.md`).
- **Sin ozono, sin nubes, sin aerosoles.**
- **Agua condensada:** la condensación y la convección (Betts-Miller) precipitan **al instante**: no hay agua de nube en el modelo. Es la pieza que falta para unas nubes con agua propia.
- **Composición del aire:** como la Tierra (decisión del 11/10/2026, `DISENO_FASE_6.4.md` §7.3), con el aviso de que **vuelve a importar en la 5.2**.

## 2. El principio que condiciona todo: modelo propio (05/10/2026)

"El objetivo no es usar otros modelos, es construir el nuestro propio de manera exhaustiva y rigurosa." Aplicado a la radiación, la línea que propongo 🔶:
- **Datos de laboratorio y de observación** (posiciones e intensidades de las líneas de absorción, secciones eficaces del ozono, el continuo del vapor, el espectro de una estrella): **sí**. Son medidas del mundo, como CRU o ERA-Interim para validar o las constantes del CLM5 que ya se usan.
- **Las tablas precalculadas de otro esquema de radiación** (por ejemplo, los coeficientes de RRTMG): **no**, porque serían la física de otro modelo dentro de M3N.
- **El continuo del vapor (MT_CKD)** está en medio: es un modelo semiempírico, pero es la referencia con la que se ajustan todos y no hay alternativa medida completa. Propuesta: usarlo como dato, declarado (§3.1). **Decisión de Carlos.**

## 3. Las piezas, con sus opciones

### 3.1 Los gases en el infrarrojo y en la luz (la "radiación por bandas")

**El problema:** la absorción de un gas cambia miles de veces en pocos números de onda (las líneas). Calcularla línea a línea en cada columna y cada paso es inabordable; los modelos agrupan.

**Opciones:**
- **(A) Semi-gris de pocas bandas** (ventana atmosférica + resto; una o dos bandas de luz), con espesores calibrados (c). Barato; sigue siendo una calibración terrestre y no da bien la estratosfera ni el efecto del vapor de P3N.
- **(B) Distribución k correlacionada propia, calculada por M3N a partir de las líneas de HITRAN.** El método estándar de los modelos de clima: en cada banda, la absorción se reordena por intensidad (la función de distribución k) y se integra con unos pocos puntos de cuadratura ("puntos g"); cada punto se calcula como si fuera monocromático, con dos flujos. Método: Lacis y Oinas (1991, *J. Geophys. Res.* 96, 9027 ⚠️), Fu y Liou (1992, *J. Atmos. Sci.* 49, 2139 ⚠️). Datos: HITRAN2020 (Gordon et al. 2022, *J. Quant. Spectrosc. Radiat. Transf.* 277, 107949 ✅t) y el continuo MT_CKD (Mlawer et al. 2012, *Phil. Trans. R. Soc. A* 370 ✅t; código y coeficientes en `github.com/AER-RC/MT_CKD`, versiones v3.5 y v3.6 ✅, alcanzable desde el entorno de la IA).
- **(C) Usar tablas k de otro esquema** (RRTMG: Mlawer et al. 1997, *J. Geophys. Res.* 102, 16663 ⚠️). Más rápido de montar, pero choca con el principio del §2.

**Propuesta 🔶: (B).** Es la única coherente con "modelo propio" y la única que funciona igual de bien para P3N (otra estrella, otro vapor) que para la Tierra. Es también la más cara de construir:
1. **Cálculo línea a línea propio** (perfil de Voigt) de H₂O, CO₂, O₃ (y, si se decide, CH₄, N₂O, O₂) en una rejilla de presión y temperatura, más el continuo. Una sola vez, en el PC de Carlos (horas).
2. **Tablas k** por banda y punto g, con la mezcla de gases resuelta (método a elegir, ⚠️).
3. **Validación del cálculo línea a línea y de las tablas** frente a cálculos de referencia publicados: RFMIP (Pincus, Forster y Stevens 2016, *Geosci. Model Dev.*, flujos línea a línea de referencia en 100 perfiles ✅t) y CIRC (Oreopoulos et al. 2012, *J. Geophys. Res.* 117, D06118 ✅t, casos observados con cálculos de referencia). Criterio a fijar antes (orientativo: error en flujos sin nubes de unos pocos W/m² y en tasas de calentamiento de ~0,1–0,2 K/día).
4. **Dos flujos con dispersión** para la luz (nubes y Rayleigh): Toon et al. (1989, *J. Geophys. Res.* 94, 16287 ⚠️) o Meador y Weaver (1980 ⚠️).

**Coste del cálculo en el modelo (estimación propia, por medir):** con ~16 bandas de infrarrojo y ~14 de luz y ~8–16 puntos g por banda, del orden de 200–400 cálculos monocromáticos por columna. Con la radiación cada hora de P3N (cada 3 pasos, como hacen muchos modelos con la radiación cada hora ⚠️), y numba, el coste por año podría parecerse al del núcleo. **Hay que medirlo con un prototipo antes de comprometer el presupuesto** de "2 h, techo una noche" (06/10/2026).

**Requisito práctico:** HITRAN se descarga desde hitran.org con una cuenta gratuita, y desde el entorno de la IA no se llega (solo GitHub y PyPI). La descarga la haría Carlos en su PC.

### 3.2 La luz de S3N

- Hoy: cuerpo negro de 5420 K (solo para el albedo de la nieve). Con bandas, **el espectro de la estrella entra en todo**: el reparto de la luz entre el ultravioleta (ozono), el visible (Rayleigh, nubes) y el infrarrojo cercano (vapor).
- **Opciones:** (a) cuerpo negro de Planck; (b) un **espectro sintético** de una atmósfera estelar real con esa temperatura, por ejemplo la biblioteca PHOENIX (Husser et al. 2013, *Astron. Astrophys.* 553, A6 ✅t). Una estrella real tiene **mucho menos ultravioleta** que un cuerpo negro de su temperatura (absorción de sus propias líneas), y el ultravioleta es lo que fabrica el ozono.
- **El ultravioleta lejano** depende de la **actividad** de la estrella (cromosfera), que no sale de la temperatura. Es una decisión de worldbuilding. Rugheimer et al. (2013) muestran que el ozono no varía de forma monótona con la temperatura de la estrella precisamente por eso (§3.3).
- **Propuesta 🔶:** (b) para el visible y el infrarrojo, más una decisión explícita sobre el ultravioleta de S3N (estrella "tranquila como el Sol" u otra cosa).
- Rayleigh: Bucholtz (1995, *Appl. Opt.* 34, 2765 ⚠️).

### 3.3 El ozono

- **Por qué importa:** calienta la estratosfera absorbiendo ultravioleta. Sin él, la tropopausa y la estratosfera salen mal (§0), y probablemente el vórtice.
- **Qué dice la literatura para una estrella como S3N** (Rugheimer et al. 2013, *Astrobiology*, doi:10.1089/ast.2012.0888, arXiv:1212.2638 ✅, leído el 12/10/2026):
  - modelo 1D con fotoquímica (55 especies, 220 reacciones), 1 bar, 21 % de O₂ y 355 ppm de CO₂, estrellas de 4250 a 7000 K;
  - **columna de ozono** (Tabla 2): Sol 5,3·10¹⁸ cm⁻²; **5500 K (G8V), la más parecida a S3N (5420 K): 3,2·10¹⁸ cm⁻²**, el 60 % de la del Sol; 5250 K (K0V): 4,1·10¹⁸;
  - con estrellas frías, estratosferas casi isotermas (sin la inversión fuerte de la Tierra) y tropopausa y estratopausa más altas;
  - el ozono **no** sigue la temperatura de forma monótona, porque depende del ultravioleta de cada estrella concreta (se usaron espectros observados con IUE).
  - Ojo: en su modelo, la Tierra con el Sol tiene ~197 unidades Dobson (cálculo propio de la conversión), por debajo de las ~300 reales. Vale la **proporción**, no la cifra.
- Segura et al. (2003, *Astrobiology* 3(4), 689–708, doi:10.1089/153110703322736024 ✅t) es el trabajo de referencia anterior (F2V, G2V, K2V).
- **Opciones:**
  - (a) **ozono prescrito**: un perfil vertical con su columna, escalado para S3N a partir de esa literatura; la forma vertical, de la misma fuente. Barato y declarado.
  - (b) **fotoquímica propia**: Chapman más los ciclos catalíticos (NOx, HOx, ClOx). Chapman solo da demasiado ozono; los ciclos exigen decidir la química de P3N (nitrógeno, vida, volcanes). Es un proyecto en sí mismo.
  - (c) esquema lineal interactivo (tipo Cariolle ⚠️): derivado de un modelo químico de la Tierra; choca con el principio.
- **Propuesta 🔶:** (a) en la 5.2, con la decisión de worldbuilding sobre el O₂ y el ultravioleta de S3N; (b) como posible fase futura.

### 3.4 Las nubes

**Cuatro piezas:**
1. **Agua de nube.** Hoy no existe (todo precipita al instante). Propuesta 🔶: agua de nube (líquida y hielo) como trazador nuevo del núcleo, con la condensación de gran escala que la forma y una conversión a lluvia (autoconversión, tipo Sundqvist et al. 1989 ⚠️).
2. **Fracción de nube.** Opciones: (a) diagnóstica por humedad relativa (Sundqvist et al. 1989, *Mon. Weather Rev.* 117 ⚠️; Slingo 1987, *Q. J. R. Meteorol. Soc.* 113, 899 ⚠️); (b) pronóstica (Tiedtke 1993, *Mon. Weather Rev.* 121, 3040 ⚠️); (c) basada en la distribución de la humedad dentro de la celda (Tompkins 2002 ⚠️). Con celdas de 5° (400–550 km), cualquier esquema es una parametrización fuerte.
3. **Nubes convectivas:** Betts-Miller no da agua de nube; hay que diagnosticarlas (por ejemplo, a partir de la lluvia convectiva, Slingo 1987 ⚠️) o cambiar de convección (no propuesto).
4. **Óptica y solapamiento:** gotas (Slingo 1989, *J. Atmos. Sci.* 46, 1419 ⚠️) y hielo (Fu 1996, *J. Climate* 9, 2058 ⚠️; Fu et al. 1998 ⚠️) por banda; **radio efectivo** de las gotas (depende de los aerosoles: decisión, porque los aerosoles de P3N no están definidos); solapamiento exponencial-aleatorio (Hogan e Illingworth 2000, *Q. J. R. Meteorol. Soc.* 126, 2903 ✅t) y McICA para combinar nubes parciales con las bandas sin multiplicar el coste (Pincus, Barker y Morcrette 2003, *J. Geophys. Res.* 108(D13) ✅t).

**Riesgo declarado:** las nubes son la mayor fuente de incertidumbre de los modelos de clima reales. Habrá parámetros (c) que se calibran en el modo Tierra (por ejemplo, el umbral de humedad o la autoconversión). La regla de siempre: se fijan con la Tierra **antes** de mirar P3N, y nunca para "arreglar" P3N.

### 3.5 Acoplamiento

- Nuevos trazadores del núcleo: agua líquida y hielo de nube (el transporte positivo y conservativo de la 6.3 ya existe).
- Radiación cada N pasos (a decidir con lo medido) y calentamiento constante entre llamadas, con el cierre de energía de siempre.
- Interruptores nuevos I17 en adelante; con todos apagados, idéntico bit a bit a la versión anterior (norma de siempre).

## 4. Validación (modo Tierra, criterio a fijar antes)

| Qué | Referencia | Valor (estado) |
|---|---|---|
| Flujos en el tope, todo cielo | CERES EBAF Ed4.0 (Loeb et al. 2018, *J. Climate* 31, 895, doi:10.1175/JCLI-D-17-0208.1), Tabla 5, julio 2005–junio 2015 | reflejada 99,1; OLR 240,1 W/m² ✅ (leído en el artículo el 12/10 y coincide con `referencias_tierra.REF_CERES`) |
| Flujos en el tope, cielo despejado | ídem | reflejada 53,3; OLR 268,1 W/m² ✅ |
| Efecto radiativo de las nubes | derivado de las dos filas anteriores | luz −45,8; infrarrojo +28,0; neto −17,8 W/m² (cálculo propio con valores ✅) |
| Radiación sin nubes, perfiles | RFMIP (Pincus et al. 2016) y CIRC (Oreopoulos et al. 2012) | por leer ⚠️ |
| Agua líquida de las nubes | MAC-LWP (Elsaesser et al. 2017, *J. Climate* 30 ✅t) | por leer ⚠️ |
| Fracción de nubes y su altura | CALIPSO-GOCCP, ISCCP | fuente y cifras por elegir y comprobar ⚠️ |
| Tropopausa, vientos, Hadley | ERA-Interim (ya en `referencias_tierra.py`) | ✅ |
| Lo de la 6.3 que debería mejorar | temperatura, banda de lluvias, hielo, estratosfera | §6.25 |

## 5. Decisiones que tocarán a Carlos (para cuando se diseñe en firme)

1. **Método de los gases:** (A), (B) propuesta o (C) (§3.1).
2. **MT_CKD como dato** (§2).
3. **Gases de P3N:** cuánto CO₂ (y si hay CH₄, N₂O) y cuánto O₂. Worldbuilding; hoy "como la Tierra", sin cifra fijada.
4. **Espectro de S3N:** cuerpo negro o espectro sintético, y **su ultravioleta** (actividad) (§3.2).
5. **Ozono:** prescrito escalado (propuesta) o fotoquímica propia (§3.3).
6. **Aerosoles y radio de las gotas** (§3.4).
7. **Ambición de las nubes:** diagnósticas o pronósticas (§3.4).
8. **Presupuesto de coste**, cuando el prototipo lo mida (§3.1).

## 6. Siguientes pasos de la investigación (sin código todavía)

1. Leer y comprobar las fuentes marcadas ⚠️, empezando por las que deciden: Lacis y Oinas 1991, Fu y Liou 1992, Mlawer et al. 1997 (como método, no como tablas), RFMIP y CIRC (qué datos de referencia se pueden descargar y en qué formato).
2. Medir en el PC el coste de un cálculo línea a línea pequeño (un gas, una banda) para dimensionar (B).
3. Elegir las referencias de nubes del modo Tierra con cifras comprobadas.
4. Con eso, el diseño en firme y las decisiones del §5.

## Referencias citadas

✅ = leído o comprobado el 12/10/2026 (se dice qué); ✅t = título, revista y autores; ⚠️ = de memoria, por comprobar.
- Bucholtz 1995 ⚠️ · Elsaesser et al. 2017 ✅t · Frierson, Held y Zurita-Gotor 2006 (ya verificada en la Fase 2.3) · Fu 1996 ⚠️ · Fu y Liou 1992 ⚠️ · Fu et al. 1998 ⚠️ · Gordon et al. 2022 (HITRAN2020) ✅t · Hogan e Illingworth 2000 ✅t · Husser et al. 2013 ✅t · Hwang y Frierson 2013 (ya citada en la 6.3) · Lacis y Oinas 1991 ⚠️ · Loeb et al. 2018 ✅ (Tabla 5) · Meador y Weaver 1980 ⚠️ · Mlawer et al. 1997 ⚠️ · Mlawer et al. 2012 (MT_CKD) ✅t · Oreopoulos et al. 2012 ✅t · Pincus, Barker y Morcrette 2003 ✅t · Pincus, Forster y Stevens 2016 ✅t · Rugheimer et al. 2013 ✅ (métodos y Tabla 2) · Segura et al. 2003 ✅t · Slingo 1987 ⚠️ · Slingo 1989 ⚠️ · Sundqvist et al. 1989 ⚠️ · Tiedtke 1993 ⚠️ · Tompkins 2002 ⚠️ · Toon et al. 1989 ⚠️ · Wild et al. 2019 (ya verificada).

## 7. Decisiones técnicas tomadas y lo construido (12/10/2026)

Carlos delega las decisiones técnicas que no son de worldbuilding ni de principio, siempre con fuentes (06/10 y 09/10/2026). Las de esta lista son de ese tipo; las del §5 siguen siendo suyas. Dependían de la opción (B) del §3.1, **aprobada por Carlos el 12/10/2026**.

| # | Decisión técnica | Por qué | Fuente |
|---|---|---|---|
| T1 | Perfil de línea de **Voigt**, calculado con la función de Faddeeva (`scipy.special.wofz`) | Es el perfil estándar de los cálculos línea a línea en la troposfera y la estratosfera (Doppler + presión) | definición de HITRAN ⚠️ |
| T2 | Dependencia de las líneas con la temperatura y la presión **como la define HITRAN** (intensidad con la energía del nivel bajo y la función de partición; anchura con n_air; desplazamiento con δ_air) | Es el convenio con el que están medidos los datos: usar otro mezclaría convenios | Šimečková et al. 2006 ⚠️ |
| T3 | Funciones de partición **TIPS** como tabla de datos (no como librería dentro de M3N) | Son datos calculados de constantes moleculares; comprobado hoy, sin red: Q(H₂O, 296 K) = 174,58, Q(CO₂) = 286,09, Q(O₃) = 3475 | Gamache et al. 2017 ⚠️; valores leídos con `hitran-api` 1.3.0.0 (licencia MIT, PyPI) |
| T4 | Para el vapor, líneas cortadas a ±25 cm⁻¹ y, si así lo exige MT_CKD, sin el "pedestal" | Tiene que ser el mismo convenio con el que se definió el continuo, o se contaría dos veces | ⚠️ por comprobar en la documentación de MT_CKD antes de usarlo; de momento es una opción del código, apagada |
| T5 | Infrarrojo de cada punto g con el **mismo solucionador de dos flujos que ya tiene M3N** (fuente de Planck lineal en τ, `DISENO_FASE_2.3.md` §2) | Ya está verificado (isotermo exacto, convergencia con las capas); se reutiliza en lugar de escribir otro | ya verificado en la 2.3 |

**Construido (prototipo; en el repositorio desde la 3.21.0, sin conectar al modelo):** `fase52_lineas.py`, el cálculo línea a línea:
- intensidad a cualquier temperatura;
- anchuras de presión y Doppler;
- perfil de Voigt;
- coeficiente de absorción de una lista de líneas.

**Probado** (`test_fase52_lineas.py`, 5 pruebas, todas pasan), con casos de solución conocida:
- el perfil tiene área 1;
- en el pico coincide con el de Lorentz cuando el Doppler es despreciable, y con el de Gauss cuando lo es la presión (error < 10⁻⁶);
- la semianchura de Gauss es exacta;
- la anchura total coincide con la fórmula de Olivero y Longbothum (1977), dentro del 0,05 %;
- la intensidad no cambia a la temperatura de referencia y sigue la ley de Boltzmann;
- la anchura Doppler del CO₂ es la de la fórmula;
- integrado, el coeficiente devuelve la suma de las intensidades, salvo las alas cortadas (< 0,4 %).

**Lo siguiente que se puede construir sin datos reales:**
1. El paso de los coeficientes a **tablas k** por banda (ordenar dentro de la banda, cuadratura en g), probado con un modelo de banda con solución analítica.
2. El **dos flujos con dispersión** para la luz (delta-Eddington ⚠️), probado con casos de referencia publicados.

**Lo que necesita a Carlos:**
- La elección de (B) y del §2.
- Más adelante, descargar en su PC las líneas de HITRAN con una cuenta gratuita, porque desde el entorno de la IA no se llega a hitran.org (comprobado: sin respuesta).

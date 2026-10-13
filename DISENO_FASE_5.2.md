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
| T4 | Para el vapor, líneas cortadas a ±25 cm⁻¹ y sin el "pedestal" (su valor de Lorentz a ±25 cm⁻¹), como exige MT_CKD | Tiene que ser el mismo convenio con el que se definió el continuo, o se contaría dos veces | ✅ hitran.org/mtckd (comprobado el 12/10/2026); en el generador de tablas k, encendido para el vapor |
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

**Construido el 12/10/2026 (3.22.0, sin conectar):** `fase52_k.py`, las **tablas k** por banda:
- ordenar los coeficientes dentro de la banda;
- la cuadratura en g;
- tablas de varias capas con la hipótesis correlacionada;
- la transmisión de la banda, línea a línea y con la tabla.

**Decisión técnica T6 (provisional):** 16 puntos de Gauss repartidos en tramos de g cada vez más estrechos hacia g = 1 (donde están los centros de las líneas). Es la práctica habitual. El reparto definitivo se fijará midiendo el error de los flujos frente al cálculo línea a línea (§4), no copiando el de otro modelo.

**Probado** (`test_fase52_k.py`, 5 pruebas):
- **En una capa, la tabla es exacta salvo la cuadratura.** Con los 16 puntos, el error de la transmisión de una línea de Lorentz es < 2·10⁻⁴ en todo el rango, de línea débil a saturada.
- **El cálculo línea a línea reproduce la fórmula exacta de Ladenburg y Reiche** para la anchura equivalente de una línea aislada, con un acuerdo de 10⁻⁸ a 10⁻¹¹ (regímenes débil, intermedio y fuerte), contando las alas que quedan fuera de la banda.
- **La hipótesis correlacionada es exacta** con una línea a dos presiones (el orden de los números de onda es el mismo en las dos capas).
- **Su límite conocido:** con líneas que no se corresponden de una capa a otra, sobreestima la transmisión (0,961 frente a 0,944 en el caso de prueba). Por eso el error se medirá contra el cálculo línea a línea, como pide el §4.

**Construido el 12/10/2026 (3.22.0, sin conectar):** `fase52_onda_corta.py`, la **luz en una capa que absorbe y dispersa** (nubes, aerosoles, aire), con dos flujos delta-Eddington.

**Decisión técnica T7 (provisional):** ecuaciones de dos flujos de Meador y Weaver (1980) con los coeficientes de Eddington y el escalado delta (Joseph, Wiscombe y Weinman 1976) ⚠️. Es el método habitual en los modelos de clima para la luz. La capa se resuelve exacta para esas ecuaciones con un sistema de 2×2, en lugar de copiar fórmulas cerradas. Lo que llega al suelo se reparte en el haz directo de verdad (sin escalar) y el difuso.

**Probado** (`test_fase52_onda_corta.py`, 21 pruebas):
- **Sin dispersión** es exactamente la ley de Beer.
- **La energía se conserva** (reflejada + absorbida en la capa + absorbida en el suelo = 1, con error < 10⁻⁹) con suelo negro, gris y blanco. Sin absorción, la capa no absorbe nada.
- **Frente a un cálculo de referencia de ordenadas discretas** con 32 direcciones (PythonicDISORT 1.8, una implementación en Python de DISORT, licencia MIT ⚠️; convergido: con 16 direcciones cambia < 10⁻⁴), en 8 casos: nubes con el sol alto, medio y bajo; nube fina y gruesa; aerosol; dispersión casi isótropa; capa muy absorbente.
  - El error es el conocido del delta-Eddington (King y Harshvardhan 1986 ⚠️): ≤ 0,02 en la reflexión y la transmisión con el sol alto o capas finas, y hasta **0,038** con el sol bajo sobre una nube (μ0 = 0,2).
  - La tabla de referencia está en la prueba, y se recalcula si PythonicDISORT está instalado (no hace falta para M3N).

**La columna entera (3.22.0):** `fase52_onda_corta.columna`. Todas las capas se resuelven a la vez y de forma exacta para las ecuaciones de dos flujos: en cada capa, su solución general más la del haz que le llega; la continuidad de los flujos en cada separación; y las condiciones de arriba y del suelo. Es un sistema de 2N ecuaciones, equivalente a «sumar capas» pero sin aproximaciones. Probado (9 pruebas más):
- con una capa es exactamente `capa`;
- partir una capa en cinco no cambia nada (10⁻¹²);
- 40 capas con una nube de espesor óptico 100, sin desbordamientos y con la energía conservada;
- tres columnas (aire, nube y aerosol, con suelo negro y gris y con el sol a 53° y a 73°) frente a ordenadas discretas de 32 direcciones: la reflexión y la transmisión, dentro de 0,04 (medido: hasta 0,034 con μ0 = 0,3); la absorción de cada capa, dentro de 0,02 (medido: hasta 0,013).

**Datos de HITRAN descargados (PC de Carlos, 12/10/2026, 20:53; `herramientas/descargar_hitran.py` con hitran-api 1.3.0.0):** en `~/Documentos/B3N/datos/hitran` (390 MB, fuera del repositorio), de 0 a 50 000 cm⁻¹, todos los isotopólogos que conoce hitran-api:

| Gas | Líneas | Hasta (cm⁻¹) | sha256 del .data |
|---|---|---|---|
| H₂O | 490 609 | 41 999,7 | `2b2b6558…a4f552` |
| CO₂ | 544 002 | 17 696,9 | `6bf399eb…6b11ed` |
| O₃ | 487 446 | 6 996,7 | `33c2b392…52d9df` |
| N₂O | 516 124 | 14 017,1 | `fe1da5e3…30c8ad` |
| CH₄ | 473 522 | 13 922,6 | `cc083881…8ea059` |
| O₂ | 22 076 | 49 999,0 | `c4412395…c035d0` |

- Comprobado que ningún gas está cortado: cada uno acaba donde acaban sus líneas en HITRAN, salvo el O₂, que llega al límite pedido; lo que queda por encima de 50 000 cm⁻¹ es ultravioleta lejano, que solo importa en la alta atmósfera.
- Están todos los isotopólogos.
- Las huellas completas están en su `procedencia.json`. La tabla TIPS está en `tips.json`, de 60 a 400 K.

**Datos que faltan** (no son líneas):
- las secciones eficaces del ozono en el ultravioleta y el visible (Hartley-Huggins y Chappuis), que calientan la estratosfera;
- la absorción inducida por colisiones (O₂–O₂, N₂–N₂);
- el continuo MT_CKD del vapor;
- las propiedades ópticas de las gotas y el hielo de las nubes.

**Lector propio y bandas (3.22.0):** `fase52_hitran.py`.
- **Lee los archivos de HITRAN** (formato .par de 160 caracteres) sin depender de hitran-api, como pide T3. Entiende los códigos de isotopólogo ('0' = 10, 'A' = 11…) y el formato raro de algunas intensidades.
- **Probado** (3 pruebas): lee bien cada columna y da exactamente los mismos valores que el lector oficial, línea a línea.
- **Decisión técnica T8 (provisional): bandas.** Se adoptan los **límites** de las bandas de RRTMG: 16 en el infrarrojo, de 10 a 3250 cm⁻¹, y 14 en la luz, de 820 a 50 000 cm⁻¹ (Mlawer et al. 1997; Iacono et al. 2008; comprobados el 12/10 con el código de RRTMG en CESM, `radconstants.F90` ✅: la primera versión ponía 2380 en lugar de 2390 cm⁻¹ entre dos bandas del infrarrojo, ya corregido). Están elegidos por dónde absorbe cada gas y son de los más usados. No se adoptan sus tablas: las de M3N se calculan con HITRAN y para la atmósfera de P3N.
- `herramientas/descargar_hitran.py` guarda además `isotopologos.json`, con la masa de cada isotopólogo (necesaria para la anchura Doppler) y su abundancia terrestre (la que ya llevan dentro las intensidades: P3N se supone con la composición isotópica de la Tierra). Para tenerlo basta volver a lanzar el programa: lo ya descargado no se repite.

**El generador de las tablas k (3.22.0):** `fase52_generador.py` y `herramientas/generar_tablas_k.py`. Hacerlo «a lo bruto» (Voigt en una ventana de ±25 cm⁻¹ con el paso que exige la anchura Doppler en la alta atmósfera) costaba **~9 ms por línea** y serían cientos de horas. Con estas decisiones técnicas provisionales, cada una con su prueba de error, baja a **5–16 s por cálculo de una banda con 40 000 líneas** en el entorno de la IA:
- **T12, rejilla adaptada (no uniforme):**
  - una base de paso 0,01 cm⁻¹;
  - alrededor de cada línea estrecha (semianchura de Voigt < 0,04 cm⁻¹), puntos con paso w/4 hasta 2w, y después cada uno un 10 % más lejos;
  - cada punto pesa lo que mide su celda, y las tablas k usan esos pesos;
  - en la alta atmósfera, de 30 a 50 veces menos puntos que una rejilla uniforme.
- **T13, núcleo y alas:**
  - Voigt a menos de 50 anchuras Doppler del centro;
  - Lorentz en el resto del núcleo, hasta max(50 semianchuras, 0,25 cm⁻¹), en la rejilla fina;
  - las alas, hasta el corte de 25 cm⁻¹, en una rejilla gruesa de 0,02 cm⁻¹ interpolada;
  - junto a cada borde, el ala exacta.
  - Opción del **pedestal** para el vapor (T4). Comprobado ✅ en hitran.org/mtckd: MT_CKD define la contribución local de cada línea como su forma de Lorentz menos su valor a ±25 cm⁻¹ del centro, y el continuo se suma a eso.
- **T14, descarte de líneas débiles:**
  - se descarta una línea si su absorción con la columna máxima posible del gas es < 10⁻⁷ cm⁻¹;
  - la suma de lo descartado, que es una cota de lo perdido, queda guardada en cada casilla de la tabla.
- **T15, rejilla de la tabla:**
  - 36 presiones de 1050 a 0,01 hPa (equiespaciadas en log p);
  - 7 temperaturas de 150 a 330 K;
  - para el vapor, fracciones molares propias 0; 0,01 y 0,04 (por encima de 100 hPa, solo la 0, copiada);
  - las 29 bandas distintas de T8.
  - El error de interpolar en ella se medirá en la validación (§4).

**Probado** (`test_fase52_generador.py`, 15 pruebas):
- **Frente al cálculo de referencia** (Voigt completo en rejilla uniforme fina), con 150 líneas al azar y con y sin pedestal, en seis casos de la superficie a 10⁻⁵ atm: la transmisión de la banda difiere < 2·10⁻⁴ (medido: hasta 10⁻⁴), y la de la tabla k < 5·10⁻⁴ (medido: hasta 3,3·10⁻⁴, incluida la cuadratura).
- **La cota del descarte** es de verdad una cota de lo perdido.
- **Una casilla con dos isotopólogos** coincide con hacerlo a mano.
- **Una banda sin líneas** da k = 0.
- **El programa completo,** probado con datos sintéticos de prueba:
  - avance;
  - cortar y continuar, porque guarda cada casilla terminada;
  - juntar las tablas.

**El continuo del vapor (3.22.0):** `fase52_continuo.py`, MT_CKD 4.3 como dato (decisión de Carlos del 12/10).
- **El archivo:** `absco-ref_wv-mt-ckd.nc` de AER, fijado por su commit (`679f0a4`) y comprobado por su huella sha256 (`69944eb8…`). Lo descarga `descargar_hitran.py`.
- **La licencia de AER** permite el uso científico y redistribuirlo con su aviso de copyright, pero no incorporarlo a software comercial. Por eso no va en el repositorio.
- **Las fórmulas** son las del módulo oficial (`mt_ckd_h2o_module.f90`, leído el 12/10), traducidas sin cambios: escalado con la densidad y la fracción molar, dependencia del propio con la temperatura, término de radiación e interpolación de 4 puntos.
- **Probado** (3 pruebas): reproduce la salida de referencia del propio MT_CKD con un acuerdo de 5·10⁻⁵ (es la precisión con la que está impresa); el término de radiación; la interpolación; y el escalado.
- **En el generador,** el continuo se suma a las líneas del vapor (sin pedestal) antes de ordenar, en la misma tabla.

**Mie propio, para las gotas de las nubes (3.22.0):** `fase52_mie.py`.
- **Decisión técnica T16 (provisional):** las propiedades ópticas de las gotas de agua se calculan con la teoría de Mie y los índices de refracción medidos del agua, en lugar de copiar una parametrización. El hielo no es esférico y tendrá su propio tratamiento.
- **El algoritmo** es BHMIE (Bohren y Huffman 1983) con un arranque más holgado de la recurrencia. Al medirlo, el arranque original no convergía del todo con el agua casi transparente y gotas grandes: 5·10⁻⁵ de error en Q_ext a x = 100.
- **La población de gotas** sigue la distribución gamma de Hansen y Travis (1974) ⚠️, con su radio efectivo y su varianza.
- **Probado** (14 pruebas):
  - 10 casos frente a miepython 3.3 (Prahl, licencia MIT ⚠️; valores en la prueba), con un acuerdo de 10⁻⁹;
  - el límite de Rayleigh;
  - la paradoja de la extinción (Q_ext → 2);
  - una nube de r_e = 10 µm: en el visible da Q_ext ≈ 2,09, ω ≈ 1 y g ≈ 0,865, lo típico de las nubes de agua; a 3 µm el agua absorbe (ω ≈ 0,52).
- ~~Faltan los índices de refracción del agua y del hielo~~: hechos el 13/10 (T17, abajo).

**Índices de refracción del agua y del hielo (3.23.0):** `fase52_indices.py` y el archivo `datos/indices_refraccion_agua_hielo.npz` (0,8 MB, en el repositorio).
- **De dónde salen:** la base de datos de refractiveindex.info (Polyanskiy 2024, *Sci. Data* 11, 94 ⚠️), de **dominio público** (CC0), en GitHub, con el commit fijado (`c5c2f18`) y la huella sha256 de cada archivo. `herramientas/preparar_indices_refraccion.py` rehace el `.npz` y da el mismo archivo, byte a byte.
- **Decisión técnica T17 (provisional):**
  - **Agua líquida:** Rowe, Fergoda y Neshyba 2020 (*J. Geophys. Res. Atmos.* 125, e2020JD032624 ✅: resumen leído el 13/10; índices a 240, 253, 263 y 273 K de 1 a 15 000 cm⁻¹, k de la literatura entre 238 y 298 K y n por Kramers-Kronig). Como nodo de 298 K, Hale y Querry 1973 (*Appl. Opt.* 12, 555 ⚠️), el dato clásico a 25 °C. La temperatura se interpola entre 240 y 298 K.
  - **Por qué con temperatura:** las nubes de agua subfundida son comunes y en el infrarrojo el índice cambia mucho. Medido en los datos: k a 12 µm vale 0,31 a 240 K y 0,20 a 298 K. Rowe et al. muestran que ignorarlo sesga los flujos de esas nubes (hacia arriba entre 250 y 580 cm⁻¹ y hacia abajo entre 710 y 940 cm⁻¹).
  - **Por qué Hale y Querry y no Segelstein 1981 a 298 K:** comparado con los datos de Rowe et al., Hale y Querry sigue su tendencia con la temperatura (por ejemplo, n a 10 µm: 1,202 a 253 K, 1,209 a 273 K y 1,218 a 298 K); Segelstein queda fuera de ella (1,193).
  - **Por debajo de 0,667 µm**, donde el agua casi no absorbe, solo Hale y Querry. **Por encima de 200 µm**, solo Rowe et al. (y por encima de 273 K, el valor de 273 K): ahí el salto llega a 0,13 en n, pero por debajo de 50 cm⁻¹ está menos del 0,1 % de la emisión térmica (medido).
  - **Hielo:** Warren y Brandt 2008 (*J. Geophys. Res.* 113, D14220), a 266 K, sin temperatura (limitación declarada ⚠️). Comprobado frente a la tabla impresa de los autores en 9 longitudes de onda: idéntico.
  - Interpolación lineal en log λ para n y log k, y lineal en la temperatura.

**Óptica de las gotas por banda (3.23.0):** `fase52_optica_nubes.py`, con Mie propio y T17.
- Para una población de gotas (gamma de Hansen y Travis), da por banda el coeficiente de extinción másico β (m²/kg; el espesor óptico de una nube es β por su agua líquida), el albedo de dispersión simple ω y el factor de asimetría g.
- **Decisión técnica T18 (provisional):** promedio "fino" en la banda (Edwards y Slingo 1996 ⚠️), pesado con la fuente: la luz de la estrella en la onda corta y Planck a la temperatura de la nube en el infrarrojo. Es exacto para capas finas; el error con nubes gruesas se medirá en la validación (§4).
- **Probado** (18 pruebas, `test_fase52_indices.py`): el archivo de datos; el hielo igual que la tabla de Warren y Brandt; el agua igual que sus datos en los nodos; la interpolación en temperatura; los saltos de las costuras, medidos; en el visible, β ≈ 157 m²/kg con r_e = 10 µm (Q_ext ≈ 2,09, la fórmula de Stephens); las fórmulas del promedio; y que el efecto de la temperatura va en el sentido que describen Rowe et al.
- **Pendiente:** el radio efectivo de las gotas (depende de los aerosoles, ficha §6.4) y el hielo de las nubes (no esférico).

**Coste estimado en el PC (orientativo):** del orden de **una noche** (~6–9 h con 11 hilos), dominada por el vapor. `--prueba` lo mide primero con 4 bandas y 3 alturas. Antes hay que volver a lanzar `descargar_hitran.py` para que escriba `isotopologos.json`; no repite la descarga.

**El infrarrojo por bandas y el solapamiento de los gases (3.23.0, 13/10/2026):** `fase52_infrarrojo.py`, `fase52_referencia.py`, `fase52_mapa.py`, `herramientas/medir_solapamiento.py` y el generador nuevo.

- **Decisión técnica T19 (provisional): factor de difusividad 1,66** (Elsasser 1942 ⚠️): el flujo difuso de una capa de espesor óptico vertical τ se transmite como e^(−1,66 τ). Es lo habitual en el infrarrojo (RRTMG, ecRad ⚠️). Su error frente al cálculo con varias direcciones se medirá en la validación (§4).
- **Decisión técnica T20 (provisional):** la emisión que cae fuera de las bandas (por debajo de 10 cm⁻¹ y por encima de 3250 cm⁻¹) se suma a la primera y a la última banda, para que las bandas sumen exactamente σT⁴ y la energía cierre. Medido: fuera queda menos del 0,05 % de σT⁴ entre 150 y 330 K.
- **El dos flujos** es el mismo de M3N (T5): fuente de Planck lineal en τ. Probado: igual que las recurrencias de `fase30_multicapa`; exacto en una atmósfera isoterma; partir cada capa en 7 no cambia nada (10⁻¹²).

**Decisión técnica T9: el solapamiento de los gases, medido.** Cuando varios gases absorben en la misma banda, hay que decidir cómo se combinan sus tablas. Se midieron las opciones de la literatura frente al cálculo línea a línea (`herramientas/medir_solapamiento.py`):
- **Datos de la medida:** líneas reales de HITRAN2016 de H₂O, CO₂ y N₂O (las que trae el repositorio público PyRADS, Koll y Cronin; solo los isotopólogos principales, de 0 a 5000 cm⁻¹; solo para medir, no se distribuyen) y MT_CKD 4.3; en el PC de Carlos se repetirá con HITRAN2020 completo (con O₃ y CH₄).
- **Columnas:** tres columnas de prueba con los niveles de M3N y la composición sorteada de P3N (`fase52_referencia.py`): «ref» (suelo a 288 K, 21 kg/m² de agua precipitable), «tropical» (300 K) y «polar» (255 K, casi seca). Los agrupamientos fijos se hacen siempre con «ref»: las otras dos miden si sirven fuera de ella.
- **Métodos comparados**, todos con 16 puntos g (T6):
  - **PM:** la mezcla ordenada capa a capa (k correlacionado «premezclado»); con **+fr**, con la fracción de Planck de cada punto g;
  - **RORR:** cada gas con su tabla y solapamiento aleatorio con reordenación y reagrupamiento (Lacis y Oinas 1991 ⚠️; formulación de Amundsen et al. 2017, *A&A* 598, A97 ✅);
  - **SUMA:** cada gas con su tabla, sumadas como si estuvieran correlacionadas (lo ingenuo);
  - **MAPA:** cada número de onda se asigna una vez a un grupo según lo hondo que llega su luz en la columna «ref» (el nivel donde τ desde arriba llega a 1), y cada gas se promedia en cada grupo conservando su transmisión difusa media pesada con Planck: es el método de ecCKD (Hogan y Matricardi 2022, *JAMES* 14, e2022MS003033 ✅t; promedio «transmission» leído en su código, Apache 2.0 ✅). La fracción de Planck de cada grupo es exacta.
- **Resultado** (bandas 0 a 7 y 13, las que tienen solapamiento en estos datos; suma de los errores absolutos de cada banda, W/m²; peor error de calentamiento de una banda, K/día):

| Método | «ref»: al espacio · al suelo | «tropical» | «polar» | peor calentamiento |
|---|---|---|---|---|
| SUMA | 4,78 · 3,83 | 5,04 · 1,67 | 2,50 · 6,58 | 0,18 |
| RORR | 1,90 · 1,47 | 2,08 · 1,19 | 1,13 · 1,71 | 0,18 |
| PM | 1,68 · 0,84 | 1,75 · 0,58 | 1,09 · 1,01 | 0,18 |
| PM+fr | 0,81 · 0,30 | 0,93 · 0,15 | 0,27 · 0,61 | 0,03 |
| **MAPA** | **0,35 · 0,21** | **0,33 · 0,32** | 0,61 · 1,10 | **0,04** |

- **Lectura:**
  - Sumar las tablas sin más se equivoca hasta 5 W/m²: el solapamiento hay que tratarlo.
  - La fracción de Planck de cada punto g importa tanto como el solapamiento: sin ella, en la banda ancha de 10–350 cm⁻¹ el flujo al espacio se equivoca en −0,75 W/m² (con ella, −0,04 con MAPA).
  - MAPA y PM+fr son los mejores, del orden de una décima de W/m² por banda. MAPA es mejor en «ref» y «tropical»; en la columna «polar», muy seca, su error crece: el peor caso es la banda de 500–630 cm⁻¹, +0,54 W/m² al suelo, frente a −0,33 de PM+fr (con 32 grupos, +0,48: el número de grupos no lo arregla; es el agrupamiento hecho con una columna húmeda).
- **Decisión T9 (provisional): MAPA.** Con errores del mismo orden que la mejor alternativa:
  - **sirve para cualquier composición**, porque los gases se suman: el CO₂ se puede cambiar por worldbuilding sin rehacer las tablas (con PM habría que rehacerlas: una noche de cálculo cada vez) y el ozono puede variar de un sitio a otro;
  - **cuesta lo mismo que un solo gas** en el modelo (RORR obliga a ordenar en cada capa y en cada llamada: Amundsen et al. lo descartan para modelos 3D por lento);
  - la fracción de Planck de cada grupo es exacta y depende solo de la temperatura.
- **Decisión técnica T21 (provisional): el vapor, con su cantidad.** Para el vapor, que cambia mucho de un sitio a otro, la cantidad con la que se promedia en cada grupo es la de su propia fracción molar en esa casilla de la tabla. Por eso la dimensión del vapor pasa a ser logarítmica: 1·10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³, 4·10⁻³, 1,5·10⁻² y 6·10⁻². Medido: promediar el vapor con la cantidad de la columna «ref» empeora la polar (en la banda de 500–630 cm⁻¹, de +0,54 a −1,08 W/m² al suelo); con su propia cantidad queda como MAPA.
- **T15 revisada: temperaturas cada 15 K** (150–330 K, 13 en lugar de 7), interpolando ln k linealmente en T. Medido con `medir_solapamiento.py --interpolar` en la banda de 820–980 cm⁻¹ (la del continuo del vapor, la más sensible a la temperatura), error al suelo por interpolar en la rejilla:

| Rejilla e interpolación | «ref» | «tropical» | «polar» |
|---|---|---|---|
| cada 30 K, lineal en k | 0,62 | 0,82 | 0,02 |
| cada 30 K, ln k lineal en T | 0,28 | 0,29 | 0,00 |
| cada 30 K, ln k lineal en 1/T | 0,19 | 0,17 | 0,01 |
| **cada 15 K, ln k lineal en T** | **0,05** | **0,13** | **0,00** |

  (W/m²; sin interpolar, el mapa solo: 0,01, 0,05 y 0,00). El coste del generador crece en proporción (13/7).
- **Limitaciones declaradas:**
  - con 16 grupos y un solo agrupamiento, el error crece en columnas muy secas (≲ 0,5 W/m² por banda en la peor medida);
  - si hiciera falta, la mejora conocida es la partición en «hipercubo» de ecCKD (cada gas con sus propios grupos) o ajustar las tablas contra el línea a línea, como hace ecCKD.
- **Las tablas de la luz** se harán con el mismo método cuando esté el espectro de S3N (la ponderación es la de la estrella), con el agrupamiento por la profundidad a la que τ llega a 0,25, como ecCKD.

**El generador, adaptado (3.23.0):** `herramientas/generar_tablas_k.py` hace ahora las tablas del infrarrojo con MAPA, en dos fases:
- los mapas de cada banda (`mapas/mapa_bNN.npz`, con la columna «ref»; como mucho 3 a la vez por memoria);
- las casillas: cada gas, presión, temperatura (y fracción de vapor), promediado en los grupos.
- Salida: `tabla_k_<GAS>.npz` y `mapas_infrarrojo.npz` (peso de cada grupo y su fracción de Planck de 100 a 350 K).
- Las tablas de la 3.22.0 (una distribución k por casilla) no se llegaron a calcular: no se pierde nada.
- `fase52_infrarrojo.TablasInfrarrojo` lee las tablas y calcula los flujos de una columna: interpolación en ln p, T y ln de la fracción de vapor, de ln k.
- **Probado de punta a punta** (`test_generar_tablas_k.py`, con líneas sintéticas): tablas del generador → flujos, frente al línea a línea: −0,036 W/m² al espacio y −0,041 al suelo con 288 K (en bandas de ~10 y ~29 W/m²); −0,018 y −0,033 con 255 K.

**La luz: Rayleigh, datos y ozono (3.23.0, 13/10/2026, tarde):**
- **Decisión técnica T22 (provisional): la dispersión de Rayleigh, gas a gas** (`fase52_rayleigh.py`). La sección eficaz de cada gas sale de su índice de refracción medido y su factor de King, y la del aire es la suma pesada por la composición. Así sirve para el aire de P3N, con un 18,9 % de O₂.
  - Fuentes: N₂, Peck y Khanna 1966; O₂, Zhang et al. 2008 y, por debajo de 0,288 µm, la tabla de Smith et al. 1976 (la fórmula de Zhang extrapolada se queda un 3,5 % corta a 0,2 µm); Ar, Peck y Fisher 1964; CO₂, Bideau-Mehu et al. 1973 (de la base de refractiveindex.info, dominio público). Factores de King de Bates 1984 (como Bodhaine et al. 1999).
  - Comprobado: la mezcla de la Tierra reproduce el índice del aire de Peck y Reeder 1972 dentro del 0,1 % de 0,25 a 1,5 µm, y la sección a 550 nm, 4,5112·10⁻²⁷ cm², la de Bodhaine et al. 1999 (4,5105·10⁻²⁷) dentro del 0,02 %.
  - Para P3N, la sección por molécula es un 0,1–0,3 % mayor que en la Tierra (más N₂ en lugar de O₂).
  - El argón de P3N se supone como el de la Tierra (0,934 %): su efecto es despreciable aquí, pero es una suposición ⚠️.
- **Los datos de la luz** (`herramientas/descargar_luz.py`, para el PC de Carlos, ~50 MB) y sus lectores (`fase52_luz.py`):
  - el espectro de la estrella: modelos PHOENIX (Husser et al. 2013) de 5300, 5400 y 5500 K, log g = 4,5 y [Fe/H] = −0,5 (como τ Ceti; S3N tiene ~5420 K), interpolados en la temperatura; con el reparto de la luz entre las bandas. Por debajo de 0,2 µm (fuera de las bandas de M3N) no hace falta el ultravioleta lejano medido;
  - el ozono en el ultravioleta y el visible: Serdyuchenko et al. 2014, 11 temperaturas de 193 a 293 K, de 213 a 1100 nm (de 200 a 213 nm faltan datos: limitación ⚠️);
  - la absorción inducida por colisiones de HITRAN2024: N₂–N₂, O₂–O₂ y O₂–N₂.
  - Lector propio de FITS (sin astropy). Si dos bloques CIA cubren el mismo tramo, el lector se para en lugar de sumarlos a ciegas.
- **El ozono de P3N, propuesta 🔶** (`fase52_ozono.py`; es la decisión 5 del §5, de Carlos):
  - la forma del perfil de la Tierra (atmósfera estándar de EE. UU. 1976, Anderson et al. 1986: máximo de 7,8 ppm hacia 5,7 hPa);
  - escalada a ~180 DU de columna: los ~300 DU de media de la Tierra por el cociente 0,60 entre una estrella G8 V y el Sol de Rugheimer et al. 2013.
  - Limitaciones: sin latitud ni estación; la forma no cambia con la estrella; el efecto del 10 % menos de O₂ está por buscar.

**Las nubes: lo construido y lo que falta decidir (3.23.0, 13/10/2026, tarde).**

*Construido (sin conectar):*
- **Solapamiento entre capas y McICA, decisión técnica T23 (provisional)** (`fase52_nubes_solape.py`):
  - exponencial-aleatorio (Hogan e Illingworth 2000): dos capas nubladas seguidas se solapan con α = e^(−Δz/L), y las separadas por cielo despejado, al azar;
  - L según la latitud, la forma de Shonk et al. 2010 que usa el IFS: L = 0,75 + 2,149 cos²(latitud) km (2,9 km en el ecuador, 0,75 en los polos). Es una calibración terrestre (c), declarada;
  - McICA (Pincus et al. 2003): cada punto g de cada banda ve una subcolumna sorteada con esas reglas (el algoritmo del generador de ecRad, Apache 2.0, leído el 13/10). En media es exacto, con ruido sin sesgo.
  - Sin heterogeneidad del agua dentro de la nube (limitación ⚠️; ecRad la añade).
  - Probado (8 pruebas): los límites aleatorio y máximo; la cobertura entre ellos y decreciente con α; con 40 000 subcolumnas, la frecuencia de nube de cada capa, de cada par y de la columna coincide con la teórica dentro del error del muestreo.
- **La óptica de las gotas en las 16 bandas del infrarrojo** (`datos/optica_gotas_infrarrojo.npz`, de `herramientas/generar_optica_gotas.py`; `fase52_optica_nubes.TablaGotas`): β, ω y g para radios efectivos de 2 a 40 µm y temperaturas de 240 a 298 K (T17, T18). La de la luz espera al espectro de S3N.
- **El radio efectivo de las gotas, decisión técnica T24 (provisional):** Martin et al. 1994, r_e = (3 LWC / (4π ρ_w N k))^(1/3), con las constantes del IFS (ecRad, leído el 13/10):
  - océano: 50 núcleos de condensación por cm³ → ~51 gotas/cm³, k = 0,77;
  - tierra: 900 → ~313 gotas/cm³, k = 0,69;
  - la corrección por la lluvia de Wood 2000, y r_e entre 4 y 30 µm.
  - Con 0,3 g/m³ de agua: 12,3 µm sobre el océano y 6,9 µm sobre tierra (lo típico de la Tierra).
  - ⚠️ Son números de la Tierra actual. P3N tiene aire natural preindustrial (ficha §6.4): sobre tierra seguramente hay menos gotas (y más grandes). La cifra preindustrial está por buscar.

*Lo que falta y decide Carlos (la decisión 7 del §5, «ambición de las nubes»):*
1. **El agua de nube**, propuesta 🔶: dos trazadores nuevos del núcleo (agua líquida y hielo de nube), transportados como el vapor (el transporte positivo y conservativo de la 6.3 ya existe). Se forman por condensación de gran escala y se convierten en lluvia o nieve con una autoconversión (tipo Sundqvist et al. 1989 ⚠️). Hoy todo lo que condensa llueve al instante.
2. **La fracción de nube**:
   - (a) diagnóstica por la humedad relativa, con un umbral calibrado en el modo Tierra (Sundqvist et al. 1989 ⚠️), propuesta 🔶: es la más sencilla y la que mejor casa con (1);
   - (b) pronóstica (Tiedtke 1993 ⚠️): más física, bastante más trabajo y más parámetros.
3. **Las nubes de la convección**, propuesta 🔶: diagnosticadas a partir de la lluvia convectiva del esquema de Betts y Miller (Slingo 1987 ⚠️), porque ese esquema no produce agua de nube.
4. **La óptica del hielo** (no es esférico, así que Mie no sirve), propuesta 🔶: usar como dato una base de propiedades de cristales medidos o calculados (por ejemplo, la de Yang et al. 2013, la de MODIS ⚠️; por comprobar que se pueda descargar), en lugar de copiar una parametrización (Fu 1996 ⚠️). Es el mismo criterio que con HITRAN.
- **Riesgo declarado:** el umbral de humedad y la autoconversión son calibraciones (c). Se fijan en el modo Tierra antes de mirar P3N, como siempre.

**Lo que queda pendiente de la 5.2, en orden:**
- **El infrarrojo por bandas** con las tablas k y el dos flujos que ya tiene M3N (T5).
- **Los datos:**
  - las líneas de HITRAN, que se descargan en el PC de Carlos;
  - el continuo MT_CKD;
  - las propiedades ópticas de las gotas y el hielo de las nubes, que tienen sus propias fuentes, por decidir.
- **Las nubes** (§3.4) y el **acoplamiento** (§3.5).
- **La validación del modo Tierra** (§4), con su criterio fijado antes.

**Lo que necesita a Carlos:**
- La elección de (B) y del §2.
- Más adelante, descargar en su PC las líneas de HITRAN con una cuenta gratuita, porque desde el entorno de la IA no se llega a hitran.org (comprobado: sin respuesta).

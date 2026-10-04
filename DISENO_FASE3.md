# M3N — Fase 3: hielo marino

**Estado:** implementado en `fase2b_atmosfera.py` (interruptor I9, `hielo_marino`), sobre la v2.2c, y validado (sección 6).
**Fecha:** 02/10/2026 · **Versión prevista:** `v2.3`

---

## 0. Resumen

La Fase 3 de la tabla de fases es la *"retroalimentación hielo-albedo + calor latente"*, con un criterio de validación: *"sin bucle de enfriamiento/calentamiento descontrolado"*.

Se implementa como **hielo marino termodinámico con espesor real**, no como un umbral de temperatura:

- **El océano guarda energía (entalpía).** Si es positiva, es agua líquida a una temperatura. Si es negativa, es hielo de un espesor concreto. Congelar exige quitar calor latente y fundir, aportarlo, así que la energía se conserva por construcción.
- **El hielo aísla:** el calor del océano sube a través de él por conducción. La temperatura de su superficie sale del balance con la atmósfera de la Fase 2b.
- **El hielo refleja más luz** (albedo 0,65 frente a ~0,05–0,2 del agua). Esa es la retroalimentación hielo-albedo.
- **Nieve en tierra:** solo como diagnóstico. Sin precipitación (Fase 5) no puede haber nieve real.

---

## 1. Lo que enseñó el intento anterior (conversación perdida, 22-27/09/2026)

| Problema de entonces | Causa | Cómo se evita aquí |
|---|---|---|
| Un umbral de 0 °C en tierra congeló el continente entero, ecuador incluido | Albedo en escalón por temperatura más difusión: la inestabilidad del casquete pequeño (North 1981) | El hielo es un **espesor** que crece y se funde con energía real, y el albedo cambia gradualmente con el espesor. En tierra no hay umbral |
| −225 °C en la tierra polar sin difusión | No había atmósfera que diera calor de noche | Resuelto por la Fase 2b |
| La extensión del hielo dependía de la longitud | Se decidía con una foto instantánea | El espesor integra la energía en el tiempo |
| 8,9 % de desequilibrio de energía | Calor latente mal contabilizado | Entalpía: el calor latente es energía almacenada y se conserva por construcción |

---

## 2. Física

### 2.1 Entalpía del océano (Wagner y Eisenman 2015)

Para cada celda de océano, con C₀ la capacidad de la capa de mezcla (~48 m de agua):

```
E = C0 · (T0 − Tf) − ρi · Lf · h

E ≥ 0  →  agua líquida:   T0 = Tf + E / C0,   h = 0
E < 0  →  hielo:          h = −E / (ρi · Lf), T0 = Tf (agua bajo el hielo a punto de congelarse)
```

En cada paso, todos los flujos que llegan a la celda (atmósfera, transporte oceánico y océano profundo) cambian E, y E decide el estado.

### 2.2 Superficie del hielo: modelo de capa cero (Semtner 1976)

La temperatura de la superficie del hielo, T_s, cumple:

```
luz absorbida + infrarrojo que baja − σ·Ts⁴ − H(Ts) + (k_i / h)·(Tf − Ts) = 0
```

- Se resuelve por el método de Newton. H es el calor sensible con la estabilidad de Louis, Tiedtke y Geleyn (1982), como en la Fase 2b.
- **Máximo de 0 °C:** si con T_s = 0 °C sigue sobrando energía, el sobrante funde el hielo por arriba (CAM6/CICE: *"when F_TOP(Ts=0) ≥ 0, the surface melts at Ts = 0"*).
- La energía neta que entra por arriba (que, si no hay fusión, es igual a la conducción) cambia la entalpía. Si es negativa, el hielo crece por abajo; si es positiva, se funde.
- El calor sensible del hielo calienta o enfría la capa límite, igual que en el resto de superficies.

### 2.3 Calor del océano profundo

M3N solo tiene la capa de mezcla del océano. En la realidad, el agua profunda aporta calor bajo el hielo. Se representa con un flujo fijo de **4 W/m²** en la base del hielo (Wagner y Eisenman 2015; Maykut y Untersteiner 1971 usan ~2 W/m² en el Ártico central). **Se retira la misma cantidad total de todo el océano**, repartida por igual por área (hielo incluido): el calor se redistribuye, no se crea.

**Corrección v2.4.1:** hasta la v2.4 se retiraba solo del océano **libre**. En el caso límite de un océano casi entero helado, las pocas celdas libres recibían retiradas de miles de W/m², y con el océano entero helado no se retiraba nada, así que se creaba energía. Repartirlo por todo el océano conserva la energía siempre (`test_modelo.py` lo comprueba) y no cambia el orden de magnitud del flujo neto bajo el hielo, que es 4 W/m² × (1 − fracción de océano helado).

### 2.4 Albedo del hielo

- **Hielo desnudo: 0,65** (SHEBA, Perovich et al. 2002; open water 0,07 en el mismo estudio).
- Sin nieve encima hasta la Fase 5 (decisión de Carlos). Con nieve sería 0,85.
- **Hielo fino más oscuro:** el albedo pasa linealmente del del agua al del hielo entre 0 y 0,5 m de espesor. Que el hielo de más de 0,5 m tiene el albedo "pleno" sigue a CAM6; la forma lineal de la transición es una suposición propia y queda documentada como tal.
- Sin charcos de fusión (albedo 0,15–0,40). Es una limitación, menor en P3N, donde el hielo apenas recibe sol.

### 2.5 Nieve en tierra: solo diagnóstico

Decisión de Carlos: nada de nieve en tierra hasta la Fase 5, porque sin precipitación no puede haberla. El modelo devuelve `nieve_permanente_posible`: celdas de tierra donde, **si nevara**, la nieve no se fundiría nunca, porque el "mes" más cálido tiene media < 0 °C (clima de casquete glaciar EF de Köppen; año de P3N dividido en 12 meses de 22,5 días). No afecta a la física.

---

## 3. Parámetros

| Parámetro | Valor | Fuente |
|---|---|---|
| Temperatura de congelación del agua de mar | −1,8 °C | T_f = −0,054 °C/psu × S, con S ≈ 34 psu (CAM6 usa −1,8 °C) |
| Conductividad del hielo | 2,034 W/m/K | Hielo puro (Maykut y Untersteiner 1971; CAM6/CICE) |
| Densidad del hielo | 917 kg/m³ | — |
| Calor latente de fusión | 3,34×10⁵ J/kg | CAM6 |
| ρ_i·L_f | 3,06×10⁸ J/m³ | Wagner y Eisenman (2015): 9,5 W·año/m³ = 3,0×10⁸ |
| Albedo del hielo desnudo | 0,65 | Perovich et al. (2002), SHEBA |
| Espesor de albedo pleno | 0,5 m | CAM6 (forma de la transición: suposición) |
| Flujo del océano profundo | 4 W/m² (compensado) | Wagner y Eisenman (2015) |
| Espesor mínimo para la conducción | 0,01 m | Numérico (evita dividir por cero) |

---

## 4. Método numérico

- La entalpía viaja en la capa 0 de la columna del océano. Bajo el hielo esa capa está en T_f; los flujos de cada paso la desplazan, y al final del paso se reparte en (T₀, h).
- La difusión oceánica (D_oc) actúa sobre esa capa: el transporte de calor hacia los polos llega por debajo del hielo y lo funde o frena su crecimiento.
- La superficie del hielo se resuelve **antes** del infrarrojo de la atmósfera, así que lo que emite el hielo y lo que absorbe el aire usan la misma temperatura. Conserva la energía.
- **Convergencia:** con hielo activado, el equilibrio se juzga por la **temperatura real de la superficie** (la del hielo donde lo hay) y, en el agua libre y el hielo **fino (< 1 m)**, por el cambio de entalpía expresado en grados de la capa de mezcla (1 cm de hielo equivale a 0,015 °C). Así la **extensión** del hielo también tiene que estabilizarse. El **espesor** del hielo grueso no entra en el criterio; ver la sección 4.1.

### 4.1 Por qué el espesor no entra en el criterio de convergencia

La primera versión exigía también que la entalpía dejara de cambiar (1 cm de hielo al año como tolerancia). Con ese criterio, ninguna simulación convergía en 60 años: el hielo polar seguía engordando y llegaba a 12–18 m.

**No es un error, es la física de un polo sin verano.** En la Tierra, el hielo perenne del Ártico se estabiliza en 2–4 m porque en verano se funde por arriba. P3N tiene 1,7° de inclinación del eje, así que la superficie del hielo polar está todo el año entre −38 y −43 °C y nunca se funde. Sin fusión, el hielo crece hasta que el calor que deja escapar por conducción iguala al que le llega desde abajo:

```
h_eq = k_i · (Tf − Ts) / F_base
```

Es la relación que usan McKay (2000) y Warren et al. (2002) para el hielo de la "Tierra bola de nieve". Con k_i = 2,034, F_base = 4 W/m² y las temperaturas de superficie que da el modelo, h_eq va de 6 a 21 m (mediana de 18,5 m). En la simulación de 60 años, el espesor mediano era de 15,3 m y se acercaba a ese valor, así que el modelo se comporta como predice la teoría.

El tiempo para alcanzarlo es del orden de ρ_i·L_f·h_eq² / (k_i·ΔT), unos 50–60 años terrestres, y la aproximación es asintótica. Esperar a que el espesor se estabilice costaría más de un siglo simulado. Sin embargo, por encima de ~1 m el espesor apenas afecta al clima: pasar de 15 a 18,5 m cambia la conducción en ~0,9 W/m², y eso mueve la temperatura de la superficie unos 0,15–0,3 K (la emisión del hielo a 235 K responde con ~2,9 W/m²/K, más el calor sensible). Por eso el criterio usa la temperatura y la extensión, y el espesor final se documenta como **"en camino al equilibrio"**.

**Limitación conocida:** en la realidad, un hielo tan grueso fluye como un glaciar ("sea glaciers", Goodman y Pierrehumbert 2003; Pollard y Kasting 2005) y se adelgaza o se extiende hacia aguas más cálidas. M3N no tiene dinámica del hielo, solo termodinámica.
- **Aceleración (1): salto geométrico** de la Fase 2b, aplicado a la entalpía (agua y hielo juntos).
- **Corrección v2.4.1:** el año en que actúa el salto del hielo grueso nunca se da por convergido, porque el estado acaba de cambiar.
- **Aceleración (2): salto del hielo grueso.** En los años 6, 12 y 18, para cada celda con más de 1 m de hielo se estima el calor que llega a su base con lo medido ese año (F_base = conducción media − ρ_i·L_f·dh/dt) y se lleva el espesor a su equilibrio h_eq = k_i·(T_f − T_s)/F_base, limitado a entre la mitad y el doble del espesor actual. No es física nueva: solo acorta el camino. Después se siguen simulando años y el criterio de convergencia es el mismo.
  - Sin este salto, P3N tardaba 42 años (33 minutos) en converger, con el hielo a 11,7 m de media y lejos aún de su equilibrio.
  - Con el salto, tarda 22 años (18 minutos), con 15,2 m de media. Clima global: 15,32 °C frente a 15,30 °C. Polos: −42,0 °C frente a −42,5 °C (con el hielo más cerca de su equilibrio, el polo es algo más frío, como predice la sección 4.1).

---

## 5. Prueba de biestabilidad

Con el albedo del hielo, un planeta puede tener **dos climas estables**: uno templado con casquetes polares y otro totalmente helado, la "Tierra bola de nieve" (North 1981; Wagner y Eisenman 2015). Para saber si P3N los tiene, se simula dos veces:

- arrancando **sin hielo**;
- arrancando con **todo el océano cubierto de 5 m de hielo** y el aire frío.

Si las dos acaban en el mismo estado, P3N tiene un solo clima posible con estos parámetros. Si no, la biestabilidad es física, no un fallo, y es información para el worldbuilding: P3N podría haber estado, o estar, en cualquiera de los dos estados según su historia.

---

## 6. Validación

| Prueba | Criterio | Resultado |
|---|---|---|
| V0. Todo apagado = v2.2 | Idéntico bit a bit | ✅ Diferencia 0,0 (5 años, como la referencia) |
| Hielo apagado = v2.2c | Idéntico | ✅ Diferencia 0,0 en el aire a 2 m, el suelo y la troposfera (10 años, como la referencia) |
| Energía (incluido el calor latente) | Diferencia en lo alto de la atmósfera < 10⁻³ | ✅ 1,7×10⁻⁵ (P3N); 1,4×10⁻⁵ (modo Tierra) |
| Biestabilidad | Mismo estado desde "sin hielo" y desde "todo helado"; si no, documentarlo | ✅ Mismo estado (60 años cada una): 15,31 °C, 8,7–9,1 % del océano helado, perfiles zonales idénticos. **P3N es monoestable** |
| Espesores | Del orden del hielo perenne real (2–4 m en el Ártico); si no, explicar por qué | ⚠️ 12–18 m y creciendo hacia el equilibrio teórico (mediana de 18,5 m). Explicado en la sección 4.1: polo sin verano |
| Modo Tierra | Extensión de hielo marino del orden de la observada | ⚠️ Máximo ártico bien; deshielo de verano y hielo antártico insuficientes. Ver la sección 6.1 |

### 6.0 Clima de P3N con hielo y masa de S3N

Con el hielo activado y la masa de S3N de la v2.2c (0,874 masas solares), el aire a 2 m medio global baja de 15,96 a **15,30 °C** (convergencia en 22 años, 18 minutos). El 8,7–9,1 % del océano está helado según la estación, con −42/−43 °C en los polos. Como 15,30 °C está dentro del rango sorteado (15–16 °C), **la masa no se recalibra** (decisión de Carlos, 02/10/2026). Como referencia de sensibilidad: 0,877 → 16,02 °C y 0,878 → 16,25 °C (unos 0,24 °C por cada milésima de masa solar).

**Con las correcciones de la v2.4.1** (compensación del océano profundo repartida por todo el océano; ningún año con salto del hielo cuenta como convergido), mismo mapa `prueba1` (03/10/2026): **15,29 °C** (−0,004 °C respecto a la v2.4), océano helado 8,75–9,1 % (igual), polos −42,7 / −43,2 °C, ecuador 29,8 °C. Converge en **25 años** en vez de 22: los saltos del hielo grueso llegan en los años 6, 12 y 17, y el estado necesita unos años más para asentarse después del último. Desequilibrio de energía final: 1,2×10⁻⁵ W/m². Los resultados prácticamente no cambian; lo que cambia es que ahora el equilibrio está garantizado.

### 6.1 Modo Tierra con hielo

Misma configuración que la calibración de D en la v2.2c (órbita, máscara de tierra y gravedad terrestres; D_atm = 2,4 y D_oc = 0,12), con el hielo activado. Converge en 40 años, con 13,8 °C de media global (14,5 °C sin hielo). Extensión del hielo marino, contando como helada toda la celda (la rejilla de 5° no permite concentraciones parciales):

| | M3N | Observado (NSIDC, media 1981–2010) |
|---|---|---|
| Ártico, máximo | 14,2 millones de km² | 15,6 (marzo) |
| Ártico, mínimo | 12,0 | ~6 (septiembre) |
| Antártico, máximo | 4,7 | ~18–19 (septiembre) |
| Antártico, mínimo | 1,6 | ~3 (febrero) |
| Espesor ártico medio | 9 m | 2–3 m |

**Lectura:**

- **El máximo ártico es del orden correcto.** La retroalimentación hielo-albedo y el calor latente funcionan.
- **El Ártico no pierde bastante hielo en verano y es demasiado grueso.** La causa principal no está en el hielo: el aire del polo norte ya está a −34 °C de media anual **sin hielo** (observado: unos −18 °C). Faltan las nubes y el vapor de agua, que en el invierno ártico devuelven mucho infrarrojo hacia la superficie (Fase 5). También faltan los charcos de fusión (albedo de verano ~0,5 en vez de 0,65) y las grietas de agua abierta entre el hielo.
- **El hielo antártico es mucho menor que el real.** Se explica sobre todo por lo que el modelo no tiene: en la realidad, los vientos empujan el hielo antártico hacia el norte (dinámica del hielo) y el océano austral está estratificado por el agua dulce. Un modelo solo termodinámico con una capa de mezcla de 48 m que recibe el transporte de calor hacia los polos no puede reproducirlo.
- **Efecto sobre la Antártida:** con hielo marino, el polo sur pasa de −14,7 a −33,6 °C. El océano austral helado ya no calienta el aire que llega al continente. Esto acerca el modelo al valor observado (unos −49 °C en el polo sur, a 2835 m).

**Conclusión:** el esquema de hielo se comporta bien en lo que depende de él. Los sesgos del modo Tierra vienen de fases todavía no implementadas (nubes y humedad en la Fase 5, dinámica). No se ajusta ningún parámetro del hielo para compensarlos, porque eso escondería el error en el sitio equivocado.

---

## Referencias

- Maykut, G. A. y Untersteiner, N. (1971). Some results from a time-dependent thermodynamic model of sea ice. *J. Geophys. Res.* 76.
- NSIDC, *Sea Ice Today* (medias de 1981–2010 de la extensión del hielo ártico y antártico).
- North, G. R. (1981). Energy balance climate models. *Rev. Geophys. Space Phys.* 19.
- Perovich, D. K. et al. (2002). Seasonal evolution of the albedo of multiyear Arctic sea ice. *J. Geophys. Res.* 107(C10).
- Semtner, A. J. (1976). A model for the thermodynamic growth of sea ice in numerical investigations of climate. *J. Phys. Oceanogr.* 6, 379–389.
- McKay, C. P. (2000). Thickness of tropical ice and photosynthesis on a snowball Earth. *Geophys. Res. Lett.* 27(14), 2153–2156.
- Pollard, D. y Kasting, J. F. (2005). Snowball Earth: A thin-ice solution with flowing sea glaciers. *J. Geophys. Res.* 110, C07010.
- Goodman, J. C. y Pierrehumbert, R. T. (2003). Glacial flow of floating marine ice in "Snowball Earth". *J. Geophys. Res.* 108(C10), 3308.
- Warren, S. G. et al. (2002). Snowball Earth: Ice thickness on the tropical ocean. *J. Geophys. Res.* 107(C10), 3167.
- Wagner, T. J. W. y Eisenman, I. (2015). How climate model complexity influences sea ice stability. *J. Climate* 28, 3998–4014.
- NCAR, *CAM6 Scientific Guide*, cap. 6: Sea Ice Thermodynamics.

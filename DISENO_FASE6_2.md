# Fase 6.2 de M3N: ecuaciones primitivas en seco (núcleo dinámico propio, 20 capas)

**Estado: DISEÑO v1 (05/10/2026), sin código todavía.** Parte de las decisiones de la 6.1 que confirmó Carlos: rejilla C, filtro polar desde 60°, leapfrog con RAW y semiimplícito obligatorio.

Etiquetas:
- ✅ verificado en una fuente;
- ✅² verificado en dos fuentes independientes;
- ⚠️ pendiente;
- 🔶 elección de la IA, provisional.

Fuentes de las ecuaciones discretas:
- **IFS Documentation Cy25r1, Parte III, cap. 2** (ECMWF): la discretización vertical de Simmons y Burridge (1981), ecs. 2.11–2.26. Leído a través de un resumen automático, que confundió un signo; por eso se coteja con la fuente siguiente.
- **Código del núcleo espectral del GFDL** (en la distribución Isca, licencia GPL), archivos `press_and_geopot.F90` y `spectral_dynamics.F90`, rutina `four_in_one`. Se usa **solo como fuente de las fórmulas** (principio de Carlos del 05/10): nada de ese código entra en M3N.

---

## 0. Qué cambia respecto a la 6.1

La 6.1 tenía una capa de fluido con altura h. La 6.2 tiene la atmósfera de M3N tal como es:
- **20 capas en σ = p/p_s**, con los mismos seminiveles que ya usa la v3.0. Hay 21 seminiveles, con σ = 0 arriba y σ = 1 en el suelo, y el tope está a p = 0.
- **Variables:** u y v en cada capa, T en cada capa, y la presión en superficie p_s (una por columna).
- **Altitud:** el geopotencial Φ ya no es una variable; se calcula a partir de T (hidrostática).

Seguimos en seco: sin vapor, ni radiación, ni suelo. El forzamiento es el idealizado de Held y Suarez. La física de M3N entra en la 6.3.

## 1. Ecuaciones continuas (Phillips 1957 ✅; forma vectorial invariante)

Con v = (u, v), ζ la vorticidad relativa, f el parámetro de Coriolis, K = |v|²/2 y σ̇ = dσ/dt:

1. **Momento:** ∂v/∂t = −(ζ + f) k × v − ∇(K + Φ) − R T ∇ln p_s − σ̇ ∂v/∂σ + F
2. **Termodinámica:** ∂T/∂t = −v·∇T − σ̇ ∂T/∂σ + κ T ω/p + Q/c_p
3. **Masa (presión en superficie):** ∂p_s/∂t = −∫₀¹ ∇·(p_s v) dσ
4. **Velocidad vertical:** p_s σ̇(σ) = −σ ∂p_s/∂t − ∫₀^σ ∇·(p_s v) dσ′
5. **Hidrostática:** ∂Φ/∂ln σ = −R T, con Φ(σ=1) = Φ_s (la orografía)
6. **ω = dp/dt** = p_s σ̇ + σ (∂p_s/∂t + v·∇p_s)

**Energía total que se conserva** (sin forzamiento): E = ∫∫ (K + c_p T) p_s dσ dA/g + ∫ Φ_s p_s dA/g.

## 2. Discretización vertical: Simmons y Burridge (1981) en σ pura

Notación: los seminiveles son σ_{k±1/2}, con k = 1 la capa de arriba; Δσ_k = σ_{k+1/2} − σ_{k−1/2}; Δp_k = Δσ_k p_s.

### 2.1 Hidrostática ✅²

- **Entre seminiveles:** Φ_{k−1/2} = Φ_{k+1/2} + R T_k ln(σ_{k+1/2}/σ_{k−1/2}), partiendo de Φ_{N+1/2} = Φ_s.
- **En el centro de la capa:** Φ_k = Φ_{k+1/2} + α_k R T_k.
- **El coeficiente α_k:**
  - α_k = 1 − (σ_{k−1/2}/Δσ_k) · ln(σ_{k+1/2}/σ_{k−1/2}) para k > 1;
  - **α_1 = ln 2** en la capa de arriba, porque σ_{1/2} = 0.

  Fuentes: IFS ec. 2.22–2.23; GFDL `press_and_geopot.F90`, líneas 176–179 ("alpha = 1 − p_half(k)·(ln p_half(k+1) − ln p_half(k))/(p_half(k+1) − p_half(k))").
- **Consecuencia útil (cálculo propio):** en σ pura, ln(σ_{k+1/2}/σ_{k−1/2}) y α_k son constantes en todo el globo. Por eso Φ_k − Φ_s depende solo de las temperaturas de la columna, no de p_s.

### 2.2 Fuerza de presión ✅²

- **La forma de partida** (IFS ec. 2.24; GFDL `four_in_one`, x1):
  (R T ∇ln p)_k = R T_k / Δp_k · [ln(p_{k+1/2}/p_{k−1/2}) ∇p_{k−1/2} + α_k ∇Δp_k]
- **En σ pura se simplifica** (cálculo propio, coincide con x1 = 1/p_s del código):
  [ln(σ_{k+1/2}/σ_{k−1/2}) σ_{k−1/2}/Δσ_k + α_k] ∇ln p_s = **1 · ∇ln p_s**

Es decir, la fuerza de presión en la capa k es −∇Φ_k − R T_k ∇ln p_s.

### 2.3 Masa y movimiento vertical ✅²

- **Presión en superficie:** ∂p_s/∂t = −Σ_k ∇·(Δσ_k p_s v_k).
- **Flujo de masa vertical hacia abajo a través de cada seminivel:** W_{k+1/2} = −Σ_{j≤k} ∇·(Δp_j v_j) − σ_{k+1/2} ∂p_s/∂t, con W = 0 arriba y abajo. Fuentes: GFDL, variable `wg`; IFS ec. 2.17–2.18.

### 2.4 Término de conversión de energía κTω/p ✅²

Es el GFDL `four_in_one` con A = 0 y B = σ:
- (ω/p)_k = v_k·∇ln p_s − (1/Δp_k) · [ln(σ_{k+1/2}/σ_{k−1/2}) · Σ_{j<k} ∇·(Δp_j v_j) + α_k ∇·(Δp_k v_k)]

  Comprobación término a término con el código ✅:
  - el código hace dt_T −= κT·x5, con x5 = x4 − v·∇ln p_s;
  - x4 = (dmean_tot · dlog_3 + dmean · dlog_1)/Δp;
  - dmean = ∇·(Δp v) en σ pura, y dmean_tot es su suma para las capas de encima;
  - dlog_3 = ln(σ_{k+1/2}/σ_{k−1/2});
  - dlog_1 = ln p_{k+1/2} − ln p_k = α_k.

  Así que (ω/p)_k = −x4 + v·∇ln p_s, que es exactamente la fórmula de arriba. En la capa 1 la suma de encima es 0, y α_1 = ln 2.
- Δp_j ∇·v_j + Δσ_j v_j·∇p_s = ∇·(Δp_j v_j).

### 2.5 Advección vertical ✅²

[X]_k = (1/(2Δp_k)) · [W_{k+1/2}(X_{k+1} − X_k) + W_{k−1/2}(X_k − X_{k−1})] (IFS ec. 2.19), para u, v y T.

**Propiedad que da Simmons y Burridge (1981, por su título ✅):** con estas formas se conservan la energía total y el momento angular, si las formas horizontales son coherentes con ellas. Esa coherencia horizontal es nuestra parte (§3).

## 3. Discretización horizontal: rejilla C (la de la 6.1, ampliada)

- **Dónde va cada variable:** T_k, p_s y Φ en los centros; u_k en las caras este y oeste; v_k en las caras norte y sur. Todo como en la 6.1.
- **Flujos de masa por cara:** U_k = u_k · \overline{Δp_k}^cara · L_cara, como F_u y F_v en la 6.1, con Δp_k = Δσ_k p_s y p_s promediado aritméticamente a la cara.
- **Momento:**
  - vorticidad potencial en las esquinas, q_k = (ζ_k + f)/\overline{Δp_k}^esquina;
  - el término q·flujo, en la forma de Sadourny que conserva la energía (la misma demostración de la 6.1, ahora capa a capa);
  - gradiente de la función de Bernoulli B_k = K_k + Φ_k (con K ponderado por áreas, como en la 6.1);
  - el término R T ∇ln p_s en la cara, con T promediado: −\overline{R T_k}^cara · δ ln p_s / dx. 🔶 Hay que demostrar que este promedio es el que cierra la energía junto con el término v·∇ln p_s de §2.4. La alternativa es la forma de Arakawa y Suarez (1983) ⚠️. **El árbitro será la prueba numérica de conservación exacta** (§6, prueba 1), igual que en la 6.1.
- **Temperatura:** en forma de flujo, ∂(Δp T)/∂t = −∇·(Δp v \overline{T}^cara) − (flujo vertical) + Δp · κT ω/p. Así la entalpía c_p T Δp solo cambia por la conversión de energía, que tiene que cancelarse con el trabajo de la fuerza de presión (§3, prueba 1).
- **Filtro polar:** el de la 6.1 (desde 60°), aplicado a las tendencias de u, v, T y p_s. No toca la media zonal, así que la masa sigue siendo exacta.

## 4. Avance en el tiempo: leapfrog semiimplícito con RAW

### 4.1 Idea

- Las ondas de gravedad (hasta ~350 m/s) obligarían a pasos de ~150 s (6.1, §5).
- El esquema semiimplícito separa la parte **lineal** que las produce, alrededor de un estado de referencia en reposo (T_ref(σ), p_ref), y la calcula como media entre los instantes n−1 y n+1 (Robert, Henderson y Turnbull 1972 ✅, la idea; Simmons, Hoskins y Burridge 1978 ⚠️, la estabilidad). El resto se calcula como siempre.
- **Estado de referencia 🔶:** isotermo, T_ref = 300 K, p_ref = P0, con peso β = 0,5 (centrado).
  - Es lo que usa el núcleo del GFDL: ✅ `spectral_dynamics.F90`, alpha_implicit = 0,5 (línea 192) y `ref_temperature_implicit(:) = 300.` (línea 476).
  - **Dato independiente sobre RAW:** el mismo código trae RAW desactivado por defecto, con el comentario "0.5 is the desired value, but this appears unstable. Requires further testing" (línea 203), y usa un coeficiente de Robert ν = 0,04 (línea 191). Esto cuadra con nuestro hallazgo de la 6.1 (§5): RAW pierde estabilidad. Nuestro análisis da además el remedio: con ν = 0,05 la amplificación a ω·Δt = 0,7 baja de 1,0056 a 1,0010. 🔶 Se mantiene ν = 0,2 con ω_s·Δt ≤ 0,45, y en la 6.2 se medirá si conviene bajar ν.
  - **Por qué tan caliente:** un T_ref más frío que la atmósfera real puede volver inestable el esquema (SHB78 ⚠️).

### 4.2 Las ecuaciones lineales

Sea P_k = (G·T)_k + R T_ref ln p_s, el "geopotencial lineal". G es la matriz de la hidrostática de §2.1.
- **Para u y v:** δ_t v = −∇\overline{P}^t + (resto explícito).
- **Para T y ln p_s:** δ_t (T, ln p_s) = −(τ, ν)·\overline{D}^t + (resto), con D = ∇·v en cada capa. τ sale de linealizar §2.4 alrededor de T_ref; ν = Δσ.
- **Combinadas:** δ_t P = −M·\overline{D}^t + (resto), con M = G τ + R T_ref ν^T, una matriz de 20 × 20.
- **Eliminando v:** (I − (2βΔt)² M ∇²) P^{n+1} = (lado derecho conocido).

### 4.3 Cómo se resuelve: exacto y barato (diseño propio)

1. **Modos verticales:** M = E Λ E⁻¹. Los valores propios λ_m = c_m² son los cuadrados de las velocidades de las ondas de gravedad de cada modo vertical (20 modos). Se calcula una sola vez.
2. **En la horizontal**, ∇² = div ∘ S ∘ grad en la rejilla C, con S el filtro polar. La rejilla es uniforme en longitud y tanto S como los operadores actúan fila a fila, así que la **FFT en longitud separa cada número de onda zonal k**.
3. **Para cada (modo m, onda k)** queda un sistema **tridiagonal en latitud** de 36 incógnitas. Se resuelve de forma exacta, sin iteraciones: son 20 × 37 sistemas pequeños por paso.

Coste estimado: despreciable frente al resto ⚠️, a medir.

### 4.4 Estabilidad con RAW (cálculo propio, oscilador modelo)

Factor de amplificación por paso, con la parte rápida (ω_f) implícita centrada y la lenta (ω_s) explícita:

| ω_s·Δt \ ω_f·Δt | 0 | 0,5 | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|---|
| 0 | 1,0000 | 0,9987 | 0,9935 | 0,9808 | 0,9631 | 0,9548 |
| 0,3 | 0,9998 | 0,9973 | 0,9902 | 0,9766 | 0,9604 | 0,9533 |
| 0,45 | 1,0001 | 0,9972 | 0,9888 | 0,9745 | 0,9590 | 0,9525 |
| 0,5 | 1,0004 | 0,9974 | 0,9884 | 0,9739 | 0,9585 | 0,9523 |

- **Las ondas de gravedad quedan estables con cualquier paso** y ligeramente amortiguadas (≤5 %/paso).
- **El límite pasa a la parte explícita:** ω_s·Δt ≲ 0,45, la misma cifra que en la 6.1.
- **Estimación del paso** (ω_s ≈ U/Δx_ef + f, con Δx_ef la anchura de la celda a 60°):

| | Δx_ef | ω_s (U = 100 m/s) | Paso máximo |
|---|---|---|---|
| Tierra | 278 km | 3,6·10⁻⁴ + 1,5·10⁻⁴ s⁻¹ | ≈ 880 s |
| P3N | 204 km | 4,9·10⁻⁴ + 1,8·10⁻⁴ s⁻¹ | ≈ 670 s |

**Decisión 🔶:** la dinámica da **2 pasos por cada paso de la física**: 450 s en la Tierra y 496 s en P3N. Además, se vigila el número de Courant en cada paso y se avisa si se supera 0,35.

## 5. Forzamiento de Held y Suarez (1994) ✅²

Valores comprobados en el artículo (vía MITgcm y Thatcher y Jablonowski 2016) y en el código de referencia `hs_forcing.F90`:

| Constante | Valor |
|---|---|
| T_0 | 315 K |
| T_estratosfera | 200 K |
| ΔT_y | 60 K |
| Δθ_z | 10 K |
| σ_b | 0,7 |
| k_f | 1/día |
| k_a | 1/40 días |
| k_s | 1/4 días |

- **Temperatura de equilibrio:** T_eq = max(200, [315 − 60 sen²φ − 10 ln(p/P0) cos²φ](p/P0)^κ).
- **Relajación de T:** con k_T = k_a + (k_s − k_a) max(0, (σ − σ_b)/(1 − σ_b)) cos⁴φ.
- **Rozamiento:** −k_v v, con k_v = k_f max(0, (σ − σ_b)/(1 − σ_b)).

En la Tierra se toman R = 287, c_p = 1004 y g = 9,80616 (HS94). La difusión que necesite el núcleo para estabilizarse se trata en §7.

## 6. Pruebas, en orden (cada una tiene que pasar antes de la siguiente)

1. **Conservación exacta del sistema semidiscreto** (sin filtro y con un paso infinitesimal):
   - masa a 10⁻¹⁵;
   - **energía total a 10⁻¹²**, con un estado aleatorio y montañas.

   Es el árbitro de las formas de §3.
2. **Reposo isotermo sobre montañas:** con T uniforme y p_s = P0·exp(−Φ_s/(R T)), la fuerza de presión tiene que ser **exactamente cero** (cálculo propio: en σ pura, ∇Φ_k = ∇Φ_s, que se cancela con R T ∇ln p_s). Es la prueba clásica del error de presión de las coordenadas σ, que aquí tiene que dar 0 salvo redondeo.
3. **Velocidades de las ondas de gravedad:** los valores propios de M frente a la teoría de una atmósfera isoterma. Se espera que el modo externo vaya a ≈ √(γ R T) ≈ 347 m/s a 300 K ⚠️ (onda de Lamb; hay que verificarlo en Vallis).
4. **Equilibrio de Jablonowski y Williamson (2006):** su estado inicial es una solución estacionaria (equilibrio de viento térmico) y tiene que mantenerse 30 días. Fórmulas ✅, cotejadas con el código del modelo MPAS (`mpas_init_atm_cases.F`, líneas 462–468, 671–690 y 850–866), que coinciden término a término con lo que recordaba de JW06:
   - **Coordenada:** η = p/P0 (con p_s = 1000 hPa uniforme, η = σ); η_v = (η − 0,252)·π/2.
   - **Viento:** u = 35 · cos^{3/2}(η_v) · sen²(2φ); v = 0.
   - **Temperatura de fondo:** T̄ = 288 · η^{RΓ/g} para η ≥ 0,2, con Γ = 0,005 K/m. Por encima (η < 0,2) se suma 4,8·10⁵ · (0,2 − η)⁵.
   - **Temperatura:** T = T̄ + ¾ (η π u₀/R) sen η_v cos^{1/2} η_v · {[−2 sen⁶φ (cos²φ + 1/3) + 10/63] · 2u₀ cos^{3/2} η_v + [8/5 cos³φ (sen²φ + 2/3) − π/4] · aΩ}.
   - **Geopotencial en superficie:** Φ_s = u₀ cos^{3/2}((1 − 0,252)π/2) · {[…] u₀ cos^{3/2}(…) + […] aΩ}, con los mismos corchetes.
   - **Perturbación para la prueba 5:** u′ = 1 m/s · exp(−(r/R)²), con R = 0,1 radios del planeta y centro en (20° E, 40° N).
5. **Onda baroclínica de JW06:** la misma, con una pequeña perturbación. A 5° solo se compara cualitativamente: la onda tiene que crecer hacia el día 7–9 ⚠️.
6. **Held y Suarez, 1200 días** (se descartan los 200 primeros):
   - medias zonales de u y T comparadas con las figuras de HS94: corrientes en chorro de ~30 m/s hacia 45° ✅ (la altura exacta ⚠️);
   - temperatura y momento angular estables.
7. **Coste y convergencia:** 5° frente a 2,5° en Held y Suarez, si el coste lo permite.

## 7. Difusión horizontal: el punto que más rigor exige

- **Por qué hace falta algo:** todos los núcleos reales disipan algo de energía en la escala de la rejilla. La cascada de enstrofía la acumula ahí y, sin disipación, el modelo se llena de ruido.
  - Held y Suarez la incluyen en su propuesta: "∇⁸ con un tiempo de amortiguamiento de 0,1 días en la escala más pequeña" para su modelo espectral ⚠️ (por verificar en el artículo).
  - El GISS usaba un filtro de Shapiro ✅.
- **Principio de Carlos (no esconder errores ajustando parámetros):**
  - La difusión no se usa para tapar inestabilidades. Primero se comprueba que el núcleo **sin difusión** conserva lo que debe (pruebas 1–4).
  - Solo entonces se añade la hiperdifusión mínima justificada por la escala de la rejilla, con su coeficiente expresado como un **tiempo de amortiguamiento de la onda más corta** (convención de HS94).
  - Se documenta cuánta energía disipa.
- **Lo que usan otros al correr Held y Suarez** (como dato de calibración, no como sustituto):
  - **CESM/CAM (NCAR) ✅:** hiperdifusión de cuarto orden (∇⁴) con un tiempo de amortiguamiento de **0,5 días en la escala más pequeña**. Dicen textualmente que se eligió por ser "una de las amortiguaciones más débiles que aún daba estructuras suaves en el campo de vorticidad" a T42 y T85.
  - **DCPAM (Japón) ✅:** "laplaciano de 4.º orden" (∇⁸, por la descripción ⚠️) con 0,1 días en el número de onda máximo, 1200 días de integración, y resoluciones de T21 a T170 con 20 capas.
- 🔶 **Propuesta:** ∇⁴ sobre u, v y T (sobre la parte rotacional y divergente del viento, y sobre T a lo largo de las superficies σ), con un tiempo de amortiguamiento de 0,5 días en la escala de 2Δx, el mismo criterio que CESM. Nuestra rejilla de 5° equivale más o menos a T21–T24, así que se comprobará si basta con algo más débil (criterio de CESM: la amortiguación más débil que deje suave la vorticidad). Se aplicará en la forma que respeta la conservación de la masa.

## 8. Estructura del código prevista

- **`fase6_nucleo.py`** (nuevo):
  - clase `NucleoSeco` (rejilla, operadores, tendencias, semiimplícito y paso);
  - reutiliza `Rejilla` y `_filtro_fourier` de `fase6_aguas_someras.py`.
- **`fase6_forzamientos.py`:** Held y Suarez, y el estado de JW06.
- **`test_fase6_2.py` y `validar_fase6_2.py`.**
- **Interruptor I16 `nucleo_dinamico`, desactivado.** No se toca nada de M3N hasta la 6.3.

## 9. Lo que queda por verificar antes de programar (orden de trabajo)

1. ~~Las fórmulas de JW06~~ ✅ (MPAS).
2. ~~T_ref = 300 K en el GFDL~~ ✅. Queda leer su fundamento en SHB78 ⚠️.
3. ~~La hiperdifusión que se usa con HS94~~ ✅ (CESM, DCPAM).
4. La onda de Lamb en un sistema hidrostático.
5. Las formas horizontales de §3 (se deciden con la prueba 1).

## Referencias

- Simmons, A. J., y Burridge, D. M. (1981). Mon. Wea. Rev. 109, 758–766 ✅ (cita); ecuaciones vía la IFS Documentation y el código del GFDL ✅².
- ECMWF, IFS Documentation Cy25r1, Part III: Dynamics and Numerical Procedures, cap. 2 ✅.
- Held, I. M., y Suarez, M. J. (1994). BAMS 75, 1825–1830 ✅; constantes ✅².
- Jablonowski, C., y Williamson, D. L. (2006). QJRMS 132, 2943–2975 ✅ (cita).
- Robert, A., Henderson, J., y Turnbull, C. (1972). Mon. Wea. Rev. 100, 329–335 ✅ (cita).
- Simmons, A. J., Hoskins, B. J., y Burridge, D. M. (1978). Mon. Wea. Rev. 106, 405–412 ⚠️.
- Arakawa, A., y Suarez, M. J. (1983). Mon. Wea. Rev. 111, 34–45 ⚠️.
- Williams, P. D. (2009). Mon. Wea. Rev. 137, 2538–2546 ✅.
- Phillips, N. A. (1957). J. Meteor. 14, 184–185 ✅.

---

## 10. Resultados (05/10/2026): código `fase6_nucleo.py`, `fase6_forzamientos.py`, `test_fase6_2.py` (6 pruebas)

### Prueba 1: conservación exacta ✅

- **Energía:** dE/dt sale 1,2·10⁻¹⁴ veces cada contribución por separado, con un estado aleatorio y una montaña, sin filtro.
- **Masa:** 10⁻¹⁷.
- **La demostración (§3 bis, abajo)** funcionó a la primera.

### §3 bis: por qué se conserva la energía (demostración propia)

1. **El trabajo de ∇Φ en las caras** (Σ_f F_f (−δΦ)) es, sumando por partes, Σ_i Φ_k D_k.
   - Desarrollando Φ_k con la hidrostática de SB81, la parte de Φ_s se cancela con Φ_s ∂p_s/∂t.
   - El resto, Σ_j R T_j dln_j Σ_{k<j} D_k + Σ α_k R T_k D_k, es exactamente lo contrario de la conversión de energía de §2.4 (parte 1).
2. **El término R T ∇ln p_s** realiza en cada cara un trabajo Y_f = A_f Δp_f u_f R \overline{T}^f δln p_s/dx.
   - En la ecuación de la entalpía se reparte mitad y mitad entre las dos celdas de la cara (parte 2 de la conversión).
   - Por construcción, ambas sumas coinciden, sea cual sea el promedio de T en la cara.
3. **La advección vertical de u y v** (forma de SB81) y el término K·(cambio de Δp por el flujo vertical) se cancelan, siempre que W en las caras sea la media aritmética de las celdas y K se pondere por áreas (como en la 6.1).
4. **Coriolis y vorticidad:** como en la 6.1, capa a capa.
5. **La entalpía** va en forma de flujo, así que el transporte solo la redistribuye.

### Prueba 2: reposo isotermo sobre montañas ✅

Con la montaña del caso 5 (2000 m), las tendencias son ≤ 3·10⁻¹⁶: **ningún viento falso**. Es el error clásico de las coordenadas σ, y aquí es exactamente cero.

### Prueba 3: ondas de gravedad ✅

- **Modo externo:** 344,6 m/s frente a √(γRT) = 347,2 m/s a 300 K, la onda de Lamb: −0,75 % con 20 capas.
- **Modos internos:** 222, 140, 95… m/s.

### Semiimplícito ✅

- Estable con dt = 900 s, unas 6 veces el límite explícito.
- La masa sigue exacta; la energía varía 3·10⁻⁷ en 5 días por el amortiguamiento previsto de RAW sobre las ondas rápidas.
- Coste: **39 ms por paso** a 5° con 20 capas, en el entorno de la IA (con otros dos procesos en marcha).

### Prueba 4: equilibrio de JW06, 30 días (dt = 450 s) ✅

- **Desviación del viento:** |u − u₀| oscila entre 0,15 y 0,4 m/s (de 35) sin ir a más; v ≤ 0,4 m/s.
- **Presión en superficie:** se mantiene entre 999,8 y 1000,2 hPa.
- **Simetría:** la parte no zonal se queda en 7·10⁻¹⁵, el redondeo.
- **Conservación:** dE ≤ 10⁻⁸ y la masa exacta.
- **Origen de esa pequeña oscilación:** es el desequilibrio de la discretización. Disminuye con orden 2 al refinar la horizontal (9,7·10⁻⁶ → 2,4·10⁻⁶ → 6,3·10⁻⁷ con 160 capas) y también al refinar la vertical (20 → 40 capas: ÷3,6). Es decir, **las fórmulas están bien implementadas y el esquema es consistente** ✅.

### Prueba 5: onda baroclínica de JW06 (5°, dt = 450 s) ✅ cualitativo

| Día | p_s mínima (hemisferio norte) | Dónde | max\|v\| |
|---|---|---|---|
| 4 | 999,0 hPa | 47,5° N | 1,4 m/s |
| 6 | 995,9 | 47,5° N | 6,7 |
| 7 | 990,2 | 52,5° N | 14,8 |
| 8 | 979,4 | 52,5° N | 32,1 |
| 9 | 960,3 | 57,5° N | 56,2 |
| 11 | 941,5 | 62,5° N | 70,9 |

- La onda crece de forma explosiva entre los días 7 y 9, como describe JW06. Sin difusión y a 5° se mantiene estable 15 días.
- La masa sigue exacta.
- **dE baja hasta −2·10⁻⁶ en 15 días.** Es el filtro polar y RAW actuando sobre la energía que la onda envía a la escala de la rejilla; no es un error, pero es lo que §7 tendrá que tratar.
- ⚠️ Comparar los valores con las figuras de JW06, cuyas referencias de alta resolución no he podido leer, y con nuestra propia corrida a 2,5° (en marcha).

### Prueba 5 a 2,5° (dt = 225 s)

| Día | p_s mínima a 2,5° | p_s mínima a 5° |
|---|---|---|
| 8 | 973,8 hPa | 979,4 hPa |
| 9 | 952,3 | 960,3 |
| 10 | 937,6 | 952,5 |
| 12 | 911,4 | 946,0 |

- La borrasca crece antes y se hace más profunda con más resolución, que es la tendencia esperable.
- A 5° el sistema está subresuelto a partir del día ~9, pero la fase de crecimiento coincide.

### Prueba 6 (Held y Suarez): primer intento y diagnóstico (05/10)

**Intento 1: niveles de M3N, sin difusión.**
- Revienta el día 26. Ya el día 10 hay T_min = 156 K y |v| = 97 m/s, en las capas 0–2 (σ < 0,02).
- El número de Courant vertical, |W|Δt/Δp, sube de 0,1 a más de 4.

**Intento 2: 20 capas σ iguales, sin difusión.**
- Las 20 capas iguales son la especificación de la prueba ✅ (DCPAM: "20 layers, equally divided in σ").
- Aguanta más de 40 días, pero **el campo de viento se llena de ruido de rejilla**:
  - el 30 % de la energía de v está en las ondas zonales con k ≥ 24, en todas las latitudes;
  - la segunda diferencia meridional de v es de ~28 m/s, con un v rms de 16 m/s;
  - el Courant vertical de ~2 es un síntoma del ruido, no la causa.

**Conclusión: son dos problemas distintos.**
1. **El techo alto de M3N** (la capa de arriba va de 0 a 7,35 hPa, con el centro a 3,7 hPa; corregido el 06/10, antes decía "0,7 hPa" por un desliz de unidades) necesita una **capa esponja**, una disipación de Rayleigh en el tope, que es lo habitual en los modelos con el techo alto ⚠️ (fuente por verificar). **Queda para la 6.3**, porque la prueba publicada no usa esos niveles.
2. **La acumulación de energía y enstrofía en la escala de la rejilla** es lo que motiva la hiperdifusión de §7. No es un fallo del esquema: todas las implementaciones de HS94 la llevan (CESM y DCPAM ✅).

### Hiperdifusión implementada (§7, 05/10)

- **Esquema:** ∇⁴ sobre T (laplaciano escalar de la rejilla C) y sobre u y v (laplaciano vectorial, grad(div) − rot(rot)). En el polo, la vorticidad es la circulación de la fila de u más cercana dividida por el área del casquete.
- **Implícita y exacta:** (I + Δt ν₄ L²) x_nuevo = x, separando cada onda zonal (los operadores no dependen de la longitud). Así es estable con cualquier Δt, incluso junto al polo.
  - Las matrices se construyen aplicando los operadores a cada onda.
  - Comprobación: coinciden con el cálculo directo hasta 4–6·10⁻¹⁶ ✅.
- **Coeficiente:** ν₄ = 1/(τ·(4/Δx_ec²)²) con τ = 0,5 días, el criterio de CESM.
  - Resultado: ν₄ = 1,38·10¹⁷ m⁴/s.
  - Comprobación cruzada: CESM usa 1,17·10¹⁶ a T42, y pasar a la mitad de resolución multiplica el coeficiente por 16, lo que da ≈1,9·10¹⁷. Es el mismo orden ✅.
- **Uso en el tiempo:** se aplica sobre el salto de 2Δt del leapfrog, antes del filtro RAW. No actúa sobre p_s, así que la masa sigue exacta.
- **Lo que no hace:** no conserva exactamente la entalpía ponderada por Δp ni el momento angular (la rotación de sólido rígido se amortigua a 3·10⁻⁵/día, despreciable). La energía cinética que disipa no se devuelve como calor, como en HS94.

**Primeros 25 días con difusión (20 capas iguales):**
- ruido de Nyquist en T de 0,01 K (antes, 2–7 K);
- Courant vertical de 0,03–0,08;
- |v| ≤ 63 m/s en el ajuste inicial;
- T_min = 187 K.

Está en marcha la simulación completa de 1200 días.

**Rectificación (05/10, 17:40): con la hiperdifusión, los niveles de M3N también son estables.**
- Held y Suarez con los 20 niveles de M3N (capa de arriba de 0 a 7,35 hPa), τ = 0,5 días y dt = 450 s: **60 días estables y razonables**.
  - Courant vertical ≤ 0,04;
  - |v| ≤ 60 m/s solo en el ajuste inicial, y luego ~30 m/s;
  - T_min = 183–192 K;
  - ruido de Nyquist ≤ 0,1 K.
- **El fallo del intento 1 era también ruido de rejilla**, que en las capas altas, muy finas y poco densas, se manifestaba antes. No era una falta de capa esponja.
- **La capa esponja (DISENO_FASE6_3 §1.6) pasa de "imprescindible" a "evaluar en simulaciones largas"**, mirando si aparecen vientos irreales en las capas de arriba o reflexión de ondas.
- La forma de referencia, por si hace falta, ya está verificada en el código de MiMA/Isca (`damping_driver.f90`, "PK02-like sponge", líneas 594–634) ✅:
  - k = k_máx·((p_sp − p)/p_sp)² para p < p_sp (p_sp = 50 Pa por defecto);
  - opción de devolver la energía cinética perdida como calor: dT = −(u·du + v·dv)/c_p, con la corrección de medio paso.

### Prueba 6 (Held y Suarez) con 20 capas iguales y τ = 0,5 días: resultados parciales (05/10, días 200–360)

- **Estable** hasta el día 420 por lo menos; la simulación completa de 1200 días sigue en marcha, con puntos de control cada 20 días.
- **Media parcial en el hemisferio norte** (16 muestras, una cada 10 días):
  - **Chorro en altura:** máximo de **36,7 m/s hacia 27,5–32,5° a σ = 0,22**.
  - **Viento en superficie:** del oeste, con máximo de 7 m/s hacia 32°; alisios entre 8° y 18°.
  - **Transporte de calor por remolinos (v′T′ a σ = 0,87):** máximo hacia 32°.
  - **Transporte de momento por remolinos (u′v′ a σ = 0,25):** positivo entre 18° y 32° y negativo entre 42° y 58°.
- **Comparación con HS94:** según lo que recuerdo del artículo (⚠️, las cifras exactas no se han podido leer), el chorro debería quedar a ~30 m/s hacia ~45°. **La zona de borrascas entera sale unos 10–15° más cerca del ecuador.**
- **Explicación con fuente ✅:** Lu, Chen, Leung, Burrows, Yang, Sakaguchi y Hagos (2015), *J. Climate* 28, 6763–6782 (resumen leído en eesm.science.energy.gov). Con resolución gruesa, el chorro de las borrascas queda más cerca del ecuador, porque los modelos gruesos son "sobredisipativos" y las ondas no llegan a las latitudes donde rompen. A ~50 km la posición empieza a converger.
- **Conclusión provisional:** es coherente con un efecto conocido de los 5° y de la difusión, no con un fallo del núcleo. Las pruebas 1–5 son exactas o convergen con su orden teórico.
- **Pendiente para cerrarlo con rigor:**
  1. la media completa de los días 200–1200 y su incertidumbre (comparando las dos mitades);
  2. la sensibilidad a τ: repetirla con τ = 2 días, más débil, y ver si el chorro se desplaza hacia el polo;
  3. si se puede, 2,5° en el PC de Carlos;
  4. las cifras exactas de HS94.

### Prueba 6 (Held y Suarez): resultado completo en el entorno de la IA (05/10, 1200 días, τ = 0,5 d)

Media de los días 200 a 1200, con 100 muestras (una cada 10 días). Estable los 1200 días.

| | Hemisferio N | Hemisferio S |
|---|---|---|
| Chorro en altura | 37,2 m/s en 28,9°, σ = 0,22 | 37,1 m/s en 28,9°, σ = 0,22 |
| Viento del oeste en superficie | 3,3 m/s en 33,5° | 3,5 m/s en 33,5° |
| Calor por remolinos (v′T′, σ = 0,87) | máx. 14,7 K·m/s en 33,2° | 15,2 K·m/s en 33,2° |
| Momento por remolinos hacia el polo (u′v′, σ = 0,25) | máx. 32,9 m²/s² en 23,9° | 34,9 m²/s² en 24,0° |

- **Simetría:** los dos hemisferios, que la prueba hace idénticos, difieren como mucho 0,7 m/s en u medio. Eso indica que la media está bien convergida.
- **Temperatura de la capa baja:** 307 K en el ecuador y 265 K en el polo.

**Interpretación**
- **El núcleo es estable, conservativo y simétrico.** Produce una circulación con chorro, borrascas, alisios y vientos del oeste en superficie.
- **Toda la zona de borrascas queda unos 10–15° más cerca del ecuador** que en HS94 (~45° ⚠️, cifra por verificar).
- Es coherente con el efecto de la resolución gruesa y de la disipación que documentan Lu et al. (2015) ✅.

**Para cerrar la prueba (en el PC de Carlos, con `held_suarez.py`):**
1. la repetición con muestreo diario y la incertidumbre entre las dos mitades;
2. la sensibilidad a la difusión (τ = 2 días). Si el chorro se desplaza claramente hacia el polo, queda demostrado que la diferencia es de disipación y resolución, no del núcleo.

### Prueba 6 (Held y Suarez): resultado DEFINITIVO en el PC de Carlos (05/10, `held_suarez.py` + `analizar_held_suarez.py`)

1200 días, media diaria de los días 200 a 1200 (1000 muestras). Cada simulación tardó unos 25 minutos en el PC de Carlos. **El día 1200 de la corrida con τ = 0,5 coincide al dígito con la del entorno de la IA:** el núcleo es reproducible en dos máquinas distintas.

| | τ = 0,5 días (CESM) | τ = 2 días (difusión 4 veces más débil) |
|---|---|---|
| Chorro N | 37,3 m/s en 28,9°, σ = 0,22 | 40,6 m/s en 31,5°, σ = 0,22 |
| Chorro S | 37,2 m/s en 29,0° | 40,4 m/s en 31,3° |
| Incertidumbre (1.ª mitad / 2.ª mitad) | ±0,3 m/s, ±0,2° | ±0,2 m/s, ±0,1° |
| Viento del oeste en superficie | 3,3–3,5 m/s en 33,5° | 4,6–4,7 m/s en 35,3° |
| Calor por remolinos (máx.) | 33,1° | 33,5° |
| Momento por remolinos hacia el polo (máx.) | 33,7–35,0 m²/s² en 24° | 33,3–36,0 m²/s² en 25,3° |

- **Estabilidad:** las dos son estables durante 1200 días. Con τ = 2, T sigue entre 186 y 311 K y no aparece ruido.

**Conclusiones**
1. **El núcleo es correcto en lo que se puede comprobar exactamente** (pruebas 1 a 5). En Held y Suarez es estable, simétrico y reproducible, y produce chorros, borrascas, alisios y vientos del oeste en superficie.
2. **La incertidumbre estadística es despreciable** (±0,2°). La diferencia con HS94 (~45° ⚠️) es sistemática.
3. **Con 4 veces menos difusión, todo se desplaza ~2,4° hacia el polo y se refuerza ~3 m/s.** Es la dirección que predicen Lu et al. (2015) ✅, de modo que **la disipación explica una parte**.
4. **Lo que queda** (~13°) es coherente con la resolución de 5°, que equivale más o menos a T21, donde las borrascas apenas se resuelven. Lu et al. hablan de convergencia hacia ~50 km.
5. **Para demostrarlo del todo haría falta repetirla a 2,5°.** Costaría unas 8 veces más por día simulado (~3,5 h en el PC de Carlos) y exige que `held_suarez.py` acepte la resolución como opción (cambio pequeño). Queda como opcional.

**Decisión de diseño 🔶 que sale de aquí:** τ = 0,5 días es el valor conservador. La simulación con τ = 2 es estable y menos disipativa, y para la 6.3 conviene usar **la difusión más débil que siga siendo estable con la física completa**, medida y no supuesta.

**Figuras (PC de Carlos, `hs_tau0.5.png` y `hs_tau2.png`):**
- **Estructura cualitativa como la de HS94:** dos chorros simétricos en la troposfera alta (σ ≈ 0,2–0,3) que bajan hasta la superficie como vientos del oeste; vientos del este en superficie en los trópicos y en latitudes altas; vientos débiles del este sobre el ecuador en altura; y una tropopausa tropical fría (~190 K).
- **Con τ = 2** los chorros quedan algo más altos de latitud y más intensos.

**v3.1-pre4:** `held_suarez.py` acepta `--filas 72` (2,5°, dt = 225 s, misma τ en la onda más corta de esa rejilla), para la prueba de resolución. En el entorno de la IA va a ~1 día simulado por minuto; en el PC de Carlos se esperan ~3–4 h para 1200 días ⚠️.

### Prueba 6: prueba de resolución a 2,5° (PC de Carlos, 05–06/10) — CIERRE DE LA 6.2

| | 5°, τ = 0,5 | 5°, τ = 2 | **2,5°, τ = 0,5** |
|---|---|---|---|
| Chorro (N / S) | 37,3 / 37,2 m/s en 28,9 / 29,0° | 40,6 / 40,4 m/s en 31,5 / 31,3° | **35,4 / 34,9 m/s en 40,6 / 40,4°** |
| Incertidumbre (mitades) | ±0,2° | ±0,1° | ±0,3° (N), ±1,1° (S) |
| Viento del oeste en superficie | 3,3–3,5 m/s en 33,5° | 4,6–4,7 m/s en 35,3° | **7,1–7,6 m/s en 42,0–42,7°** |
| Calor por remolinos (máx.) | 33,1° | 33,5° | 36,4–36,7° |
| Momento por remolinos hacia el polo (máx.) | 34–35 m²/s² en 24° | 33–36 m²/s² en 25° | **67–68 m²/s² en 31°** |

**Conclusión: la 6.2 queda CERRADA ✅**
- **Al reducir la celda a la mitad, el chorro se desplaza ~11,5° hacia el polo (hasta ~40,5°) y se debilita hacia ~35 m/s.**
  - El viento del oeste en superficie se duplica y se coloca en ~42°.
  - El transporte de momento por remolinos se duplica.
  - Todo se acerca a HS94 (~30 m/s hacia ~45° ⚠️, cifra del artículo pendiente de verificar).
- **El desplazamiento era sobre todo un efecto de la resolución,** como describen Lu et al. (2015) ✅, y en menor medida de la disipación (τ = 2 a 5°: +2,4°). No era un fallo del núcleo.
- **El núcleo da la convergencia esperada al refinar la rejilla.**

**Consecuencias para M3N**
- **A 5°** (≈ 408 km en P3N), la zona de borrascas de P3N quedará previsiblemente unos 10° más cerca del ecuador de lo que daría un modelo de alta resolución, con unos vientos del oeste en superficie la mitad de intensos. Es una **limitación conocida y medida**, que se documentará en la validación de la 6.3 y en H3N.
- **Opciones a decidir más adelante** (con Carlos):
  - (a) aceptarlo;
  - (b) simular a 2,5°, unas 8 veces más caro;
  - (c) una corrección de los remolinos no resueltos, solo si la validación en modo Tierra lo justifica.
- ✅ **Decidido por Carlos (06/10/2026, reafirmado el 08/10):** de momento la rejilla se queda en 5°. La de 2,5° "ya habrá tiempo de pensarlo cuando estén decididas todas las características del planeta y podamos permitirnos una sola simulación larguísima" (`DISENO_FASE6_3.md` §6.3). La opción (c) solo se estudiaría si la validación en modo Tierra lo justificara.

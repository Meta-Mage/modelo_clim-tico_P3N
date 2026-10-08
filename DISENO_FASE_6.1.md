# Fase 6.1 de M3N: ecuaciones de aguas someras en la esfera (núcleo dinámico propio, primer ladrillo)

> **Nomenclatura (08/10/2026, M3N 3.12.1):** este documento se llamaba `DISENO_FASE6_1.md`. Fases y versiones van con la numeración normalizada; la correspondencia con los nombres anteriores está en `NOMENCLATURA.md` y `CHANGELOG.md`.

**Estado: código y validación hechos (05/10/2026). Las elecciones marcadas con 🔶 las ha tomado la IA de forma provisional y Carlos tiene que confirmarlas.**

Archivos:
- `fase6_aguas_someras.py`: el modelo y los casos de prueba;
- `validar_fase6_1.py`: el informe de validación;
- `test_fase6_1.py`: 9 pruebas automáticas.

Etiquetas: ✅ verificado; ⚠️ pendiente.

---

## 1. Qué es y para qué sirve

Una sola capa de fluido sobre una esfera en rotación: la versión más simple de la atmósfera que tiene ya los ingredientes difíciles de la dinámica, que son
- Coriolis;
- la fuerza de presión;
- la advección;
- las ondas de gravedad y de Rossby;
- el problema de los polos.

Sirve para construir y comprobar con soluciones conocidas la rejilla, los operadores y el avance en el tiempo que usará el núcleo de 20 capas (6.2). No se acopla todavía a M3N.

## 2. Discretización

| Pieza | Elección | Fuente |
|---|---|---|
| Rejilla | La de M3N, 72 × 36 de 5°. La resolución es un parámetro, solo para medir la convergencia. | — |
| Escalonado | 🔶 **Rejilla C**: h en el centro de cada celda, u en las caras este y oeste, v en las caras norte y sur. En los polos no hay v. | Arakawa y Lamb 1977 ✅ |
| Ecuaciones | Forma vectorial invariante: vorticidad potencial q en las esquinas, función de Bernoulli B = g(h+hs) + K en los centros. | Libros de la bibliografía (Vallis, Durran) |
| Término de Coriolis y vorticidad | Forma que conserva la energía, con los flujos de masa por cara. | Sadourny 1975 ⚠️ (la cita), con demostración propia (§3) |
| Masa | Forma de flujo: lo que sale de una celda entra en la vecina. | — |
| Polos | 🔶 Filtro de Fourier sobre las tendencias por encima de 60°: la onda zonal k se multiplica por S = min(1, cos φ / (cos 60° · sen(kΔλ/2))). La media zonal no se toca. | Remedio clásico de las rejillas latitud-longitud (planetWRF ✅) |
| Tiempo | Leapfrog con el filtro RAW (α = 0,53, ν = 0,2); el primer paso, Euler. | Williams 2009 ✅ |

## 3. Conservación: lo que está garantizado y por qué

- **Masa: exacta.** Cada flujo por una cara se resta de una celda y se suma a la vecina; el filtro polar no toca la media zonal. Medido: 10⁻¹⁶ (redondeo).
- **Energía total: exacta en el sistema semidiscreto** (sin filtro polar y con un paso de tiempo infinitesimal). Demostración propia, en tres partes que se cancelan:
  1. **Fuerza de presión contra la divergencia.** Con h en las caras igual a la media aritmética y K = (1/2A)·Σ_caras A_cara u²/2, el trabajo del gradiente de B en las caras es exactamente lo contrario de B·(cambio de h) en los centros.
  2. **Coriolis y vorticidad.** En cada esquina, la contribución de u es q/4·(U1+U2)(V1+V2) y la de v es la misma con signo contrario (U y V son los flujos de masa por cara). Se anula esquina a esquina.
  3. **Comprobación numérica:** con un estado con ruido y una montaña, la suma de las tres contribuciones es 10⁻¹⁴ de cada una por separado (prueba automática).
  4. **Lo que sí la cambia:** el paso de tiempo y el filtro polar. Medido: dE ~ 10⁻⁹ a 10⁻⁶ en las simulaciones.
- **Enstrofía potencial:** este esquema no la conserva exactamente. Medido: dZ ~ 10⁻³ en 14 días. Hay un esquema que conserva las dos cosas (Arakawa y Lamb 1981 ⚠️) y queda como mejora posible si en la 6.2 aparece ruido de escala de rejilla.

## 4. Resultados de la validación (5°, dt = 150 s, filtro desde 60°)

### Caso 2 (flujo zonal estacionario, solución exacta; error de h a los 5 días)

| Eje del flujo | l1 | l2 | l∞ |
|---|---|---|---|
| α = 0 | 1,6·10⁻⁴ | 1,8·10⁻⁴ | 3,3·10⁻⁴ |
| α = 0,05 | 3,1·10⁻⁴ | 5,5·10⁻⁴ | 2,8·10⁻³ |
| α = π/4 | 2,7·10⁻³ | 3,7·10⁻³ | 1,3·10⁻² |
| α = π/2 − 0,05 | 2,5·10⁻⁴ | 3,3·10⁻⁴ | 1,5·10⁻³ |
| α = π/2 (por encima de los polos) | 2,3·10⁻⁴ | 3,1·10⁻⁴ | 1,4·10⁻³ |

- **Convergencia (α = π/4):** l2 baja de 3,7·10⁻³ a 5° a 8,2·10⁻⁴ a 2,5° y a 1,97·10⁻⁴ a 1,25°, es decir, **orden 2,17 y luego 2,06** (se espera 2) ✅.
- **Un error mío que corrigió esta prueba:** en el caso 2 con el eje girado hay que girar también el eje de rotación del planeta, f = 2Ω(−cos λ cos φ sen α + sen φ cos α) (Williamson ✅ vía SWEET). Sin eso salía un error del 12 %.
- Comparación con modelos publicados: ⚠️ pendiente. Las tablas de Williamson (1992) usan otros modelos y resoluciones, y su escaneo no se lee bien.

### Caso 5 (montaña, 15 días)

- Masa exacta.
- dE = −1,4·10⁻⁷ y dZ = −8,4·10⁻⁴.
- Estable, con la superficie libre entre 5035 y 5954 m.

### Caso 6 (onda de Rossby-Haurwitz, 14 días)

- **Masa exacta y energía con dE = −1,9·10⁻⁶.**
- **Fórmula de la altura:** la transcribí del escaneo, así que la verifiqué por mi cuenta. El desequilibrio inicial baja como 1/4 al dividir la celda por 2 (0,18 → 0,048 → 0,012 → 0,003), así que la fórmula es la correcta ✅.
- **Estudio de convergencia (05/10):** el caso 6 se ha repetido a 5°, 2,5° y 1,25° (dt = 150, 75 y 37,5 s).

| Rejilla | Velocidad de la onda (°/día) | Amplitud de la onda 4 hacia 50° N (relativa a la inicial) | Diferencia de h a los 14 días con la de 1,25° (l2) |
|---|---|---|---|
| 5° | 11,409 | máximo 1,19 el día 12; 1,14 el día 14 | 1,8·10⁻² |
| 2,5° | 11,446 | máximo 1,14 el día 6; 1,06 el día 14 | 1,2·10⁻³ |
| 1,25° | 11,422 | máximo 1,13 los días 6–8; 1,06 el día 14 | — |

- **Velocidad de la onda:** ya está convergida, en **11,43 ± 0,03°/día**. A 5° la desviación es de solo −0,02°/día. La fórmula barotrópica no divergente da 12,20°/día, y la diferencia del −6 % es física: en aguas someras el fluido puede converger y divergir, y la onda se frena (Williamson et al. 1992 dice que el caso 6 solo resuelve "aproximadamente" las aguas someras; lo cita así el manuscrito de GMD 2014-126 ✅). No he podido leer un valor publicado de aguas someras: Thuburn y Li (2000, Tellus A) no ha sido accesible. Queda como comparación pendiente y no bloquea nada, porque el valor convergido es una propiedad de las ecuaciones y no de nuestra rejilla.
- **Amplitud:** la subida de ~13 % hacia los días 6–8 seguida de una bajada **es física**: sale igual a 2,5° y a 1,25° ✅. La rejilla de 5° la retrasa y la exagera un poco (máximo de 1,19 el día 12), un efecto de la resolución que es esperable y queda medido.

### Caso 1 (advección de una campana)

- **El esquema centrado de segundo orden dispersa mucho una campana de solo ~8 celdas de diámetro:** l2 = 0,34 en un día a 5° y ~1 a los 12 días, con valores negativos de hasta −450 m.
- **Lo que pasa, en cada prueba:**
  - el error no cambia con dt (de 300 a 60 s) ni con el filtro, así que es error espacial;
  - converge con la resolución: 0,34 a 5°, 0,105 a 2,5° y 0,031 a 1,25°.
- **Conclusión:** para h, que es un campo suave, no importa. Pero confirma lo que ya preveía el borrador (§3.4): **el vapor y los demás trazadores necesitarán un esquema de volúmenes finitos con limitador**, que no dé valores negativos.

## 5. Hallazgo: el límite de estabilidad del filtro RAW

- **El síntoma:** el caso 6 reventaba al noveno día con dt = 300 s, aunque el análisis lineal decía que el paso máximo era de unos 430 s.
- **Cómo encontré la causa:**
  - Calculé los autovalores del sistema linealizado (método de Arnoldi, con el producto por el jacobiano por diferencias finitas).
  - El modo más rápido es una onda de gravedad hacia 62,5°, en el borde del filtro, con ω = 2,36·10⁻³ s⁻¹. Sus autovalores son imaginarios puros: el operador espacial es neutro y no crece.
  - La culpa la tiene el **filtro RAW**. Calculé el factor de amplificación de leapfrog + RAW (α = 0,53, ν = 0,2) sobre una oscilación:

| ω·dt | 0,1 | 0,3 | 0,45 | 0,5 | 0,7 | 0,9 |
|---|---|---|---|---|---|---|
| Amplificación por paso | 0,99997 | 0,99982 | ≈1 | 1,0004 | 1,0056 | 1,063 |

- **El cálculo cuadra con el síntoma:** con dt = 300 s, ω·dt = 0,71 y la onda crece un 0,56 % por paso. En 9 días son 2592 pasos, es decir, un factor de e¹⁴: justo lo que se vio.
- **Comparación:** el leapfrog sin filtro es neutro hasta 1, y Robert-Asselin (α = 1) es estable hasta ~0,9, pero amortigua mucho más.
- **Decisión 🔶:** se mantiene RAW, porque conserva mucho mejor la amplitud de las ondas lentas, y en la 6.1 se usa **dt = 150 s** (ω·dt = 0,35). Queda como prueba automática.
- **Consecuencia para la 6.2:** con las ondas de gravedad explícitas, el paso sería de ~150 s, muy lejos de los 992 s de la física de M3N. El esquema **semiimplícito** previsto en el borrador, que trata implícitamente las ondas de gravedad, no es opcional: es lo que permite llegar a 992 s. Después de él, el límite lo pone el viento: 0,45 / (velocidad / Δx).
- ⚠️ No he encontrado en la fuente que Williams (2009) documente esta pérdida de estabilidad; es un resultado de cálculo propio y reproducible (`test_limite_de_estabilidad_de_raw`).

## 6. Coste

A 5°, un paso cuesta ~1 ms (numpy, sin numba). Los 14 días del caso 6, con 8064 pasos, tardan ~10 s en el entorno de la IA.

## 7. Decisiones de Carlos (05/10/2026, 16:12)

Carlos confirma las cuatro propuestas: **rejilla C**, **filtro polar desde 60°**, **RAW con paso corto en la 6.1 y semiimplícito obligatorio en la 6.2**, y **se pasa a la 6.2** después de cerrar los ⚠️ de la 6.1. Pide que todo sea óptimo, pero sobre todo preciso y muy riguroso.

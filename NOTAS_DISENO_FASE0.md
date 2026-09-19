# Notas de diseño — Fase 0 (rejilla vectorizada)

Fecha: 20/09/2026. Resume las decisiones tomadas al vectorizar el modelo
sobre una rejilla lat/lon real, para quien retome esto más adelante sin
tener que releer toda la conversación en la que se hizo.

## Resolución de trabajo
72 columnas × 36 filas (5° por celda). Misma convención de "centro de
celda" que usa C3N en `mapa.js` (`lat = 90 - (fila+0.5)*grados_por_fila`,
`lon = (col+0.5)*grados_por_col - 180`), para que el futuro puente
C3N → M3N no tenga que traducir coordenadas.

## Qué se reutilizó sin tocar, y por qué
- `declinacion_solar()`: no depende de la posición, se llama igual que antes.
- `angulo_horario()`: solo hace aritmética, ya funciona sobre arrays de NumPy sin ningún cambio.
- `i_toa()`, `i_atm()`, `i_abs()`, `t_eq()`: aritmética pura, funcionan igual sobre arrays.
- `precalcular_orbita()`: se reutiliza tal cual. Aprovecha que `angulo_horario()` es lineal en la longitud (`angulo(hora, lon) = angulo(hora, 0) + lon` en radianes), así que el ángulo horario por columna se obtiene sumando la longitud a lo que `precalcular_orbita()` ya calculaba — no hizo falta rehacer los cálculos orbitales por columna.

## Qué se creó nuevo, y por qué
- `angulo_cenital_rejilla()`: `angulo_cenital()` usa `math.sin/cos/acos`, que no aceptan arrays.
- `masa_aire_rejilla()`, `trans_rejilla()`: mismo motivo (`math.sin`, `math.exp`).
- `i_inst_rejilla()`: mismo motivo (`math.cos`).
- `simular_rejilla()`: el bucle de acumulación con inercia térmica, adaptado para toda la rejilla a la vez.

## Decisiones deliberadas (no descuidos)
- **Noche = `np.nan`, no `None`.** `None` no cabe en un array de NumPy. `np.nan` se propaga solo por la cadena (`trans`, `i_atm`, `i_abs`) y se convierte a 0 al final con `np.nan_to_num`.
- **Convergencia conjunta.** Como todas las celdas se calculan a la vez, no tiene sentido que cada una "salga" del bucle por separado. Se sigue iterando hasta que la celda que MÁS ha cambiado en el año baje de la tolerancia (0.01°C).
- **No se guarda el ciclo horario completo por celda.** Guardar el año entero (~90.000 valores) para 2.592 celdas sería un volumen de memoria grande sin que se haya pedido. `simular_rejilla()` devuelve solo la temperatura final de cada celda. Si hace falta el detalle horario de alguna celda concreta más adelante, se puede añadir sin tocar esto.
- **Instante inicial de noche → 273.15 K**, igual que hace la versión de un punto (`t_eq(0)` daría 0 K, sin sentido físico).
- **Recorte defensivo en `arccos`** (`np.clip` a [-1, 1]): con miles de celdas a la vez, un error de redondeo minúsculo puede sacar `cos_cenital` de rango y romper `arccos`. Con un único punto nunca se nota.

## Validación
- Física de un instante (ángulo cenital → masa de aire → transmitancia → irradiancia absorbida): validada celda a celda contra el modelo de un punto, diferencia del orden de 1e-13 (ruido de redondeo puro).
- Simulación completa con inercia térmica: validada en 5 celdas representativas (polos, latitudes medias, ecuador, extremos de longitud) contra el modelo de un punto. Diferencia máxima 6.3e-3 °C — mayor que cero pero por debajo de la tolerancia de convergencia (0.01°C), y concentrada en los polos, que son las celdas más lentas en converger (la rejilla espera a la más lenta de todas; el punto aislado paraba antes).
- Sensibilidad al paso de tiempo (900 s vs. 450 s): diferencia máxima 0.30°C (en polos), media 0.08°C en toda la rejilla. Pequeña frente al rango de temperaturas del modelo y frente a sus propias simplificaciones físicas (invernadero de una capa, inercia térmica reducida a un número). Decisión: mantener 900 s por ahora; revisar si en el futuro la precisión importa más (p. ej. con geografía real).
- Todo lo anterior está también como suite de `pytest` en `test_fase0.py` (`pytest test_fase0.py -v`).

## Pendiente, sin urgencia
- Limpiar los `print()`/cálculo de un punto fijo que se ejecutan solos al importar `geometria.py`/`atmosfera.py`/`radiacion.py` (ruido heredado de antes de la Fase 0, no afecta al resultado).
- `visualizacion.py` está roto (llama a `simular()` con argumentos que ya no coinciden) y no se usa en ningún flujo actual — decidir si se borra.
- Cuando llegue la Fase 1 (geografía real), cada celda necesitará su propio albedo/inercia según sea agua o tierra — hoy `albedo` y `profundidad_optica` siguen siendo un único valor para toda la rejilla.

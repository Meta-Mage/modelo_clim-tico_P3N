# M3N — Fase 4: exportación de datos estructurados

> **Nomenclatura (08/10/2026, M3N 3.12.1):** este documento se llamaba `DISENO_FASE4.md`. Fases y versiones van con la numeración normalizada; la correspondencia con los nombres anteriores está en `NOMENCLATURA.md` y `CHANGELOG.md`.

**Estado:** implementado en `exportar_clima.py` y en el registro horario de `fase2b_atmosfera.py`.
**Fecha:** 03/10/2026 · **Versión:** `v2.4` (correcciones en la `v2.4.1`: hitos del año en el día correcto, codificación de valores extremos, altitud 0 en las celdas de agua, huella documentada)

---

## 0. Resumen

La Fase 4 es el **punto donde M3N empieza a hablar con C3N y H3N**. Hasta ahora el puente iba en un solo sentido: C3N escribía el mapa activo y M3N lo leía. Ahora M3N escribe su clima en un archivo que cualquier herramienta puede leer sin saber nada de Python ni de la física:

```
~/Documentos/B3N/mapa_activo/mapa_activo_m3n.json    C3N → M3N   (geografía)
~/Documentos/B3N/clima_activo/clima_activo_m3n.json  M3N → H3N   (clima)
```

Uso:

```
python exportar_clima.py
```

Simula el mapa activo de C3N, o lo recupera de la caché si ya estaba simulado, y escribe el archivo. Con la rejilla actual ocupa del orden de 40–60 MB.

---

## 1. Decisiones (Carlos, 03/10/2026)

| Decisión | Detalle |
|---|---|
| Solo JSON | Un único archivo `.json`. Las series grandes van **comprimidas dentro** (si no, los datos horarios pesarían cientos de MB) |
| Sin calendario | El tiempo se expresa como "día N / 270". Los hitos del año (perihelio, afelio, solsticios y equinoccios) van en el archivo |
| Día de 24 horas | `ROTACION_PERIODO = 86400` s. La hora de referencia es la del **meridiano 0** |
| Datos horarios en todas las celdas | Aire a 2 m y superficie, cada hora de los 270 días |
| Rejilla | 72 × 36 (5°) por ahora |

---

## 2. Registro horario (cambio en el modelo)

En el año final, el que se registra, `simular_fase2b` guarda además del resumen diario el valor **instantáneo** a cada hora en punto (hora del meridiano 0), en todas las celdas:

- `horario_aire2m`: aire a 2 m, °C, forma (270 × 24, 36, 72);
- `horario_superficie`: la temperatura real de la superficie (suelo, agua o la del hielo donde lo hay), °C.

**Detalles:**
- "Instantáneo a la hora en punto" es el criterio de las observaciones meteorológicas sinópticas: el valor en ese momento, no la media de la hora.
- El registro horario no cambia la física ni ningún otro resultado.
- Necesita que el paso de tiempo divida una hora exacta (900 s sí). Si no, los campos horarios valen `None` y la exportación avisa.
- **Caché:** la clave pasó de `simular_fase2b_v1` a `simular_fase2b_v2` en la v2.4, porque las simulaciones guardadas antes no tienen el registro horario. En la v2.4.1 pasa a `v3` e incluye además la huella del código fuente de la física. La primera simulación tras actualizar es completa (unos 20 minutos).

---

## 3. Formato `m3n-clima`, versión 1

```jsonc
{
  "formato": "m3n-clima",
  "version_formato": 1,
  "generado": "2026-10-03T19:30:00",
  "m3n":   { "version": "v2.4", "interruptores": { "...": true } },
  "mapa":  { "nombre": "prueba1", "huella": "3f2a…" },      // huella de la geografía simulada
  "rejilla": { "filas": 36, "columnas": 72, "latitudes": [...], "longitudes": [...] },
  "tiempo": {
    "dias": 270, "horas_por_dia": 24, "duracion_dia_s": 86400, "duracion_año_dias": 270.11,
    "referencia_hora": "meridiano 0 (hora solar local = hora + longitud/15)",
    "fechas_clave": [ { "nombre": "Perihelio", "dia": 255 }, … ]   // dia: índice 0..269
  },
  "astronomia": {                  // a la hora 12 del meridiano 0 de cada día
    "declinacion_solar_grados": [...270], "flujo_toa_W_m2": [...270], "distancia_S3N_UA": [...270]
  },
  "parametros": { "masa_S3N_soles": 0.874, … },
  "simulacion": { "anos_hasta_convergencia": 22, "saltos": 0, "balance_energia": 1.7e-5, "paso_tiempo_s": 900 },
  "celdas": { "tipo": ["agua", …], "altitud_m": [...], "nieve_permanente_posible": [...] },  // 2592, fila a fila
  "diario":  { "aire2m_min": SERIE, "aire2m_media": SERIE, "aire2m_max": SERIE,
               "superficie_min": SERIE, "superficie_media": SERIE, "superficie_max": SERIE,
               "hielo_espesor_m": SERIE },                       // forma [270, 36, 72]
  "horario": { "aire2m": SERIE, "superficie": SERIE }            // forma [6480, 36, 72]
}
```

**SERIE** (las series grandes):

```json
{ "codificacion": "int16-centesimas-zlib-base64", "forma": [270, 36, 72], "datos": "<base64>" }
```

Para leerla:
1. Descodificar el base64.
2. Descomprimir con zlib.
3. Leer enteros de 16 bits (little-endian), en orden C: el último índice es la columna.
4. Dividir entre 100. El resultado está en °C, o en metros en el caso del hielo.

El valor −32768 significa "sin dato". La precisión es de 0,01 °C, más que suficiente para el modelo. El lector de referencia en Python es `leer_serie()`, dentro de `exportar_clima.py`.

**Índices:**
- Fila 0 = la más al norte (87,5° N); columna 0 = 177,5° O. Es la misma rejilla de M3N y el mismo orden que C3N.
- Día `d` = 0…269, es decir, el día d+1.
- Hora `k = d·24 + h`, a las h:00 del meridiano 0. La hora solar local de una celda es `h + longitud/15`, y la longitud subsolar es `−(h − 12)·15`.

**Huella del mapa:** identifica la geografía exacta simulada. H3N la compara con la del mapa activo **actual** para saber si el clima está al día. Receta (`huella_mapa()` en `exportar_clima.py`):
1. SHA-256 de los bytes del tipo de celda (36 × 72, enteros de 8 bits, orden C),
2. seguidos de los bytes de la altitud en metros redondeada a 3 decimales (36 × 72, float64 little-endian, orden C),
3. los 16 primeros caracteres hexadecimales del resultado.

La altitud es la que lee `puente_c3n.py`, es decir, con las celdas de agua a 0 m (v2.4.1).

**Momentos del registro:**
- El **registro horario** toma el valor a las h:00 en punto (antes del paso que empieza a esa hora): de 00:00 a 23:00.
- El **registro diario** (mín., media, máx.) se calcula con los valores al final de cada paso: de 00:15 a 24:00.

La diferencia, un paso de 15 minutos, es despreciable, pero se documenta.

**Fin del año:** el año de P3N dura 270,11 días: 270 días completos más 11 pasos de 15 minutos (2 h 45 min). Esos pasos se simulan, pero no se registran como un día 271. Al empezar el año siguiente, el ángulo horario vuelve a la hora 0 del meridiano 0, así que cada año "repite" unas 2,75 horas. Es un efecto despreciable de la discretización en días enteros.

**Versión del formato:** si en el futuro se añaden variables (precipitación en la Fase 5, biomas en la Fase 7), se añaden claves nuevas sin cambiar las existentes y la versión sigue igual. Si cambia algo que ya existía, sube la versión, y los lectores rechazan las versiones que no conocen.

---

## 4. Validación

| Prueba | Criterio | Resultado |
|---|---|---|
| El registro horario no cambia el resto | Resultados diarios idénticos a la v2.3 con el mismo mapa | ✅ Diferencia 0,0 en mínima, media y máxima del aire, suelo y hielo (mapa `prueba1`, 22 años en los dos casos) |
| Ida y vuelta del formato | `leer_serie(serie(x))` = x con error ≤ 0,005 | ✅ Error máximo 0,005 (el redondeo a centésimas) |
| Coherencia hora/día | La media de las 24 horas ≈ la media diaria registrada | ✅ Diferencia media 0,001 °C, máxima 0,07 °C (24 muestras frente a 96 pasos) |
| Ciclo diario | En tierra, la máxima horaria cae por la tarde en hora solar local | ✅ Mediana a las 15,8 h solares locales (entre 15,5 y 16,2 h), como en la Tierra, donde suele darse entre las 14 y las 16 h |
| Tamaño | Manejable | 49,4 MB con la rejilla actual |
| Hitos del año (v2.4.1) | Cada hito cae en el registro del día que lo contiene | ✅ Solsticio de verano en el día 138 (índice 137) y equinoccio de otoño en el 208 (índice 207); en la v2.4 salían un día tarde. Los otros cuatro no cambian |
| Exportación completa (v2.4.1) | Con todas las correcciones, el clima sigue siendo el mismo | ✅ Mapa `prueba1`: aire a 2 m medio 15,29 °C (15,30 en la v2.4), océano helado 8,75–9,1 %, convergencia en 25 años; huella del mapa `2a6fe3e0cef76ac6` (cambia respecto a la v2.4 porque el agua tiene ahora altitud 0) |

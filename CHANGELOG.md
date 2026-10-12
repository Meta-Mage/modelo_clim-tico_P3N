# Registro de cambios de M3N

Formato [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/), con versiones [SemVer](https://semver.org/lang/es/). Las reglas y las correspondencias con los nombres antiguos están en `NOMENCLATURA.md`. Va de la versión más reciente a la más antigua.

En cada entrada:
- **Fase** es la del plan de trabajo (NOMENCLATURA §2).
- **Clima exportado** dice si cambia el clima que M3N exporta por defecto para H3N.
- **Pruebas** es el número de pruebas automáticas que pasan; "no consta" si no está documentado.

## [3.20.0] — 2026-10-12

- Fase: 6.4 (resultados de la simulación de P3N e investigación de los vientos del oeste) y 6.3 (lo que la validación no miraba).
- Clima exportado: no cambia (el modelo no cambia: con los valores por defecto, idéntico bit a bit a la 3.19.0).
- Pruebas: 144.

### Añadido
- `held_suarez.py --planeta p3n`: la prueba de Held y Suarez con el radio, la gravedad y la rotación de P3N y el mismo forzamiento (`DISENO_FASE_6.4.md` §8.4); `held_suarez.constantes_planeta`. La salida dice el planeta y su geometría. `analizar_held_suarez.py` dice el planeta de cada simulación y compara P3N con la Tierra de la misma configuración. Por defecto, todo como antes.
- `herramientas/control_tierra_p3n.py`: aplica al modo Tierra las mismas funciones que `analizar_p3n.py` aplica a P3N (todas las células, la función de corriente a σ 0,8 hacia los polos, las bandas de viento de la capa baja y el viento ecuatorial a σ ~0,25).
- `DISENO_FASE_6.4.md` §8: la simulación definitiva de P3N (PC de Carlos, 11–12/10/2026: 7,85 h, equilibrio en 15 años, 30 de climatología): los bloqueantes A1–A6 se cumplen; el clima de P3N; lo que no era lo esperado; el control con el modo Tierra; y la investigación de los vientos del oeste desplazados al polo (decisión de Carlos del 12/10/2026), con cómo se leerá fijado antes.
- `DISENO_FASE_6.3.md` §6.26: sin célula polar y con viento del oeste en la alta troposfera ecuatorial también en el modo Tierra (no lo medía `validar_i16.py`); limitaciones declaradas, sin reabrir el cierre.
- 3 pruebas (`test_held_suarez_planeta.py`).

### Cambiado
- `DISENO_FASE_6.4.md`: el estado, la prueba de arranque (hecha), el "1 bar" y el relieve (decididos el 11/10, §7), el criterio §4.7 (aprobado) y el plan del §5, al día.
- `README.md`: la Fase 6.3, cerrada; la 6.4, en curso; las limitaciones de la circulación del núcleo; las herramientas nuevas.
- `NOMENCLATURA.md`: commit de la v3.19.0 (`2788ad5`, leído de GitHub) y paquete 031.

## [3.19.0] — 2026-10-11

- Fase: 6.3 (cierre) y 6.4 (diseño y herramientas).
- Clima exportado: no cambia (el oficial sigue siendo el modelo de 2 capas; con los valores por defecto, el modelo con I16 es idéntico bit a bit a la 3.18.1).
- Pruebas: 141.

### Añadido
- `herramientas/criterio_cierre_6_3.py`: aplica el criterio de cierre de la Fase 6.3 (`DISENO_FASE_6.3.md` §6.22) a un `resumen.txt` y un `informe.txt`: tabla de bloqueantes, estructura e intensidades, y lo que necesita explicación.
- **Cierre de la Fase 6.3** (`DISENO_FASE_6.3.md` §6.25, aprobado por Carlos el 11/10/2026): el criterio de §6.22 aplicado a la climatología de la 3.18.1 del modo Tierra (PC de Carlos, 9,71 h, equilibrio en 8 años). Los siete bloqueantes se cumplen; 15 de 19 criterios de estructura; cada fallo y cada intensidad fuera de 0,5–2, explicados con fuente y ligados a una limitación declarada.
- `DISENO_FASE_6.4.md`: el diseño de la Fase 6.4 (P3N con el núcleo dinámico, la 4.0.0), con las decisiones de Carlos del 11/10/2026: inclinación 1,7°, océano con valores provisionales y sensibilidad, composición del aire como la Tierra, "1 bar" de aire seco, relieve más adelante, criterio para la 4.0.0 sin resultados exigidos y vista rápida anotada como idea; y la prueba del arranque de P3N (§4.3).
- `analizar_p3n.py`: análisis del clima de P3N con I16 sin suponer su forma (todas las células, chorros, bandas de viento y máximos de lluvia; estaciones por la distancia a la estrella; bloqueantes; orientación de Held y Hou; comparación con el clima oficial de 2 capas; sensibilidades) (`DISENO_FASE_6.4.md` §5.1).
- `clima_dinamico.py --oceano-difusion F --oceano-capa-mezcla F` (con `--anos`): pruebas de sensibilidad del océano en su propia carpeta, sin poder mezclarse con otra configuración ni exportarse como el clima de H3N (`DISENO_FASE_6.4.md` §4.2). `fase6_clima.simular_clima(..., opciones_nucleo=...)`. La exportación da los factores en `climatologia.sensibilidad` y el D_oc usado.
- `validar_i16.ano_extra` acepta otro calendario de temporadas (por defecto, los meses del modo Tierra, como antes).
- 10 pruebas: `test_analizar_p3n.py` (8) y `test_sensibilidad_oceano.py` (2).

### Cambiado
- `NOMENCLATURA.md`: la Fase 6.3, cerrada; la 6.4, en curso; commits de la v3.15.0 a la v3.18.1 (leídos de GitHub); paquete 030.

## [3.18.1] — 2026-10-11

- Fase: 6.3.
- Clima exportado: no cambia (solo cambia el modelo con I16).
- Pruebas: 131.

### Corregido
- El τ fijo del infrarrojo, recalibrado con el código actual por el método de la v3.0 en la columna observada: τ₀ = 2,4450, f = 0,2326 (antes 2,452 y 0,240, que con el esquema actual daban OLR 265,9 W/m² en vez de 267). Solo con I16 y `COLUMNA_OBSERVADA_I16` (decisión de Carlos del 11/10/2026; `DISENO_FASE_6.3.md` §6.24).

### Añadido
- `DISENO_FASE_6.3.md` §6.22: el criterio escrito de cierre de la Fase 6.3 (aprobado por Carlos el 11/10/2026).
- `DISENO_FASE_6.3.md` §6.23: la climatología de la 3.17.0 (efecto de C_E: +0,8 % de precipitación y −0,04 K; cadena limpia de 9,9 h) y una vista previa del criterio.
- 1 prueba (el τ recalibrado da los flujos de Wild et al. 2019).

## [3.18.0] — 2026-10-09

- Fase: 6.3.
- Clima exportado: no cambia (el modelo de 2 capas no usa nada de esto; cambia el modelo con I16).
- Pruebas: 130.

### Cambiado
- Con I16, el vapor es **razón de mezcla** (kg por kg de aire seco) en toda la física, como el núcleo de masa seca que lo transporta: saturación m_s = ε·e_s/(p − e_s) con la presión **total** (la seca más el peso del vapor de encima, como la ec. 66 de Lauritzen et al. 2018), Betts-Miller sin las conversiones q ↔ r de Isca y con la masa seca de cada capa, condensación y evaporación. La exportación sigue dando la humedad específica (m/(1 + m)) y la relativa como e/e_s. Interruptor `VAPOR_RAZON_MEZCLA_I16` (decisión de Carlos del 09/10/2026; `DISENO_FASE_6.3.md` §6.20.1 y §6.21).
- Con I16, los espesores ópticos se normalizan con la **columna media observada** de la Tierra (983,05 hPa de aire seco; Trenberth y Smith 2005) en vez de con 1000 hPa, y el **modo Tierra** tiene esa masa de aire seco (la presión al nivel del mar se calcula con el mapa). P3N sigue con 1 bar al nivel del mar. La antigua `P0` se separa en `P_REF` (referencia de θ, 1000 hPa), `P_NIVEL_MAR` (masa) y `P_COLUMNA_OBSERVADA`. A y B del vapor radiativo (I15, prototipo) recalibrados con la columna observada. Interruptor `COLUMNA_OBSERVADA_I16` (§6.20.2 y §6.21).
- Con los dos interruptores apagados, idéntico bit a bit a la 3.17.0.

### Añadido
- `calibrar_columna_observada.py`: reproduce las calibraciones del infrarrojo con la columna de 1000 hPa y las repite con la observada.
- 7 pruebas (`test_vapor_y_masa.py`). La del punto de control con I16 compara 2 "años" (antes 3): el tercero de ese estado irreal está al borde de la inestabilidad (§6.21).

### Corregido
- `DISENO_FASE_6.3.md` §6.19 y la entrada de la 3.17.0 de este registro daban mal el error de la saturación (+1,15 %): es −1,7 % en los trópicos (§6.20.1).

## [3.17.0] — 2026-10-09

- Fase: 6.3.
- Clima exportado: no cambia (el oficial sigue siendo el modelo de 2 capas, que no usa C_E; cambia el modelo con I16).
- Pruebas: 123.

### Cambiado
- Con I16, la evaporación y la sublimación sobre el agua usan C_E con la rugosidad del vapor del ECMWF (z₀_q = 0,62 ν/u*), algo mayor que C_H (C_E/C_H ≈ 1,04 en neutro); sobre tierra, C_E = C_H exactamente. Interruptor `EVAPORACION_Z0Q_I16` (decisión de Carlos del 09/10/2026; `DISENO_FASE_6.3.md` §6.19). Con el interruptor apagado, idéntico bit a bit a la 3.16.0. Deja anticuada la climatología del modo Tierra hasta la próxima simulación limpia.

### Corregido
- La exportación escribía NaN (no válido en JSON) en la deriva de una climatología muy corta; ahora escribe `null` y nunca NaN (`allow_nan=False`).

### Añadido
- `DISENO_FASE_6.3.md` §6.19: el vapor en la presión y la P₀ del modo Tierra necesitan un diseño previo (la saturación usa la fórmula de la humedad específica con la presión seca, mientras el balance trata q como razón de mezcla: ~~+1,15 % en q_s con e = 30 hPa~~ — cifra errónea, corregida en la 3.18.0: −1,7 % en los trópicos; P₀ tiene a la vez el papel de referencia de θ y de masa de aire, y la calibración del infrarrojo depende de ella).
- 2 pruebas (C_E sobre el agua y la tierra; JSON sin NaN).

## [3.16.0] — 2026-10-09

- Fase: 6.3 (y 4: claves nuevas de `m3n-clima` v1).
- Clima exportado: no cambia (el oficial sigue siendo el modelo de 2 capas; el modelo con I16 es idéntico bit a bit).
- Pruebas: 121.

### Añadido
- `exportar_clima.py --dinamico`: exporta la climatología con el núcleo dinámico (medias de N años, lluvia, nieve, evaporación, viento y humedad de la capa baja, extremos absolutos, desviación entre años, campos anuales, climatología y unidades) en `m3n-clima` v1 con claves nuevas (`DISENO_FASE_4.md` §5). Hasta la 4.0.0 va a `clima_dinamico_m3n.json`; con `--activo`, al oficial. Comprobado que H3N v0.2.1 lo lee sin cambios.
- `clima_dinamico.py` guarda en `resultado.pkl` el nombre y la huella del mapa; la exportación rechaza otro mapa.
- `DISENO_FASE_6.3.md` §6.18: análisis con fuentes de C_E frente a C_H (ECMWF; Large y Yeager 2004), de la masa del vapor (F_ps = 0: ~3,5 hPa, 0,35 % en q_s; Trenberth y Smith; Lauritzen et al. 2018) y de la presión del modo Tierra (~1 % más de aire seco que la Tierra), pendientes del visto bueno de Carlos.
- La validación completa de prueba comprueba también la exportación.

## [3.15.0] — 2026-10-09

- Fase: 6.3.
- Clima exportado: no cambia (el modelo es idéntico bit a bit; climatología, diagnósticos y documentación).
- Pruebas: 121.

### Añadido
- Climatología de 30 años fijos (normal climatológica estándar de la OMM, WMO-No. 1203 §3); `clima_dinamico.py --anos N` (N ≥ 12, WMO-No. 1203 §5.2.3) para una climatología exploratoria en su propia carpeta, marcada en el resumen.
- Comprobación de deriva a posteriori en el resumen (`fase6_clima.deriva_climatologia`): tendencia de Santer et al. (2000) del aire a 2 m global, de cada banda y de la precipitación global, con el contraste múltiple de Benjamini y Hochberg (1995).
- Medias diarias del viento (u, v, rapidez) y la humedad (específica y relativa) de la capa baja con I16, acumuladas en la climatología (`DISENO_FASE_6.3.md` §6.17).
- Fase 6.6 (arrastre de ondas de gravedad que conserva el momento, después de la 5.2) en el plan; sustituye a la capa esponja prevista (Shaw y Shepherd 2007).
- `DISENO_FASE_6.3.md` §6.17: las decisiones del 09/10 (umbrales, estratosfera, 4.0.0 provisional tras la 6.4) y su justificación con fuentes.
- 2 pruebas (tendencia de Santer y Benjamini-Hochberg; deriva) y comprobaciones del viento y la humedad en la de la climatología.

### Cambiado
- El criterio de error de la media ya no decide cuántos años se registran; el error de cada banda se informa.
- Volver a lanzar `clima_dinamico.py` con una climatología terminada rehace el resumen sin simular.

### Corregido
- `DISENO_FASE_6.3.md` §6.16: las cifras del ruido del equilibrio salieron de años sin equilibrar (σ 0,074 K); en equilibrio, σ = 0,036 K (§6.17).
- `NOMENCLATURA.md`: commit de la v3.14.0 y paquete 025.

## [3.14.0] — 2026-10-08

- Fase: 6.3.
- Clima exportado: no cambia (el modelo es idéntico; diagnósticos, trazabilidad y documentación). La caché del modelo de 2 capas se invalida una vez porque la huella del código cubre ahora más módulos: la primera `exportar_clima.py` vuelve a simular y da el mismo clima.
- Pruebas: 119.

### Añadido
- Trazabilidad de las simulaciones largas: `procedencia.json` (versión y huella del código de cada tramo, con aviso si el código cambia a mitad), sección PROCEDENCIA en `resumen.txt` y `r["climatologia"]["procedencia"]` (`DISENO_FASE_6.3.md` §6.16).
- `cache_simulacion.MODULOS_I16`: la huella del código de una simulación con I16.
- `validar_i16.py`: chorros de invierno de cada hemisferio y estratosfera de verano; la huella del código y la procedencia en el informe. `referencias_tierra.REF_ERAI["chorros_temporada"]`, calculado con `herramientas/referencias_era_interim.py`.
- `DISENO_FASE_6.3.md` §6.16: el chorro bien medido, la lectura corregida (el desplazamiento hacia el polo de los vientos de superficie no lo explica la rejilla; hipótesis del vórtice estratosférico, Polvani y Kushner 2002) y los datos de la revisión para los pendientes (ruido del equilibrio, coste, viento sin acumular, difusión del océano, P3N sin probar con I16).
- 3 pruebas: chorro como núcleo, huella del código completa, procedencia.

### Corregido
- `validar_i16.py`: el chorro de la troposfera es un núcleo (máximo local en latitud y altura, con el viento bajando encima; Manney et al. 2011). En M3N la definición anterior medía la cola del chorro de la estratosfera en el borde de la ventana (~112 hPa). Con ERA-Interim da los mismos valores que antes.
- `validar_i16.py`: el año extra se reutiliza solo si coinciden el punto de control y la huella del código.
- `cache_simulacion.MODULOS_FISICA`: faltaban los módulos del núcleo dinámico (`fase6_*`) y los de la Fase 0 que importa `fase1_geografia.py`, y `fase1_mapa_mixto.py`.
- `DISENO_FASE_6.3.md`: el título ya no dice "BORRADOR v1"; §6.15, notas de corrección en el chorro y en "no afecta a la superficie".
- `NOMENCLATURA.md` y `README.md`: la 6.3 abarca hasta la v3.14.0 (la v3.12.1 no tiene fase); commit de la v3.13.1; paquete 024; espesor del hielo con pizarra1.
- `fase30_multicapa.py`: el comentario de `GRADIENTE_CRITICO` (con I12 el ajuste seco usa g/c_p).

## [3.13.1] — 2026-10-08

- Fase: 6.3.
- Clima exportado: no cambia (solo documentación).
- Pruebas: 116.

### Cambiado
- `DISENO_FASE_6.3.md` §6.15 completado con el informe de la 3.13.0: células de Hadley, chorros, vientos de la capa baja y tropopausa frente a ERA-Interim, y su lectura.
- `NOMENCLATURA.md`: commits de la v2.4.1–v2.4.3 y de la v3.13.0.

### Corregido
- §6.15: la alta atmósfera no tenía "causa desconocida": es la limitación declarada de la radiación gris y sin ozono (`DISENO_FASE_2.3.md` §12), con una corrección a lo que allí se esperaba (en los trópicos la estratosfera sale más caliente, no más fría).

## [3.13.0] — 2026-10-08

- Fase: 6.3.
- Clima exportado: no cambia (el modelo es idéntico; solo diagnósticos, referencias y documentación).
- Pruebas: 116.

### Añadido
- Referencias de ERA-Interim 1979–2016 (Dee et al. 2011, vía TropD) en `referencias_tierra.py`: células de Hadley, chorros, vientos de la capa baja, viento zonal por filas y tropopausa por bandas, medidas con las mismas funciones que M3N y en su rejilla. Transporte de energía a través del ecuador (Donohoe et al. 2013).
- `herramientas/referencias_era_interim.py`: el programa que calcula esas referencias, con la huella de cada archivo.
- `validar_i16.py`: tropopausa interpolada entre niveles y punto más frío; chorro de la estratosfera aparte; transporte a través del ecuador.
- `DISENO_FASE_6.3.md` §6.15: resultado de la validación en el PC de Carlos y su lectura.
- 3 pruebas en `test_validar_i16.py`.

### Corregido
- `validar_i16.py`: cada célula de Hadley se busca en su lado (en diciembre–febrero, la "del sur" era la célula de Ferrel del norte); el chorro de la troposfera ya no es el nivel más alto del modelo.
- Los 18 avisos "Mean of empty slice" de las pruebas (`nieve_permanente_posible` con menos de 12 días), con el mismo resultado.
- `NOMENCLATURA.md`: commits de todas las etiquetas; la `v2.0` nunca se etiquetó; la `v2.2` está en c67729d; las etiquetas se suben una a una (también en `herramientas/etiquetas_canonicas.sh`).

## [3.12.1] — 2026-10-08

- Fase: ninguna (nomenclatura).
- Clima exportado: no cambia (el modelo es idéntico).
- Pruebas: 113.

### Añadido
- `NOMENCLATURA.md`: la norma de fases (numeración jerárquica sin letras), versiones (SemVer), etiquetas, commits, paquetes, interruptores y documentos, con todas las correspondencias.
- `CHANGELOG.md` (este archivo), con toda la historia.
- `test_nomenclatura.py`: comprueba que el repositorio cumple la norma.
- `herramientas/etiquetas_canonicas.sh`: crea las etiquetas canónicas como alias de las antiguas.
- `herramientas/diagnostico_tierra.py`: el diagnóstico del 07/10, antes suelto en Descargas.

### Cambiado
- Documentos de diseño renombrados a `DISENO_FASE_<id>.md` (tabla en NOMENCLATURA §7), con las referencias actualizadas en todo el repositorio.
- Fases: 2b → 2.2; la antigua 8 (atmósfera de N capas) → 2.3; 5a → 5.1; 5b → 5.2; 6b → 6.5; 7 y 7b → 7.1 y 7.2; el barrido final → 9.
- Versiones citadas en documentos y comentarios con su número canónico (v3.1-pre12 → v3.10.0, etc.).
- `VERSION` sin la "v" (formato SemVer); README actualizado (plan de fases, documentos, herramientas y pruebas).

### Corregido
- La cabecera de `fase2b_atmosfera.py` hablaba de "diez interruptores"; son dieciséis.
- El factor de rotación de D aparecía como 0,684 en tres sitios; el valor es 0,683 (0,68294).

## [3.12.0] — 2026-10-08 (antes: v3.1-pre16)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 106.

### Añadido
- `validar_i16.py`: validación del modo Tierra con I16 frente a la Tierra real (temperatura CRU 1961–1990, precipitación GPCP, radiación CERES, transporte de Trenberth y Caron, hielo NSIDC, vientos y células de Hadley del modelo).
- `referencias_tierra.py`: los valores observados, cada uno con su fuente.
- `DISENO_FASE_6.3.md` §6.14: primera climatología del modo Tierra (equilibrio en el año 30, 30 años).

### Cambiado
- La versión de las cabeceras de `clima_dinamico.py` sale de una sola constante.
- Documentación: las decisiones de Carlos que constaban como abiertas (estrella hasta las nubes; rejilla de 5° por ahora) y el README al día.

## [3.11.0] — 2026-10-07 (antes: v3.1-pre15)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 101.

### Añadido
- Corrector global de energía de la dinámica (como CAM; Lauritzen y Williamson 2019), `CORRECTOR_ENERGIA_I16`. Medido en modo Tierra: la dinámica perdía −1,31 ± 0,29 W/m² en equilibrio (§6.13).

## [3.10.2] — 2026-10-07 (antes: v3.1-pre14)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 100.

### Corregido
- Inestabilidad en la meseta antártica durante la noche polar. La corrección de presión de la hiperdifusión pasa a ser una interpolación vertical sin extrapolar, con peso ≤ 0,5 (§6.12). Sustituye a la 3.10.1.

## [3.10.1] — 2026-10-07 (antes: v3.1-pre13)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 100.

### Corregido
- Inestabilidad junto a los escalones de la Antártida: corrección de presión limitada a media capa (§6.11).

## [3.10.0] — 2026-10-07 (antes: v3.1-pre12)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 99.

### Añadido
- Equilibrio (`fase6_equilibrio.py`), arranque caliente y climatología de N años (`fase6_clima.py`).
- `clima_dinamico.py`: la cadena completa, con barra de progreso y punto de control (§6.10).

## [3.9.0] — 2026-10-07 (antes: v3.1-pre11)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 95.

### Añadido
- Modelo acoplado física–núcleo: interruptor I16 `nucleo_dinamico`, apagado (§6.9).

## [3.8.0] — 2026-10-07 (antes: v3.1-pre10)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 91.

### Añadido
- Opciones del núcleo (calor de rozamiento, corrección de presión).
- `fase6_superficie.py`: rugosidad del océano del ECMWF y albedo de Cox-Munk con viento local, sin conectar (§6.8).

## [3.7.0] — 2026-10-07 (antes: v3.1-pre9)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 85.

### Añadido
- Columna con presión en superficie variable y presión de capa de Simmons y Burridge.
- `fase6_acoplamiento.py`.

### Corregido
- Caché de la geometría de la convección (§6.7).

## [3.6.0] — 2026-10-07 (antes: v3.1-pre8)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 79.

### Añadido
- Núcleo compilado con numba (`fase6_nucleo_nb.py`): idéntico bit a bit y 1,4 veces más rápido (§6.6).

## [3.5.0] — 2026-10-07 (antes: v3.1-pre7)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 76.

### Añadido
- Capa límite en columna (`fase6_capa_limite.py`: arrastre de Louis y mezcla de Frierson), sin conectar.
- Perfil del coste del núcleo (§6.5–6.6).

## [3.4.0] — 2026-10-06/07 (antes: v3.1-pre6)

- Fase: 6.3.
- Clima exportado: no cambia.
- Pruebas: 69.

### Añadido
- `held_suarez.py --niveles m3n`. Con los niveles de M3N no hace falta capa esponja (§6.4).
- Decisiones §1.1–1.8 de la 6.3, con sus fuentes (§6).

## [3.3.0] — 2026-10-06 (antes: v3.1-pre5)

- Fase: 6.2 (cierre) y 6.3 (inicio).
- Clima exportado: no cambia.
- Pruebas: 69.

### Añadido
- Cierre de la 6.2: Held y Suarez a 5° y 2,5°.
- Borrador de la 6.3 y transporte del vapor positivo y conservativo (`fase6_trazadores.py`).
- ⚠️ El reparto exacto entre las versiones 3.1.0, 3.2.0 y 3.3.0 no consta entrega a entrega en los documentos. Se comprueba con `git show --stat v3.1.0 v3.2.0 v3.3.0`.

## [3.2.0] — 2026-10-05/06 (antes: v3.1-pre4)

- Fase: 6.2.
- Clima exportado: no cambia.
- Pruebas: no consta.

### Añadido
- `held_suarez.py --filas 72`, para la prueba a 2,5° (DISENO_FASE_6.2 §10).

## [3.1.0] — 2026-10-05 (antes: v3.1-pre3; incluye los borradores v3.0-pre2, v3.1-pre1 y v3.1-pre2, nunca etiquetados)

- Fase: 2.3 (revisión), 5.1, 6.1 y 6.2 (inicio).
- Clima exportado: **sí**, por la corrección de la capa de mezcla del océano: medía 43,3 m sin querer desde la v2.4.3 y ahora ~48,7 m fijos (DISENO_FASE_2.3 §14). ⚠️ Según los documentos, la corrección entró con la revisión del 05/10 (borrador v3.0-pre2), incluida aquí.
- Pruebas: no consta.

### Añadido
- Ciclo del agua (I11–I14) y prototipo de vapor radiativo (I15), apagados (DISENO_FASE_5.1).
- Fase 6.1 (aguas someras, Williamson et al. 1992) y núcleo seco de la 6.2.

### Corregido
- Revisión del código de la 3.0.0 (DISENO_FASE_2.3 §14): la capa de mezcla, la huella de la caché, la escritura atómica y diagnósticos.

## [3.0.0] — 2026-10-05 (antes: v3.0-pre1)

- Fase: 2.3 (antes "Fase 8", adelantada).
- Clima exportado: no cambia (I10 apagado).
- Pruebas: no consta.

### Añadido
- Atmósfera de N capas en sigma (I10), modo Tierra, `calibrar_v30.py` y `validar_v30.py` (DISENO_FASE_2.3).

## [2.4.3] — 2026-10-04

- Fase: ninguna (parámetros del mundo).
- Clima exportado: **sí**.

### Cambiado
- Día solar de 19,84 h (Carlos y Ozan); hora de P3N = 1/24 del día; paso de 992 s; D de la atmósfera escalado con la rotación (factor 0,683).

## [2.4.2] — 2026-10-04

- Fase: ninguna.
- Clima exportado: no cambia (la duración del día se generaliza sin cambiar su valor).

### Cambiado
- Sin "24 h" escondidas en el código; borrador del diseño del agua (DISENO_FASE_5).

## [2.4.1] — 2026-10-03

- Fase: 4.

### Corregido
- Revisión general: hitos del año, codificación, agua a 0 m, huella del código en la caché, océano profundo y salto del hielo (DISENO_FASE_4, DISENO_FASE_3 §2.3).

## [2.4.0] — 2026-10-03 (sin etiqueta propia: incluida en la v2.4.1)

- Fase: 4.

### Añadido
- Exportación `m3n-clima` v1 y registro horario (DISENO_FASE_4).

## [2.3.0] — 2026-10-02 (antes: v2.3)

- Fase: 3.
- Clima exportado: **sí** (15,30 °C con hielo).

### Añadido
- Hielo marino termodinámico (I9) (DISENO_FASE_3).

## [2.2.2] — 2026-10-02 (antes: v2.2c; sin etiqueta propia, incluida en la v2.3)

- Fase: 2.2.
- Clima exportado: **sí**.

### Corregido
- Auditoría:
  - umbral convectivo de 40,5 K;
  - albedo del océano según el sol (I8);
  - D_atm 2,4 y D_oc 0,12;
  - masa de S3N 0,874 (DISENO_FASE_2.2 §12).

## [2.2.1] — 2026-10-02 (antes: v2.2b)

- Fase: 2.2.
- Clima exportado: **sí**.

### Añadido
- Atmósfera de 2 capas, calor sensible, convección y aire a 2 m (I1–I7) (DISENO_FASE_2.2).

## [2.2.0] — (antes: v2.2; etiqueta en c67729d: el commit a0b5059 más el arreglo del solsticio)

- Fase: 2.1.

### Añadido
- Transporte horizontal de calor (difusión) y columna de suelo.

## [2.1.0] — (antes: v2.1)

- Fase: 1.

### Añadido
- Geografía real del mapa de C3N.

## [2.0.0] — 2026-09-20 (antes: v2.0; nunca se etiquetó)

- Fase: 0.

### Añadido
- Reescritura vectorizada en una rejilla de 72 × 36 a partir del modelo original de un punto (DISENO_FASE_0).

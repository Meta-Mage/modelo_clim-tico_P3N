# Modelo Climático P3N — Proyecto B3N

Este es el modelo climático del planeta P3N del proyecto de worldbuilding B3N.
Desarrollado por Carlos (Meta-Mage) y Ozan.

---

## ¿Qué necesitas instalar?

Antes de nada, instala estas tres cosas en orden:

### 1. Python
- Ve a **python.org/downloads**
- Descarga la versión más reciente (3.12 o superior)
- Durante la instalación, **marca la casilla "Add Python to PATH"** — es importante
- Siguiente, siguiente, instalar

### 2. Git
- Ve a **git-scm.com/download/win**
- Descarga e instala (siguiente, siguiente, instalar, sin cambiar nada)

### 3. VS Code
- Ve a **code.visualstudio.com**
- Descarga e instala

---

## Descargar el proyecto

Abre la terminal (busca "cmd" en el buscador de Windows o pulsa **Windows + R** y escribe `cmd`) y ejecuta:

```
git clone https://github.com/Meta-Mage/modelo_clim-tico_P3N.git
```

Esto creará una carpeta con todos los archivos del proyecto.

---

## Instalar las librerías necesarias

En la terminal, entra en la carpeta del proyecto:

```
cd modelo_clim-tico_P3N
```

E instala las librerías:

```
pip install matplotlib numpy
```

---

## Instalar ffmpeg (necesario para los vídeos)

1. Ve a **github.com/BtbN/FFmpeg-Builds/releases**
2. Descarga `ffmpeg-master-latest-win64-gpl-shared.zip`
3. Descomprímelo en una carpeta llamada `ffmpeg` dentro de tu carpeta de usuario (por ejemplo `C:\Users\TU_NOMBRE\ffmpeg`)

Luego abre `temperatura.py` en VS Code y busca esta línea:

```python
matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\sala.AULASUC-214VNRO\ffmpeg\...'
```

Cámbiala por la ruta donde tú hayas puesto ffmpeg, por ejemplo:

```python
matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\Users\TU_NOMBRE\ffmpeg\ffmpeg-master-latest-win64-gpl-shared\bin\ffmpeg.exe'
```

---

## Abrir el proyecto en VS Code

En la terminal ejecuta:

```
code modelo_clim-tico_P3N
```

O abre VS Code, ve a **File → Open Folder** y selecciona la carpeta.

---

## Cómo usar el modelo

Abre la terminal integrada de VS Code (**Ctrl + `**) y ejecuta:

```
python temperatura.py
```

El programa te preguntará qué quieres calcular:

- **Opción 1** — Temperatura en un momento, día y lugar concreto. Genera gráficos del ciclo diario y un mapa de calor animado de 24 horas.
- **Opción 2** — Visión global anual. Muestra todas las latitudes con su temperatura media, máxima, mínima y habitabilidad. Genera un gráfico comparativo y un mapa de calor animado del año completo.

Para el barrido paramétrico (búsqueda de combinaciones habitables):

```
python main.py
```

Los resultados se guardan en `resultados.html` — ábrelo en el navegador para verlos de forma interactiva.

---

## Archivos del proyecto

| Archivo | Qué hace |
|---|---|
| `parametros.py` | Todos los parámetros físicos del sistema — aquí se cambian masa, albedo, emisividad, etc. |
| `orbita.py` | Cálculos orbitales — anomalías, distancia al Sol, estaciones |
| `geometria.py` | Geometría solar — ángulo cenital, declinación, ángulo horario |
| `atmosfera.py` | Atenuación atmosférica — masa de aire, transmitancia |
| `radiacion.py` | Irradiancia — TOA, instantánea, atenuada, absorbida |
| `temperatura.py` | Simulación de temperatura con inercia térmica — punto de entrada principal |
| `main.py` | Barrido paramétrico — busca combinaciones habitables |
| `outputs/` | Carpeta donde se guardan todos los gráficos y vídeos generados |

---

## Subir cambios a GitHub

Si modificas algo y quieres guardarlo:

```
git add .
git commit -m "Descripción de lo que cambiaste"
git push
```

## Descargar cambios que haya subido el otro

```
git pull
```

---

## Contacto

- Carlos: **Meta-Mage** en GitHub
- Cualquier duda sobre el modelo, pregúntale a Carlos

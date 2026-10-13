# fase52_ozono.py -- Fase 5.2 (BORRADOR, 13/10/2026; DISENO_FASE_5.2.md §3.3): el OZONO PRESCRITO de P3N (la opcion
# (a) del §3.3, propuesta 🔶: falta la decision de Carlos). No se usa todavia en el modelo.
#
# Propuesta: un perfil vertical de fraccion molar con la FORMA del de la Tierra, escalado a la COLUMNA que dan los
# modelos fotoquimicos para una estrella como S3N:
#   - forma: la atmosfera estandar de EE. UU. de 1976 en las tablas de Anderson et al. 1986 (AFGL-TR-86-0110, tabla 1f;
#     informe del gobierno de EE. UU., de dominio publico; leida el 13/10/2026 de la copia del paquete joseki 2.7.0):
#     maximo de 7,84 ppm a 35 km (5,7 hPa) y 345,7 DU en total (calculado aqui);
#   - columna: la de la Tierra por el cociente entre una estrella de 5500 K (G8 V, la de tau Ceti) y el Sol en
#     Rugheimer et al. 2013 (Astrobiology 13, 251, tabla 2 ✅ leida el 12/10/2026): 3,2e18 / 5,3e18 = 0,60. Con unos
#     300 DU de media en la Tierra ("averaged over the entire world the ozone column has a value of about 300 DU",
#     TEMIS/KNMI, temis.nl ✅ leido el 13/10/2026), ~180 DU.
#   Limitaciones declaradas: la forma no cambia con la estrella (Rugheimer et al.: con estrellas mas frias la
#   estratosfera es casi isoterma y la capa de ozono algo distinta); sin latitud ni estacion; el 18,9 % de O2 de P3N
#   (frente al 20,9 %) no se tiene en cuenta (el ozono depende del O2 de forma no lineal; su efecto con un 10 % menos
#   de O2 esta por buscar ⚠️).

import numpy as np

N_A = 6.02214076e23
M_AIRE = 28.964e-3                     # kg/mol (aire seco)
DU = 2.6867e16                         # moleculas/cm^2 por unidad Dobson
COCIENTE_RUGHEIMER = 3.2e18 / 5.3e18   # G8 V (5500 K) frente al Sol (Rugheimer et al. 2013, tabla 2)

# (presion hPa, fraccion molar en ppm): atmosfera estandar de EE. UU. de 1976, Anderson et al. 1986, tabla 1f
PERFIL_US = np.array([
    (1013, 0.0266), (898.8, 0.0293), (795, 0.0324),
    (701.2, 0.0332), (616.6, 0.0339), (540.5, 0.0377),
    (472.2, 0.0411), (411.1, 0.0501), (356.5, 0.0597),
    (308, 0.0917), (265, 0.131), (227, 0.215),
    (194, 0.31), (165.8, 0.385), (141.7, 0.503),
    (121.1, 0.651), (103.5, 0.87), (88.5, 1.19),
    (75.65, 1.59), (64.67, 2.03), (55.29, 2.58),
    (47.29, 3.03), (40.47, 3.65), (34.67, 4.17),
    (29.72, 4.63), (25.49, 5.12), (17.43, 5.8),
    (11.97, 6.55), (8.01, 7.37), (5.746, 7.84),
    (4.15, 7.8), (2.871, 7.3), (2.06, 6.2),
    (1.491, 5.25), (1.09, 4.1), (0.7978, 3.1),
    (0.425, 1.8), (0.219, 1.1), (0.109, 0.7),
    (0.0522, 0.3), (0.024, 0.25), (0.0105, 0.3),
    (0.00446, 0.5), (0.00184, 0.7), (0.00076, 0.7),
    (0.00032, 0.4), (0.000145, 0.2), (7.1e-05, 0.05),
    (4.01e-05, 0.005), (2.54e-05, 0.0005),
])


def forma(p_hPa):
    """Fraccion molar del perfil de la Tierra (sin escalar) a la presion p, interpolada en ln p (fuera, el extremo)."""
    lp = np.log(PERFIL_US[::-1, 0])
    return np.interp(np.log(np.asarray(p_hPa, dtype=float)), lp, PERFIL_US[::-1, 1]) * 1e-6


def columna_DU(p_seminiveles_hPa, vmr_capas, g):
    """Columna (DU) de una columna de capas: Σ vmr dp N_A / (g M) (dp en Pa)."""
    dp = np.diff(np.asarray(p_seminiveles_hPa, dtype=float)) * 100.0
    return float(np.sum(np.asarray(vmr_capas) * dp) * N_A / (g * M_AIRE) / 1e4 / DU)


def perfil(p_seminiveles_hPa, g, columna_objetivo_DU):
    """Fraccion molar de ozono de cada capa (entre seminiveles de arriba abajo) con la forma de la Tierra y la
    columna pedida (DU) para esta columna (su gravedad y su presion en superficie)."""
    ph = np.asarray(p_seminiveles_hPa, dtype=float)
    pm = np.sqrt(ph[1:] * np.maximum(ph[:-1], 1e-6))
    v = forma(pm)
    return v * columna_objetivo_DU / columna_DU(ph, v, g)

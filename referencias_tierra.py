# referencias_tierra.py -- v3.12.0: valores OBSERVADOS de la Tierra para validar el modo Tierra
# (validar_i16.py; DISENO_FASE_6.3.md §6.14). Aqui solo hay valores con la fuente comprobada.
#
# Tablas por bandas de latitud: las 36 filas de M3N (5 grados), de 87,5 N a 87,5 S (rejilla.LATITUDES_GRADOS).
#
# 1) Temperatura en superficie, climatologia CRU 1961-1990 (Jones, New, Parker, Martin y Rigor 1999, "Surface
#    air temperature and its changes over the past 150 years", Rev. Geophys. 37, 173-199,
#    doi:10.1029/1999RG900002). Archivo "absolute.nc" (5x5 grados, 12 meses, variable "tem",
#    "CRU_Global_1961-1990_Mean_Monthly_Surface_Temperature_Climatology"), copia en el repositorio publico
#    NCAR/GeoCAT-datafiles (netcdf_files/absolute.nc), sha256
#    e4babac4a86af8f4bc0bbfebf28270bfee31a56beaf7d591900bffaa7031b439.
#    Es la MISMA rejilla de 5 grados que M3N. Aire a 2 m sobre tierra; sobre el mar, temperatura del agua
#    (como en el articulo). Media global 13,97 C (norte 14,59, sur 13,36), la de Jones et al. (1999).
#    CRU_T_MES[m][i]: media zonal (C) del mes m (0 = enero) en la fila i; CRU_T_ANUAL: media de los 12 meses.
# 2) Precipitacion, GPCP version 2.2, 1979-2010 (Huffman, Adler, Bolvin y Gu 2009, GRL 36, L17808,
#    doi:10.1029/2009GL040000). Archivo "V22_GPCP.1979-2010.nc" (2,5 grados, mensual, mm/dia), copia en
#    NCAR/GeoCAT-datafiles, sha256 6d37505e655f36e9d9c889d16d142260f3eda6e57f6eacebbf54d1fcc7cc5a2d.
#    Climatologia de cada mes (media de los 32 años; julio de 1987 tiene 2 celdas sin dato y se promedian los
#    demas años), media zonal y, para cada fila de 5 grados, media de sus dos filas de 2,5 grados pesada por
#    el area. GPCP_P_ANUAL: media de los meses pesada por sus dias. Media global 2,674 mm/dia.
# 3) Valores globales: ver REF_* mas abajo, cada uno con su fuente.
#
# Como se calcularon las tablas (para repetirlo): ver DISENO_FASE_6.3.md §6.14.

CRU_T_ANUAL = (
    -17.07, -16.73, -15.14, -12.09, -8.83, -4.60, 0.48, 3.12, 6.16, 10.07, 13.62, 17.09,
    21.27, 24.03, 25.66, 26.42, 26.35, 26.15, 25.97, 25.77, 24.91, 23.77, 22.25, 20.25,
    17.79, 15.02, 11.55, 7.85, 4.31, 0.92, -3.35, -10.76, -24.68, -34.66, -38.96, -46.18,
)
CRU_T_MES = (
    (-30.33, -29.93, -27.47, -25.53, -24.49, -20.22, -11.98, -8.83, -4.87, 0.53, 5.44, 10.25, 15.52, 19.82, 22.84, 24.84, 25.82, 26.04, 26.06, 26.05, 25.70, 25.30, 24.83, 23.56, 20.96, 17.57, 13.54, 9.45, 6.02, 3.36, 1.49, -2.26, -11.48, -17.85, -20.29, -25.54),
    (-30.83, -30.21, -28.01, -25.59, -23.74, -18.84, -11.34, -7.99, -4.07, 1.39, 6.10, 10.77, 16.06, 20.31, 23.26, 25.21, 26.10, 26.33, 26.36, 26.32, 25.97, 25.64, 25.12, 23.91, 21.53, 18.08, 14.10, 9.84, 6.16, 3.67, 1.01, -4.73, -16.96, -25.80, -29.78, -37.05),
    (-30.41, -29.78, -27.07, -23.22, -19.48, -13.71, -6.65, -3.43, 0.09, 4.56, 8.58, 12.69, 17.85, 21.82, 24.40, 26.05, 26.63, 26.65, 26.66, 26.60, 26.17, 25.55, 24.72, 23.13, 20.65, 17.55, 13.65, 9.71, 6.28, 3.23, -0.82, -8.53, -23.78, -35.12, -40.43, -49.13),
    (-22.18, -22.15, -20.23, -16.34, -11.97, -6.07, -0.13, 2.84, 5.55, 8.62, 11.88, 15.40, 20.25, 23.66, 25.74, 26.95, 26.97, 26.78, 26.74, 26.68, 25.96, 24.91, 23.42, 21.43, 19.14, 16.40, 12.71, 9.00, 5.36, 1.49, -3.48, -12.15, -28.32, -39.96, -44.98, -53.37),
    (-9.21, -9.60, -9.16, -6.51, -2.48, 2.02, 5.48, 7.87, 10.03, 12.67, 15.35, 18.40, 22.85, 25.58, 26.97, 27.46, 27.01, 26.72, 26.59, 26.31, 25.17, 23.55, 21.55, 19.43, 17.37, 14.89, 11.60, 7.86, 4.12, 0.28, -4.99, -13.97, -30.16, -41.80, -46.73, -54.62),
    (-1.27, -1.54, -1.54, 0.83, 5.14, 8.17, 9.99, 11.81, 13.70, 16.23, 18.79, 21.57, 25.13, 26.76, 27.42, 27.26, 26.55, 26.16, 25.87, 25.48, 24.09, 22.23, 20.09, 17.81, 15.61, 13.71, 10.32, 7.19, 3.52, -0.48, -6.10, -15.18, -31.07, -42.46, -47.35, -55.08),
    (-0.18, 0.05, 0.67, 3.92, 8.95, 11.53, 12.71, 14.31, 16.42, 19.30, 21.74, 23.80, 26.27, 27.08, 27.23, 26.85, 26.16, 25.70, 25.32, 24.92, 23.48, 21.52, 19.28, 16.98, 14.84, 12.57, 9.27, 6.29, 2.87, -1.23, -7.02, -16.49, -32.53, -44.13, -49.20, -56.90),
    (-1.25, -1.25, -0.72, 2.36, 6.88, 9.91, 11.95, 13.98, 16.54, 19.80, 22.29, 24.22, 26.44, 27.19, 27.26, 26.80, 26.18, 25.64, 25.24, 24.85, 23.51, 21.68, 19.45, 16.99, 14.58, 12.23, 9.35, 6.10, 2.30, -1.65, -7.60, -17.12, -33.03, -44.50, -49.44, -56.81),
    (-7.92, -7.40, -5.49, -1.55, 1.87, 5.01, 8.23, 10.42, 13.16, 16.67, 19.68, 22.38, 25.45, 26.85, 27.26, 26.89, 26.23, 25.78, 25.40, 25.06, 23.95, 22.46, 20.29, 17.74, 15.04, 12.65, 9.66, 6.01, 2.58, -0.94, -6.58, -15.89, -31.75, -43.31, -48.55, -56.25),
    (-17.41, -16.95, -14.57, -10.48, -7.21, -2.88, 2.38, 4.93, 7.94, 11.83, 15.52, 18.95, 22.98, 25.47, 26.75, 26.98, 26.33, 25.93, 25.65, 25.39, 24.53, 23.36, 21.54, 19.10, 16.24, 13.51, 10.36, 6.59, 3.08, -0.17, -4.65, -12.47, -26.63, -36.61, -40.90, -48.06),
    (-24.68, -23.96, -22.12, -19.36, -17.24, -12.31, -4.85, -1.70, 2.04, 6.76, 10.98, 14.97, 19.67, 23.08, 25.25, 26.36, 26.30, 26.07, 25.84, 25.70, 25.02, 24.17, 22.81, 20.60, 17.83, 14.77, 11.35, 7.45, 4.08, 1.14, -1.87, -7.24, -18.48, -26.29, -29.50, -35.70),
    (-29.23, -28.07, -26.03, -23.62, -22.20, -17.83, -10.04, -6.73, -2.60, 2.45, 7.11, 11.70, 16.74, 20.76, 23.56, 25.35, 25.88, 26.00, 25.89, 25.83, 25.35, 24.81, 23.91, 22.28, 19.63, 16.28, 12.65, 8.72, 5.38, 2.37, 0.47, -3.10, -11.96, -18.05, -20.42, -25.59),
)
GPCP_P_ANUAL = (
    0.446, 0.491, 0.671, 1.036, 1.573, 2.009, 2.404, 2.438, 2.507, 2.727, 2.729, 2.328,
    1.908, 1.741, 2.079, 3.301, 5.602, 4.560, 3.993, 3.788, 2.813, 2.101, 1.887, 2.146,
    2.474, 2.810, 2.933, 2.815, 2.944, 3.230, 2.227, 1.312, 0.818, 0.586, 0.359, 0.376,
)
GPCP_P_MES = (
    (0.329, 0.367, 0.554, 1.008, 1.592, 2.012, 2.517, 2.530, 2.691, 3.148, 3.197, 2.496, 1.554, 0.954, 0.895, 1.237, 3.317, 4.541, 4.719, 5.482, 5.014, 3.926, 2.747, 2.341, 2.308, 2.370, 2.499, 2.779, 2.699, 2.821, 2.022, 1.321, 0.729, 0.431, 0.183, 0.122),
    (0.257, 0.335, 0.513, 0.909, 1.441, 1.808, 2.169, 2.238, 2.394, 2.831, 3.143, 2.528, 1.554, 0.963, 0.777, 0.948, 2.586, 4.272, 4.838, 5.608, 5.203, 4.097, 2.933, 2.504, 2.409, 2.437, 2.594, 2.900, 2.875, 3.128, 2.305, 1.565, 0.936, 0.530, 0.240, 0.230),
    (0.179, 0.248, 0.437, 0.741, 1.212, 1.560, 1.904, 2.006, 2.154, 2.567, 2.888, 2.499, 1.596, 0.934, 0.740, 0.884, 2.895, 4.910, 5.329, 5.483, 4.448, 3.314, 2.608, 2.563, 2.626, 2.721, 2.899, 3.069, 3.687, 4.189, 2.786, 1.640, 0.928, 0.564, 0.286, 0.397),
    (0.228, 0.274, 0.417, 0.683, 1.082, 1.345, 1.568, 1.735, 2.019, 2.358, 2.593, 2.212, 1.717, 1.119, 0.850, 1.222, 4.204, 5.869, 5.514, 4.650, 3.098, 2.309, 2.098, 2.427, 2.670, 2.982, 3.260, 3.288, 3.685, 3.925, 2.584, 1.440, 0.797, 0.550, 0.323, 0.496),
    (0.253, 0.259, 0.362, 0.649, 1.081, 1.407, 1.710, 1.878, 2.124, 2.321, 2.432, 2.191, 1.992, 1.740, 1.601, 2.515, 6.292, 6.037, 4.347, 3.273, 2.042, 1.571, 1.808, 2.330, 2.799, 3.244, 3.556, 3.162, 3.206, 3.540, 2.351, 1.247, 0.780, 0.609, 0.389, 0.455),
    (0.396, 0.443, 0.585, 0.833, 1.349, 1.848, 2.288, 2.411, 2.411, 2.437, 2.454, 2.494, 2.351, 2.186, 2.570, 4.223, 7.633, 5.057, 3.279, 2.588, 1.463, 1.100, 1.387, 2.043, 2.716, 3.349, 3.605, 3.037, 2.975, 3.238, 2.214, 1.193, 0.812, 0.690, 0.468, 0.415),
    (0.723, 0.756, 0.894, 1.156, 1.735, 2.190, 2.644, 2.734, 2.574, 2.435, 2.301, 2.260, 2.350, 2.760, 3.326, 5.433, 7.489, 4.076, 2.961, 2.371, 1.262, 0.880, 1.183, 1.896, 2.588, 3.266, 3.301, 2.854, 2.937, 3.194, 2.107, 1.215, 0.869, 0.706, 0.498, 0.493),
    (0.984, 0.967, 1.125, 1.451, 2.053, 2.540, 2.758, 2.761, 2.491, 2.395, 2.116, 2.163, 2.470, 3.091, 4.039, 6.491, 7.306, 3.622, 2.757, 2.138, 1.086, 0.792, 1.028, 1.707, 2.473, 3.087, 3.114, 2.787, 2.950, 3.199, 2.093, 1.244, 0.895, 0.709, 0.519, 0.527),
    (0.879, 0.859, 1.027, 1.340, 1.907, 2.463, 2.925, 2.722, 2.595, 2.667, 2.409, 2.215, 2.289, 2.792, 3.826, 5.928, 7.024, 3.612, 2.877, 2.415, 1.318, 0.884, 1.118, 1.738, 2.343, 2.869, 2.878, 2.617, 2.814, 3.183, 2.188, 1.232, 0.858, 0.736, 0.551, 0.583),
    (0.450, 0.547, 0.854, 1.405, 1.986, 2.400, 2.810, 2.671, 2.762, 2.928, 2.815, 2.203, 1.876, 1.865, 2.888, 4.846, 6.837, 4.024, 3.223, 3.042, 2.013, 1.350, 1.469, 1.905, 2.255, 2.585, 2.601, 2.429, 2.674, 3.065, 2.217, 1.251, 0.781, 0.596, 0.403, 0.440),
    (0.344, 0.456, 0.674, 1.148, 1.775, 2.373, 2.878, 2.830, 2.953, 3.278, 3.094, 2.320, 1.644, 1.379, 2.033, 3.581, 6.479, 4.226, 3.802, 3.878, 2.930, 2.040, 1.916, 2.069, 2.273, 2.445, 2.455, 2.299, 2.492, 2.852, 2.135, 1.285, 0.798, 0.528, 0.278, 0.241),
    (0.316, 0.374, 0.592, 1.092, 1.649, 2.144, 2.657, 2.725, 2.899, 3.357, 3.324, 2.373, 1.486, 1.057, 1.325, 2.151, 4.991, 4.460, 4.325, 4.639, 4.010, 3.059, 2.407, 2.252, 2.229, 2.345, 2.419, 2.565, 2.338, 2.422, 1.733, 1.127, 0.641, 0.386, 0.163, 0.109),
)

# Calendario del modo Tierra: el dia 0 de la orbita de M3N es el solsticio de diciembre (declinacion -23,44
# grados; el perihelio cae 12,7 dias despues, el 3-4 de enero). Asi, el 1 de enero es el dia 11 (+-1 dia).
DIA_1_ENERO = 11
DIAS_MES = (31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Temperatura media global
REF_T_GLOBAL = {
    "CRU 1961-1990 (Jones et al. 1999)": 13.97,
    # Copernicus/ECMWF (ERA5): 2024 fue 15,10 C, 0,72 C por encima de 1991-2020 y 1,60 C sobre 1850-1900
    # (https://climate.copernicus.eu/copernicus-2024-first-year-exceed-15degc-above-pre-industrial-level)
    "ERA5 1991-2020 (Copernicus)": 15.10 - 0.72,
}
# Precipitacion media global, GPCP v2.3 1979-2016 (Adler et al. 2018, Atmosphere 9, 138, doi:10.3390/atmos9040138):
# 2,69 mm/dia (tierra 2,24; oceano 2,90; los tres ya usados en validar_v30.py)
REF_P_GLOBAL = 2.69
# Agua precipitable media: 24,9 kg/m2 (Trenberth y Smith 2005: 1,27e16 kg de vapor)
REF_AGUA_PRECIPITABLE = 24.9

# Balance de radiacion en el tope, CERES EBAF Ed4.0, julio 2005 - junio 2015, con la restriccion del
# desequilibrio neto (Loeb et al. 2018, J. Climate 31, 895-918, doi:10.1175/JCLI-D-17-0208.1, tabla 5). W/m2.
REF_CERES = {
    "entrada": 340.0,
    "reflejada_todo_cielo": 99.1, "reflejada_despejado": 53.3,
    "olr_todo_cielo": 240.1, "olr_despejado": 268.1,
    "neto_todo_cielo": 0.7,
}
# Superficie con cielo despejado (Wild et al. 2019, Clim. Dyn. 52; los mismos que validar_v30.py): luz absorbida
# en el suelo 214 y en el aire 73; infrarrojo hacia el suelo 314 (con la atmosfera real, calentada tambien por
# las nubes).
REF_WILD = {"sw_suelo": 214, "sw_atm": 73, "dlr": 314}
# Transporte de calor hacia los polos (Trenberth y Caron 2001, J. Climate 14, 3433): atmosfera maximo
# 5,0 +- 0,14 PW a 43 N, parecido cerca de 40 S; a 35 grados el oceano lleva el 22 % (norte) y el 8 % (sur).
REF_TRANSPORTE = {"atm_max_PW": 5.0, "atm_max_lat": 43.0, "oceano_35N": 0.22, "oceano_35S": 0.08}
# Hielo marino, extension (millones de km2), NSIDC 1981-2010 (DISENO_FASE_3.md §6.1; validar_v30.py)
REF_HIELO = {"N_max": 15.6, "N_min": 6, "S_max": 18.5, "S_min": 3}

# 3.13.0 -- VIENTOS, CELULAS DE HADLEY Y TROPOPAUSA: ERA-Interim (Dee et al. 2011, Q. J. R. Meteorol. Soc. 137,
# 553-597), medias zonales mensuales de 1979 a 2016 distribuidas con TropD (Adam et al. 2018, Geosci. Model Dev.
# 11, 4339; github.com/tropd/pytropd, ValidationData, commit c25ce1118c748f678869b9e69d19456a73fab1f4; huellas
# sha256 en herramientas/referencias_era_interim.py, que es el programa que calcula estos valores).
# Se miden con las MISMAS funciones que M3N (validar_i16.celulas_hadley, chorros, vientos_bajos, perfil_banda),
# despues de pasar ERA-Interim a la rejilla de 5 grados de M3N.
#   hadley[temporada][celula] = (maximo de psi en 10^10 kg/s, latitud de la cara, borde a 500 hPa)
#   chorros["capa_hemisferio"] = (m/s, latitud, p / 1000 hPa), de la media anual. 3.14.0: el de la troposfera
#       es el NUCLEO (maximo local en latitud y altura, validar_i16.chorros); en ERA-Interim sale igual que con la
#       definicion de la 3.13.0 (comprobado: mismos valores), porque su chorro subtropical es un maximo claro.
#   chorros_temporada[temporada]["capa_hemisferio"] (3.14.0): lo mismo en diciembre-febrero y junio-agosto. En la
#       estratosfera de verano el "maximo" es negativo: alli el viento es del este (calentamiento del ozono).
#   vientos_bajos_975hPa (comparable con la capa baja de M3N, a ~200 m) y vientos_bajos_10m = (m/s, latitud)
#   u_bandas[hPa] = viento zonal medio anual en cada fila de M3N, de 87,5 N a 87,5 S. A 975 hPa, por debajo
#       de ~65 S (Antartida) y en parte sobre las montañas, ese nivel queda bajo el suelo y el valor es extrapolado.
#   tropopausa[banda]: la de la OMM en un nivel y la interpolada, el punto mas frio entre 500 y 50 hPa y T a 25 km
#       (C). Comprobacion independiente: Seidel et al. (2001, JGR 106, 7857), con radiosondas 1961-1990, dan en el
#       ecuador el punto mas frio a ~16,9 km y ~96 hPa, a unos -81 C, y la tropopausa de la OMM a ~16,5 km.
REF_ERAI = {'fuente': 'ERA-Interim (Dee et al. 2011), medias zonales mensuales 1979-2016 de TropD '
           '(github.com/tropd/pytropd, ValidationData, commit c25ce11), en la rejilla de 5 grados de M3N',
 'hadley': {'anual': {'norte': (8.37, 15.0, 32.49), 'sur': (-11.48, -10.0, -32.14)},
            'dic-feb': {'norte': (20.8, 10.0, 29.99), 'sur': (-5.45, -20.0, -36.27)},
            'jun-ago': {'norte': (2.53, 25.0, 37.9), 'sur': (-23.91, -5.0, -28.97)}},
 'chorros': {'troposfera_norte': (27.08, 32.5, 0.175),
             'troposfera_sur': (30.76, -32.5, 0.2),
             'estratosfera_norte': (16.53, 57.5, 0.001),
             'estratosfera_sur': (33.7, -62.5, 0.005)},
 'chorros_temporada': {'dic-feb': {'troposfera_norte': (42.13, 32.5, 0.2),
                                   'troposfera_sur': (31.4, -47.5, 0.225),
                                   'estratosfera_norte': (45.58, 47.5, 0.001),
                                   'estratosfera_sur': (-6.08, -67.5, 0.01)},
                       'jun-ago': {'troposfera_norte': (21.6, 42.5, 0.2),
                                   'troposfera_sur': (40.98, -27.5, 0.2),
                                   'estratosfera_norte': (-7.07, 67.5, 0.01),
                                   'estratosfera_sur': (88.9, -47.5, 0.001)}},
 'vientos_bajos_975hPa': {'alisios_n': (-3.76, 17.5),
                          'alisios_s': (-4.73, -17.5),
                          'oeste_n': (2.33, 47.5),
                          'oeste_s': (9.18, -52.5)},
 'vientos_bajos_10m': {'alisios_n': (-3.14, 17.5),
                       'alisios_s': (-3.98, -17.5),
                       'oeste_n': (1.76, 47.5),
                       'oeste_s': (7.02, -52.5)},
 'u_bandas': {250: [1.84, 3.93, 5.89, 7.64, 9.01, 10.65, 13.22, 16.71, 20.53, 23.5, 24.78, 24.16, 21.31, 16.13,
                    9.92, 3.75, -1.19, -3.04, -2.16, 0.32, 5.08, 12.64, 21.31, 26.94, 27.4, 26.52, 27.37, 28.71,
                    28.04, 24.21, 18.16, 12.25, 7.97, 5.23, 2.97, 0.92],
              500: [1.2, 2.56, 3.68, 4.52, 5.2, 6.07, 7.75, 9.88, 12.14, 13.49, 13.01, 11.51, 9.21, 5.5, 1.06,
                    -2.24, -3.66, -3.84, -3.44, -2.86, -0.93, 3.04, 7.61, 10.83, 13.03, 15.27, 17.85, 19.79,
                    19.67, 16.91, 11.89, 6.42, 2.68, 1.8, 1.11, -0.51],
              975: [0.12, 0.35, 0.15, -0.13, 0.21, 0.5, 1.35, 2.05, 2.33, 2.27, 1.85, 0.81, -1.2, -3.0, -3.76,
                    -3.58, -2.05, -1.72, -2.54, -3.49, -4.54, -4.73, -3.87, -2.05, 0.67, 3.89, 6.77, 8.74, 9.18,
                    7.61, 3.91, -0.85, -2.26, -0.89, -0.7, -2.0]},
 'tropopausa': {'90..60': {'tropopausa_km': 8.78,
                           'T_tropopausa': -52.07,
                           'tropopausa_interp_km': 9.33,
                           'T_tropopausa_interp': -53.14,
                           'frio_km': 20.29,
                           'frio_hPa': 50.0,
                           'T_frio': -55.24,
                           'T_25km': -54.61},
                '60..30': {'tropopausa_km': 11.88,
                           'T_tropopausa': -54.72,
                           'tropopausa_interp_km': 11.68,
                           'T_tropopausa_interp': -54.28,
                           'frio_km': 16.25,
                           'frio_hPa': 100.0,
                           'T_frio': -60.24,
                           'T_25km': -53.5},
                '30..10': {'tropopausa_km': 16.55,
                           'T_tropopausa': -76.6,
                           'tropopausa_interp_km': 16.36,
                           'T_tropopausa_interp': -75.97,
                           'frio_km': 16.55,
                           'frio_hPa': 100.0,
                           'T_frio': -76.6,
                           'T_25km': -54.13},
                '10..-10': {'tropopausa_km': 16.56,
                            'T_tropopausa': -80.82,
                            'tropopausa_interp_km': 16.48,
                            'T_tropopausa_interp': -80.43,
                            'frio_km': 16.56,
                            'frio_hPa': 100.0,
                            'T_frio': -80.82,
                            'T_25km': -54.79},
                '-10..-30': {'tropopausa_km': 16.53,
                             'T_tropopausa': -76.3,
                             'tropopausa_interp_km': 16.31,
                             'T_tropopausa_interp': -75.56,
                             'frio_km': 16.53,
                             'frio_hPa': 100.0,
                             'T_frio': -76.3,
                             'T_25km': -54.11},
                '-30..-60': {'tropopausa_km': 11.0,
                             'T_tropopausa': -54.05,
                             'tropopausa_interp_km': 11.3,
                             'T_tropopausa_interp': -54.6,
                             'frio_km': 16.12,
                             'frio_hPa': 100.0,
                             'T_frio': -59.95,
                             'T_25km': -53.7},
                '-60..-90': {'tropopausa_km': 9.58,
                             'T_tropopausa': -58.83,
                             'tropopausa_interp_km': 9.11,
                             'T_tropopausa_interp': -57.78,
                             'frio_km': 19.63,
                             'frio_hPa': 50.0,
                             'T_frio': -61.0,
                             'T_25km': -57.83}}}

# Transporte de energia de la atmosfera a traves del ecuador y la banda de lluvias (Donohoe et al. 2013,
# J. Climate 26, 3597): media anual observada 0,1 PW hacia el SUR; la banda de lluvias (su centroide) esta en
# 1,65 N; cuanto mas transporta la atmosfera hacia el norte, mas al sur esta la banda (-2,4 a -3,2 grados por PW).
REF_TRANSPORTE_ECUADOR = {"atm_PW": -0.1, "itcz_centroide": 1.65}

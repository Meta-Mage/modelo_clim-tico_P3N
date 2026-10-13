# fase52_continuo.py -- Fase 5.2 (BORRADOR, 12/10/2026; DISENO_FASE_5.2.md §7): el CONTINUO del vapor de agua,
# MT_CKD 4.3 (Mlawer et al. 2012, Phil. Trans. R. Soc. A 370, 2520, doi:10.1098/rsta.2011.0295), como DATO (decision
# de Carlos del 12/10/2026: MT_CKD como dato). Es la absorcion del vapor que no explican las lineas cortadas a 25 cm^-1
# con el pedestal restado (T4); se SUMA a ellas. No se usa todavia en el modelo.
#
# Datos: absco-ref_wv-mt-ckd.nc de MT_CKD 4.3 (AER; github.com/AER-RC/MT_CKD, commit 679f0a4; sha256 69944eb8...;
# lo descarga herramientas/descargar_hitran.py). Licencia de AER: uso cientifico y de investigacion; se puede
# redistribuir con su aviso de copyright; no se puede incorporar a software comercial. Por eso no va en el repositorio:
# se descarga, como HITRAN.
#
# Formulas: las del modulo oficial (LBLRTM/src/mt_ckd_h2o_module.f90, subrutina mt_ckd_h2o_absco, leida el
# 12/10/2026), traducidas sin cambios:
#   rho_rat = (p / p_ref) (T_ref / T)                              (p_ref = 1013 mb, T_ref = 296 K, del archivo)
#   self    = self_absco_ref(nu) (T_ref/T)^self_texp(nu) * vmr       * rho_rat * R(nu, T)
#   foreign = for_absco_ref(nu)                          * (1 - vmr) * rho_rat * R(nu, T)
#   R(nu, T) = nu (1 - e^-x)/(1 + e^-x), x = nu/(T/radcn2)  (y 0,5 x nu si x <= 0,01; nu si x > 10): "radiation term"
#   coeficientes en cm^2 por molecula de AGUA; se interpolan de su rejilla de 10 cm^-1 con el polinomio de 4 puntos de
#   la subrutina myxint del mismo modulo.
# radcn2 = 1,4387752 cm K (la de LBLRTM, para reproducir MT_CKD exactamente; la de CODATA 2018 es 1,438776877: difiere
# en 1e-6).

import numpy as np

RADCN2 = 1.4387752


def leer_mtckd(ruta):
    from scipy.io import netcdf_file
    with netcdf_file(ruta, "r", mmap=False) as f:
        d = {k: np.array(f.variables[k].data, dtype=float) for k in
             ("wavenumbers", "self_absco_ref", "for_absco_ref", "self_texp")}
        d["ref_press"] = float(f.variables["ref_press"].data)
        d["ref_temp"] = float(f.variables["ref_temp"].data)
    return d


def termino_radiacion(nu, T):
    """El "radiation term" de MT_CKD (funcion myradfn)."""
    nu = np.asarray(nu, dtype=float)
    x = nu / (T / RADCN2)
    out = nu.copy()
    peq = x <= 0.01
    med = (x > 0.01) & (x <= 10.0)
    out[peq] = 0.5 * x[peq] * nu[peq]
    e = np.exp(-x[med])
    out[med] = nu[med] * (1.0 - e) / (1.0 + e)
    return out


def _interpolar(nu, rejilla, a):
    """Interpolacion de 4 puntos de MT_CKD (myxint): de la rejilla uniforme de los datos a los nu pedidos."""
    nu = np.asarray(nu, dtype=float)
    dva = rejilla[1] - rejilla[0]
    j = np.floor((nu - rejilla[0]) / dva).astype(int)          # el punto j tiene nu_j <= nu < nu_{j+1}
    j = np.clip(j, 1, len(rejilla) - 3)
    p = (nu - rejilla[j]) / dva
    c = (3.0 - 2.0 * p) * p * p
    b = 0.5 * p * (1.0 - p)
    b1, b2 = b * (1.0 - p), b * p
    return -a[j - 1] * b1 + a[j] * (1.0 - c + b2) + a[j + 1] * (c + b1) - a[j + 2] * b2


def coeficientes(nu, T, p_hPa, vmr, datos, radiacion=True):
    """(propio, ajeno) en cm^2 por molecula de agua en los numeros de onda nu (cm^-1)."""
    rho_rat = (p_hPa / datos["ref_press"]) * (datos["ref_temp"] / T)
    w = datos["wavenumbers"]
    s = datos["self_absco_ref"] * (datos["ref_temp"] / T) ** datos["self_texp"] * vmr * rho_rat
    f = datos["for_absco_ref"] * (1.0 - vmr) * rho_rat
    if radiacion:
        R = termino_radiacion(w, T)
        s, f = s * R, f * R
    return _interpolar(nu, w, s), _interpolar(nu, w, f)


def k_continuo(nu, T, p_hPa, vmr, datos):
    """k del continuo del vapor (propio + ajeno), en cm^2 por molecula de agua; 0 fuera de los datos (0-20 000 cm^-1)."""
    s, f = coeficientes(nu, T, p_hPa, vmr, datos)
    nu = np.asarray(nu, dtype=float)
    fuera = (nu < datos["wavenumbers"][1]) | (nu > datos["wavenumbers"][-3])
    k = s + f
    k[fuera] = 0.0
    return np.maximum(k, 0.0)

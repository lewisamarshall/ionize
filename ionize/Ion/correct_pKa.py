import warnings
from math import log


def correct_pKa(obj):
    """Return the pKa corrected for temperature.

    If dCp for the ion is available, returns the Clark-Glew correction, which
    is most accurate. If only dH is available, returns the van't Hoff
    correction, which is less accurate. If neither is available, returns the
    original pKa, with a warning.
    """
    if obj.T == obj._T_ref:
        return obj._pKa_ref
    elif obj.dH and obj.dCp:
        return clark_glew(obj)
    elif obj.dH and not obj.dCp:
        return vant_hoff(obj)
    else:
        warnings.warn('No data available to correct pKa for temperature.')
        return obj._pKa_ref


def vant_hoff(obj):
    T = obj.T + 273.15
    T_ref = obj._T_ref + 273.15
    if abs(T-T_ref) > 20:
        warnings.warn('Using the van\'t Hoff correction for dT > 20 deg.')
    pKa_ref = obj._pKa_ref
    dH = obj.dH
    if dH and len(dH) == len(pKa_ref):
        pKa = [p - h/(2.303 * obj._R)*(1/T_ref - 1/T)
               for p, h in zip(pKa_ref, dH)]
    else:
        warnings.warn('No dH available. Returning uncorrected pKa.')
    return pKa


def clark_glew(obj):
    T = obj.T + 273.15
    T_ref = obj._T_ref + 273.15
    if abs(T-T_ref) > 100:
        warnings.warn('Using the Clark-Glew correction for dT > 100 deg.')
    pKa_ref = obj._pKa_ref
    dH = obj.dH
    dCp = obj.dCp
    if dH and dCp and len(dH) == len(pKa_ref) == len(dCp):
        pKa = [p - h/(2.303 * obj._R)*(1/T_ref - 1/T)
               - c/(2.303 * obj._R) * (T_ref/T - 1 - log(T/T_ref))
               for p, h, c in zip(pKa_ref, dH, dCp)]
    else:
        warnings.warn('No dCp available. Returning uncorrected pKa.')
    return pKa

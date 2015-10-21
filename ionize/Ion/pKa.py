import warnings
from math import log
from ..constants import gas_constant, kelvin
import numpy as np


def acidity(self, ionic_strength=None, temperature=None):
    """Return the effective Ka values for the ion.

    Args:
        I (float): The ambiant ionic strength.

    This function correct the Ka for ionic strength, using the Dubye-Huckel
    theory to calculate activity coefficients. If no ionic strength is
    supplied, and the Ion is nested in a Solution, the solution ionic
    strength will be used. Otherwise, the ionic strength is assumed to be 0.
    """
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    # Make the effective Ka vector the same size as the Ka vector.
    Ka_eff = []

    gam_i = self.activity(None, ionic_strength, temperature)
    gam_h, = self.activity(1, ionic_strength, temperature)

    # For each acidity coefficient, get the effective
    # coefficient by multiplying by activities.
    for i, Kp in enumerate(self.mid_Ka(ionic_strength, temperature)):
        Ka_eff.append(Kp*gam_i[i+1]/gam_i[i]/gam_h)

    return np.array(Ka_eff)


def pKa(self, ionic_strength=None, temperature=None):
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    return -np.log10(self.acidity())


def mid_Ka(self, ionic_strength, temperature):
    return 10**(-self.mid_pKa(ionic_strength, temperature))


def mid_pKa(self, ionic_strength=None, temperature=None):
    """Return the pKa corrected for temperature.

    If dCp for the ion is available, returns the Clark-Glew correction, which
    is most accurate. If only dH is available, returns the van't Hoff
    correction, which is less accurate. If neither is available, returns the
    original pKa, with a warning.
    """
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    if self.enthalpy is not None and self.heat_capacity is not None:
        return _clark_glew(self, temperature)
    elif self.enthalpy is not None and self.heat_capacity is None:
        return _vant_hoff(self, temperature)
    else:
        warnings.warn('No data available to correct pKa for temperature.')
        return self.reference_pKa


def _vant_hoff(self, temperature):
    temperature = kelvin(temperature)
    reference_temperature = kelvin(self.reference_temperature)

    if abs(temperature - reference_temperature) > 20:
        warnings.warn("Using the van't Hoff correction for dT > 20 deg.")

    if len(self.enthalpy) != len(self.reference_pKa):
        raise RuntimeError('Enthalpy must have an entry for each pKa.'
                           + repr(self))

    pKa_ref = self.reference_pKa
    dH = self.enthalpy
    pKa = [p - h/(2.303 * gas_constant)*(1/reference_temperature -
                                         1/temperature)
           for p, h in zip(pKa_ref, dH)]
    return np.array(pKa)


def _clark_glew(self, temperature):
    T = kelvin(temperature)
    T_ref = kelvin(self.reference_temperature)
    if abs(T-T_ref) > 100:
        warnings.warn('Using the Clark-Glew correction for dT > 100 deg.')
    pKa_ref = self.reference_pKa
    dH = self.enthalpy
    dCp = self.heat_capacity
    if dH is not None and \
       dCp is not None and \
       len(dH) == len(pKa_ref) == len(dCp):
        pKa = [p - h/(2.303 * gas_constant)*(1/T_ref - 1/T) -
               c/(2.303 * gas_constant) * (T_ref/T - 1 - log(T/T_ref))
               for p, h, c in zip(pKa_ref, dH, dCp)]
    else:
        warnings.warn('No dCp available. Returning uncorrected pKa.')
        pKa = pKa_ref
    return np.array(pKa)

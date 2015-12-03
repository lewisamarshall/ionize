from __future__ import division
import warnings
from math import log, sqrt
import numpy as np

from ..constants import gas_constant, kelvin, pitts


def acidity(self, ionic_strength=None, temperature=None):
    """Return the effective acidity constant, Ka, for each valence state.

    The value is corrected for ionic strength and temperature using the
    Debye-Huckel approximation.
    """
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    if self.enthalpy is not None and self.heat_capacity is not None:
        acidity = self._clark_glew_acidity(temperature)
    elif self.enthalpy is not None and self.heat_capacity is None:
        acidity = self._vant_hoff_acidity(temperature)
    else:
        if temperature != self.reference_temperature:
            warnings.warn('No data available to correct pKa for temperature.')
        acidity = 10**(-self.reference_pKa)

    # Correct for the activity of ion and H+
    gam_i = self._solvent.activity(self._valence_zero(), ionic_strength, temperature)
    gam_h = self._solvent.activity(1, ionic_strength, temperature)
    acidity = acidity * gam_i[1:] / gam_i[:-1] / gam_h

    return acidity


def pKa(self, ionic_strength=None, temperature=None):
    """Return the effective pKa for the ion."""
    return -np.log10(self.acidity(ionic_strength, temperature))


def _vant_hoff_pKa(self, temperature):
    temperature = kelvin(temperature)
    reference_temperature = kelvin(self.reference_temperature)

    if abs(temperature - reference_temperature) > 20:
        warnings.warn("Using the van't Hoff correction for dT > 20 deg.")

    pKa = (self.reference_pKa -
           self.enthalpy / (2.303 * gas_constant) *
           (1/reference_temperature - 1/temperature))
    return pKa


def _vant_hoff_acidity(self, temperature):
    return 10.**(-self._vant_hoff_pKa(temperature))


def _clark_glew_pKa(self, temperature):
    T = kelvin(temperature)
    T_ref = kelvin(self.reference_temperature)

    if abs(T-T_ref) > 100:
        warnings.warn('Using the Clark-Glew correction for dT > 100 deg.')

    pKa = (self._vant_hoff_pKa(temperature) -
           self.heat_capacity /
           (2.303 * gas_constant) * (T_ref/T - 1 - log(T/T_ref))
           )

    return pKa


def _clark_glew_acidity(self, temperature):
    return 10.**(-self._clark_glew_pKa(temperature))

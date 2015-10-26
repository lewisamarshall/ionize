from math import copysign, sqrt
import warnings
import numpy as np

from ..constants import pitts

def mobility(self, pH=None, ionic_strength=None, temperature=None):
    """Return the effective mobility of the ion at a given pH and I.

    Args:
        pH (float): The ambiant pH.

        I (float): The ambiant ionic strength.

    If an actual mobility from Onsager-Fouss is available, it is used,
    otherwise, the Robinson-Stokes mobility estimate is used.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].effective_mobility()

    Otherwise, always call with a pH argument.
    """
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    ionization_fraction = self.ionization_fraction(pH, ionic_strength,
                                                   temperature)
    actual_mobility = self.actual_mobility(ionic_strength, temperature)

    effective_mobility = np.sum(ionization_fraction * actual_mobility)

    return effective_mobility

def actual_mobility(self, ionic_strength=None, temperature=None):
    if ionic_strength is None and \
            temperature is None and \
            self.context() is not None:
        try:
            return self.context().actual_mobility(self)
        except:
            warnings.warn('Context failed to return an actual mobility.')
    return self.robinson_stokes_mobility(ionic_strength, temperature)

def absolute_mobility(self, temperature):
    _, _, temperature = \
        self._resolve_context(None, None, temperature)

    if self._nightingale_function:
        absolute_mobility = \
            (self._nightingale_function(temperature).tolist() *
             10.35e-11 / self._solvent.viscosity(temperature) *
             self.valence)

        if not (self.nightingale_data['min'] <
                temperature <
                self.nightingale_data['max']):
            warnings.warn('Temperature outside range'
                          'for nightingale data.')
    else:
        absolute_mobility =\
            (self._solvent.viscosity(self.reference_temperature) /
             self._solvent.viscosity(temperature)*self.reference_mobility)

    return absolute_mobility

def robinson_stokes_mobility(self, ionic_strength=None, temperature=None):
    """Return the Robinson-Stokes correction to fully ionized mobility.

    This correction is appropriate if a generic ionic strength is known,
    but the specific ions in solution are unknown.
    """
    # Currently using the ionic strength where Bahga 2010
    # uses twice the ionic strength. This appears to work, and follows the
    # SPRESSO implimentation.
    # Likely typo in paper.
    # TODO: Recheck this typo.
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    reference_temperature = self.reference_temperature
    d = self._solvent.dielectric(temperature)
    d_ref = self._solvent.dielectric(reference_temperature)

    A = 0.2297*((reference_temperature+273.15)*d_ref/(temperature+273.15)/d)**(-1.5)
    B = 31.410e-9 * ((reference_temperature+273.15)*d_ref/(temperature+273.15)/d)**(-0.5) *\
        self._solvent.viscosity(reference_temperature)/self._solvent.viscosity(temperature)

    absolute_mobility = self.absolute_mobility(temperature)
    actual_mobility = (absolute_mobility -
                       (A * absolute_mobility +
                        B * np.sign(self.valence)) *
                       sqrt(ionic_strength) /
                       (1. + pitts * sqrt(ionic_strength))
                       )

    return actual_mobility

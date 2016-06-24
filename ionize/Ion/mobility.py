from __future__ import division
from math import copysign, sqrt
import warnings
import numpy as np

from ..constants import pitts, reference_temperature, kelvin, faraday, \
    elementary_charge, lpm3, gpkg, onsager_fuoss

def mobility(self, pH=None, ionic_strength=None, temperature=None):
    """Return the effective mobility of the ion in m^2/V/s.

    If a context solution is available, mobility uses the full Onsager-Fuoss
    correction to mobility. Otherwise, the Robinson-Stokes model is used.

    :param pH
    :param ionic_strength
    :param temperature
    """
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    ionization_fraction = self.ionization_fraction(pH, ionic_strength,
                                                   temperature)
    actual_mobility = self.actual_mobility(ionic_strength, temperature)

    effective_mobility = np.sum(ionization_fraction * actual_mobility)

    return effective_mobility


def actual_mobility(self, ionic_strength=None, temperature=None):
    """Return the mobility for each charge state.

    Uses the Onsager-Fuoss correction if a context solution is available,
    otherwise uses the Robinson-Stokes mobility.

    :param ionic_strength
    :param temperature
    """
    try:
        return self.onsager_fuoss_mobility()
    except AttributeError:
        # warnings.warn("Insufficient information for Onsager-Fuoss "
        #               "correction to mobility. Returning the Robinson-"
        #               "Stokes approximation.")
        return self.robinson_stokes_mobility(ionic_strength, temperature)


def absolute_mobility(self, temperature=None):
    """Return the mobility of each charge state at infinite dilution.

    :param temperature
    """
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
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    dielectric = self._solvent.dielectric(temperature)
    viscosity = self._solvent.viscosity(temperature)

    alpha = (5.799e5 * abs(self.valence) /
             (kelvin(temperature) * dielectric)**(3./2.)
             )
    beta = (3.022588e-9 * abs(self.valence) / viscosity /
            (kelvin(temperature) * dielectric)**(1./2.))

    mobility = self.absolute_mobility(temperature)
    mobility -= (alpha * mobility +
                 beta * np.sign(self.valence)) * \
        (sqrt(2 * ionic_strength) / (1. + pitts * sqrt(2 *ionic_strength)))

    return mobility


def onsager_fuoss_mobility(self):
    """Return the Onsager-Fuoss corrected mobility of each charge state."""
    _, ionic_strength, temperature = \
        self._resolve_context(None, None, None)

    dielectric = self._solvent.dielectric(temperature)
    viscosity = self._solvent.viscosity(temperature)
    interaction = _interaction(self, self.context())

    alpha = (1.98074e6 * abs(self.valence) * interaction /
             (kelvin(temperature) * dielectric)**(3./2.)
             )
    beta = (3.022588e-9 * abs(self.valence) / viscosity /
            (kelvin(temperature) * dielectric)**(1./2.))

    mobility = self.absolute_mobility()
    mobility -= (alpha * mobility +
                 beta * np.sign(self.valence)) * \
        (sqrt(2 * ionic_strength) / (1. + pitts * sqrt(2 *ionic_strength)))

    return mobility


def _interaction(ion, solution):
    ions = [ion_ for ion_ in solution.ions
            if solution.concentration(ion_) > 0] + \
           [solution._hydroxide,  solution._hydronium]

    if ion not in ions:
        ions = ions + [ion]
    ion_index = ions.index(ion)
    if ion_index != 0:
        start_index = len(np.concatenate([ion_.valence
                                          for ion_ in ions[:ion_index]]))
    else:
        start_index = 0
    end_index = start_index + len(ion.valence)

    omega = np.concatenate([ion.absolute_mobility() / ion.valence
                            for ion in ions]) / faraday

    if np.any(omega == 0.):
        raise RuntimeError('Onsager-Fuoss approximation '
                           'diverges for non-mobile ions. ')

    valences = np.concatenate([ion.valence for ion in ions])
    concentrations = np.concatenate([solution.concentration(ion) *
                                     ion.ionization_fraction(solution.pH)
                                     for ion in ions])

    potential = concentrations * valences**2. / (2. * solution.ionic_strength)

    h = potential * omega / (omega + omega[:, np.newaxis])
    d = np.diag(np.sum(h, 1))
    B = 2 * (h + d) - np.identity(len(omega))

    r = np.zeros([len(omega), 6])
    r[:, 0] = (valences - (np.sum(valences * potential) /
                           np.sum(potential / omega)) / omega)

    for i in range(1, 6):
        r[:, i] = np.dot(B, r[:, i-1])

    factor = np.dot(onsager_fuoss, np.transpose(r))

    return factor[start_index:end_index]

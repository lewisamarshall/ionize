import numpy as np
from math import sqrt
from ..constants import faraday, kelvin, reference_temperature, onsager_fuoss


def interaction(self, ion):
    """Return the Onsager-Fuoss correction to the mobilities of ions.

    This function returns a list of all corrected actual mobilities.
    These mobilities are automatically assigned to the correct ions when
    a solution is initialized."""

    ions = [ion_ for ion_ in self.ions if self.concentration(ion_) > 0] + \
           [self._hydroxide,  self._hydronium]

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
    concentrations = np.concatenate([self.concentration(ion) *
                                     ion.ionization_fraction(self.pH)
                                     for ion in ions])

    potential = concentrations * valences**2. / (2. * self.ionic_strength)

    h = potential * omega / (omega + omega[:, np.newaxis])
    d = np.diag(np.sum(h, 1))
    B = 2 * (h + d) - np.identity(len(omega))

    r = np.zeros([len(omega), 6])
    r[:, 0] = (valences - (np.sum(valences*potential) /
                           np.sum(potential/omega)) / omega)

    for i in range(1, 6):
        r[:, i] = np.dot(B, r[:, i-1])

    factor = np.dot(onsager_fuoss, np.transpose(r))

    return factor[start_index:end_index]

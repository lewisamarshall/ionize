from math import sqrt
from ..constants import pitts
import numpy as np


def activity(self, valence=None, ionic_strength=None,
                         temperature=None):
    """Return activity coefficients of a charge state at ionic strength I."""
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    if valence is None:
        valence = self._valence_zero()
    else:
        try:
            valence = np.array([v for v in valence])
        except:
            valence = np.array([valence])

    # There are two coefficients that are used repeatedly.
    # Specified in Bahga.
    A = (self._solvent.debye_huckel(temperature)*sqrt(ionic_strength) /
         (1. + pitts * sqrt(ionic_strength))
         )
    # TODO: check if this is right
    B = 0.1*ionic_strength  # Matching STEEP implementation.

    # Use them to calculate the activity coefficients.
    # gamma = [10**(v**2*(B-A)) for v in valence]
    gamma = 10**((valence**2)*(B-A))

    return gamma

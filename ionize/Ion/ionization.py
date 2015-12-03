from __future__ import division
import numpy as np


def ionization_fraction(self, pH=None, ionic_strength=None, temperature=None):
    """Return the fraction of time the ion is in each valence state.

    Value is returned as a numpy array. This array will not sum to 1 due to
    the fraction of ion in the uncharged state.
    """
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    assert pH is not None, 'Calculation requires a pH.'

    # Compute the concentration of H+ from the pH.
    cH = 10**(-pH)/self._solvent.activity(1, ionic_strength, temperature)

    # Calculate the numerator of the function for ionization fraction.
    i_frac_vector = (self.acidity_product(ionic_strength, temperature) *
                     cH ** self._valence_zero())

    # Filter out the neutral fraction
    i_frac = i_frac_vector[self._valence_zero() != 0] / i_frac_vector.sum()

    return i_frac


def charge(self, pH=None, ionic_strength=None, temperature=None, moment=1):
    """Return the time-averaged charge of the ion.

    :param moment: Control which moment average is returned. Default is 1.
    """
    fraction = self.ionization_fraction(pH, ionic_strength, temperature)
    return np.sum(fraction * self.valence**moment)


def acidity_product(self, ionic_strength=None, temperature=None):
    """Return the products of the acidity.

    This vector, commonly referred to as L, is useful in computing the
    equilibrium pH in a solution, and to compute the ionization fraction of an
    ion.
    """
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    Ka = self.acidity(ionic_strength, temperature).tolist()
    index_0 = list(self._valence_zero()).index(0)
    Ka.insert(index_0, 1)

    Lp = np.cumprod(Ka)
    Lpp = np.cumprod(Ka[::-1])[::-1]
    L = np.where(self._valence_zero() >= 0,
                 Lp[self._valence_zero() == 0] / Lp,
                 Lpp / Lpp[self._valence_zero() == 0])

    return L

import numpy as np


def ionization_fraction(self, pH=None, ionic_strength=None, temperature=None):
    """Return the ionization fractions of an ion.

    Args:
        pH (float): The ambiant pH.

        I (float): The ambiant ionic strength.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].ionization_fraction()

    Otherwise, always call with a pH argument.
    """
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    # Compute the concentration of H+ from the pH.
    cH = 10**(-pH)/self.activity(1, ionic_strength, temperature)

    # Calculate the numerator of the function for ionization fraction.
    i_frac_vector = (self.acidity_product(ionic_strength, temperature) *
                     cH ** self._valence_zero())

    # Filter out the neutral fraction
    i_frac = i_frac_vector[self._valence_zero() != 0] / i_frac_vector.sum()

    return i_frac


def acidity_product(self, ionic_strength=None, temperature=None):
    """Return the L products of acidity constants.

    Args:
        I (float): The ambiant ionic strength.

    This function uses ionic strength to correct the Ka of ions. If no ionic
    strength is supplied, and the Ion is nested in a Solution, the solution
    ionic strength will be used. Otherwise, the ionic strength is assumed to be
    0.

    L is used by Solution selfects to calculate equilibrium pH.
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

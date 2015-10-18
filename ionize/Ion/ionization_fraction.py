import numpy as np


def ionization_fraction(self, pH, ionic_strength=0., temperature=25.):
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

    # Get the vector of products of acidity constants.
    L = self.L(ionic_strength, temperature)
    # Compute the concentration of H+ from the pH.
    cH, = 10**(-pH)/self.activity(1, ionic_strength, temperature)

    # Calculate the numerator of the function for ionization fraction.
    i_frac_vector = [Lp * cH ** z for (Lp, z) in zip(L, self._valence_zero())]

    # Calculate the vector of ionization fractions
    # Filter out the neutral fraction
    denom = sum(i_frac_vector)
    i_frac = [i/denom for (i, z) in zip(i_frac_vector,
                                        self._valence_zero()) if z]

    return np.array(i_frac)

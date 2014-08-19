def ionization_fraction(self, pH=None, I=None):
    """Return the ionization fractions of an ion.

    Args:
        pH (float): The ambiant pH.

        I (float): The ambiant ionic strength.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].ionization_fraction()

    Otherwise, always call with a pH argument.
    """
    if pH is None:
        assert self._pH, 'requires an input pH'
        pH = self._pH

    if I is None:
        if self._I:
            I = self._I
        else:
            I = 0.0

    # Get the vector of products of acidity constants.
    L = self.L(I)
    # Compute the concentration of H+ from the pH.
    cH = 10**(-pH)/self.activity_coefficient(I, [1])[0]

    # Calculate the numerator of the function for ionization fraction.
    i_frac_vector = [Lp * cH ** z for (Lp, z) in zip(L, self.z0)]

    # Calculate the vector of ionization fractions
    # Filter out the neutral fraction
    denom = sum(i_frac_vector)
    i_frac = [i/denom for (i, z) in zip(i_frac_vector, self.z0) if z]

    return i_frac

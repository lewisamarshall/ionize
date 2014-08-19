def molar_conductivity(self, pH=None, I=None):
    """Retun the molar conductivity of the ion based on the pH and I.

    Provides conducitivity in Siemens per meter per mole.

    Args:
        pH (float): The ambiant pH.

        I (float): The ambiant ionic strength.

    If an actual mobility from Onsager-Fouss is available, it is used,
    otherwise, the Robinson-Stokes mobility estimate is used.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].molar_conductivity()

    Otherwise, always call with a pH argument.
    """
    if pH is None:
        assert self._pH, 'requires an input pH'
        pH = self._pH

    if I is None:
        if self._I:
            I = self._I
        else:
            I = 0

    if self.actual_mobility:
        actual_mobility = self.actual_mobility
    else:
        actual_mobility = self.robinson_stokes_mobility(I)

    i_frac = self.ionization_fraction(pH, I)

    m_conductivity = (self._Lpm3 * self._F *
                      sum(z * f * m for (z, f, m)
                          in zip(self.z, i_frac, actual_mobility)))

    return m_conductivity

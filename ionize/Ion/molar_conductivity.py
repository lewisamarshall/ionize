from ..constants import faraday, lpm3


def molar_conductivity(self, pH=None, ionic_strength=None, temperature=None):
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
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    actual_mobility = self.actual_mobility(ionic_strength, temperature)

    i_frac = self.ionization_fraction(pH, ionic_strength, temperature)

    m_conductivity = (lpm3 * faraday *
                      sum(z * f * m for (z, f, m)
                          in zip(self.valence, i_frac, actual_mobility)))

    return m_conductivity

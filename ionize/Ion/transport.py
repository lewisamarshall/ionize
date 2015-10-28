from ..constants import faraday, lpm3, boltzmann, kelvin, elementary_charge
import numpy as np


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

    m_conductivity = (lpm3 * faraday *
                      sum(self.valence *
                          self.ionization_fraction(pH,
                                                   ionic_strength,
                                                   temperature) *
                          self.actual_mobility(ionic_strength,
                                               temperature)
                          )
                      )

    return m_conductivity


# TODO: Use charge and ionization fraction to simplify math
def diffusivity(self, pH=None, ionic_strength=None, temperature=None):
    """Return the diffusivity of the species at a specified pH and temperature.

    The diffusivity is returned in units of m^2/s.
    """
    pH, ionic_strength, temperature = self._resolve_context(pH,
                                                            ionic_strength,
                                                            temperature)
    actual_mobility = self.actual_mobility(ionic_strength, temperature)
    ionization_fraction = self.ionization_fraction(pH,
                                                   ionic_strength,
                                                   temperature)

    diffusivity = np.sum(actual_mobility *
                         ionization_fraction /
                         self.valence *
                         boltzmann * kelvin(temperature) /
                         elementary_charge) / \
        np.sum(ionization_fraction)
    return diffusivity

from __future__ import division
from ..constants import faraday, lpm3, boltzmann, kelvin, elementary_charge
import numpy as np


def molar_conductivity(self, pH=None, ionic_strength=None, temperature=None):
    """Retun the molar conductivity of the ion in S/m/M

    :param pH
    :param ionic_strength
    :param temperature
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


def diffusivity(self, pH=None, ionic_strength=None, temperature=None):
    """Return the diffusivity of the ion in m^2/s.

    :param pH
    :param ionic_strength
    :param temperature
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

"""Create the Aqueous class to hold the properties of water."""
from math import log10, log, pi, sqrt
from .constants import gas_constant, reference_temperature, \
                       kelvin_conversion, elementary_charge, avagadro,\
                       boltzmann, permittivity, lpm3


class Aqueous(object):

    """Access the properties of water."""

    _reference_dissociation = 1.E-14    # Water equilibrium constant [mol^2]
    _reference_pKw = -log10(_reference_dissociation)
    _enthalpy = 55.815e3                # enthalpy of dissociation water
    _heat_capacity = -224.              # heat capacity of water

    def dielectric(self, temperature):
        """Return the dielectric constant of water at a specified temperature.

        The temperature should be specified in Celcius. Correlation is based on
        the CRC handbook.
        """
        temperature = self.temperature_kelvin(temperature)
        dielectric_ = 249.21 - 0.79069*temperature + 0.72997e-3*temperature**2
        return dielectric_

    def viscosity(self, temperature):
        """Return the viscosity of water at the specified temperature.

        Correlation is based on Fox and McDonald's Intro to FLuid Mechanics.
        """
        temperature = self.temperature_kelvin(temperature)
        viscosity_ = 2.414e-5 * 10**(247.8 / (temperature - 140))
        return viscosity_

    def dissociation(self, temperature):
        """Return the dissociation constant of water."""
        reference_temperature_ = self.temperature_kelvin(reference_temperature)

        temperature = self.temperature_kelvin(temperature)

        enthalpy_contribution = (self._enthalpy / 2.303 / gas_constant) * \
            (1. / reference_temperature_ - 1. / temperature)

        cp_contribution = (self._heat_capacity/2.303/gas_constant) * \
            (reference_temperature_ / temperature - 1. +
             log10(temperature / reference_temperature_))

        pKw = self._reference_pKw - enthalpy_contribution - cp_contribution

        dissociation_ = 10.0**(-pKw)
        return dissociation_

    def debye_huckel(self, temperature):
        """Return the Debye-Huckel constant, in M^-(1/2)."""
        debye_huckel_ = elementary_charge**3. * sqrt(avagadro) / 2**(5./2.) / pi / \
            (self.dielectric(temperature) * permittivity *
             boltzmann * self.temperature_kelvin(temperature))**(3./2.)

        # Before returning answer, use log 10, convert from meter**3 to liter
        return debye_huckel_ / log(10.) * sqrt(lpm3)

    def temperature_kelvin(self, temperature):
        """Convert a temperature from Celsius to Kelvin."""
        return temperature + kelvin_conversion

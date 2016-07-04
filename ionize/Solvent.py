"""Create the Aqueous class to hold the properties of water."""
from __future__ import division
from math import log10, log, pi, sqrt, exp
from .constants import gas_constant, reference_temperature, \
                       kelvin, elementary_charge, avogadro,\
                       boltzmann, permittivity, lpm3, pitts


class Solvent(object):

    reference_dissociation = None
    enthalpy = None
    heat_capacity = None

    def __new__(cls, *args, **kwargs):
        raise TypeError('Solvents may not be instantiated.')

    @classmethod
    def reference_pKs(self):
        """Return the reference value of pKs at the reference temperature."""
        return -log10(self.reference_dissociation)

    @classmethod
    def dielectric(self, temperature):
        raise NotImplementedError

    @classmethod
    def viscosity(self, temperature):
        raise NotImplementedError

    @classmethod
    def dissociation(self, ionic_strength, temperature):
        """Return the dissociation constant of water."""
        if temperature == reference_temperature:
            dissociation_ = self.reference_dissociation
        else:
            reference_temperature_k = kelvin(reference_temperature)

            temperature_k = kelvin(temperature)

            enthalpy_contribution = (self.enthalpy / log(10.) /
                                     gas_constant) * \
                (1. / reference_temperature_k - 1. / temperature_k)

            cp_contribution = (self.heat_capacity / log(10.) /
                               gas_constant) * \
                (reference_temperature_k / temperature_k - 1. +
                 log10(temperature_k / reference_temperature_k))

            pKs = (self.reference_pKs() -
                   enthalpy_contribution -
                   cp_contribution)

            dissociation_ = 10.0**(-pKs)

        # correct for ionic strength
        dissociation_ /= self.activity(1., ionic_strength=ionic_strength,
                                       temperature=temperature)**2.
        return dissociation_

    @classmethod
    def debye(self, ionic_strength, temperature):
        """Return the Debye length of the solvent."""
        dielectric = self.dielectric(temperature)
        lamda = (dielectric * permittivity * boltzmann * kelvin(temperature) /
                 elementary_charge**2. /
                 (ionic_strength * lpm3) / avogadro) ** .5
        return lamda

    @classmethod
    def debye_huckel(self, temperature):
        """Return the Debye-Huckel constant, in M^-(1/2)."""
        dh = elementary_charge**3. * sqrt(avogadro) / 2.**(5./2.) / pi / \
            (self.dielectric(temperature) * permittivity *
             boltzmann * kelvin(temperature))**(3./2.)

        # Before returning answer, use log 10, convert from meter**3 to liter
        return dh / log(10.) * sqrt(lpm3)

    @classmethod
    def bjerrum(self, temperature):
        """Return the Bjerrum length of the solvent."""
        dielectric = self.dielectric(temperature)
        lamda = elementary_charge**2 /\
                4 / pi / dielectric / permittivity/ boltzmann / kelvin(temperature)
        return lamda

    @classmethod
    def ionic_strength(self, pH=None, temperature=None):
        """Return the ion contribution to ionic strength."""
        try:
            cH = 10**-pH
        except TypeError:
            cH = sqrt(self.dissociation(ionic_strength=0.,
                                        temperature=temperature))

        if temperature is None:
            temperature = self.reference_temperature

        cOH = self.dissociation(ionic_strength=0., temperature=temperature)/cH
        return (cH + cOH)/2.

    @classmethod
    def pKs(self, ionic_strength, temperature):
        """Return the pKs for the solvent."""
        return -log10(self.dissociation(ionic_strength, temperature))

    @classmethod
    def activity(self, valence, ionic_strength, temperature):
        """Return activity coefficients of a charge state."""

        # Specified in Bahga.
        A = (self.debye_huckel(temperature) * sqrt(ionic_strength) /
             (1. + pitts * sqrt(ionic_strength))
             )
        B = 0.1 * ionic_strength

        gamma = 10**((valence**2)*(B-A))

        return gamma


class Aqueous(Solvent):

    """Access the properties of water."""

    reference_dissociation = 1.E-14    # Water equilibrium constant [mol^2]
    enthalpy = 55815.                # enthalpy of dissociation water
    heat_capacity = -224.              # heat capacity of water

    @classmethod
    def dielectric(self, temperature):
        """Return the dielectric constant of water at a specified temperature.

        The temperature should be specified in Celcius. Correlation is based on
        the CRC handbook.
        """
        temperature = kelvin(temperature)
        dielectric_ = 249.21 - 0.79069*temperature + 0.72997e-3*temperature**2
        return dielectric_

    @classmethod
    def viscosity(self, temperature):
        """Return the viscosity of water at the specified temperature.

        Correlation is based on Fox and McDonald's Intro to Fluid Mechanics.
        """
        temperature = kelvin(temperature)
        viscosity_ = 2.414e-5 * 10**(247.8 / (temperature - 140))
        return viscosity_

    @classmethod
    def henry_CO2(self, temperature):
        """Returns the henry's law constant for CO2."""
        temperature = kelvin(temperature)
        reference = kelvin(reference_temperature)
        H = 0.034 * exp(2400. * (1./temperature - 1./reference))
        return H

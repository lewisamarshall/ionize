"""Module containing the Ion class."""
import warnings
from math import copysign
import json
import numpy as np

from ..Aqueous import Aqueous
from ..BaseIon import BaseIon
from ..constants import reference_temperature


class Ion(BaseIon):

    r"""Describe an ion dissolved in aqueous solution.

    Args:
        name (str): The chemical name of the ion.

        z (list): A list of valence states for the ion, as integers.

        pKa_ref (list): The pKa of each valence at the refernce temperature,\
        as floats.

        absolute_mobility_ref (list): The signed absolute mobility of each\
        valence at the reference temperature, as floats, in units\
        of m^2/V/s. Expect O(10^-8).

        dH (list): The enthalpy of dissociation of each valence, at the\
        reference temperature, as floats.

        dCp (list): The change in heat capacity of dissociation of each\
        valence, at the reference temperature, as floats.

        nightingale_data (function): A function describing absolute\
        mobility as a function of temperature, for special ions.

        T (float): The temperature to use to calculate the properties of the
        ions, in degrees C.

        T_ref (float): The reference temperature for the reference properties,
        in degrees C.

    Attributes:
        z (list): A list of valence states for the ion, as integers.

        pKa (list): The pKa of each valence at the refernce temperature,\
            as floats.

        absolute_mobility (list): The signed absolute mobility of each\
            valence at the reference temperature, as floats, in units\
            of m^2/V/s. Expect O(10^-8).

        T (float): The temperature to use to calculate the properties of the\
            ions, in degrees C.

    Raises:
        None

    Example:
        To to initialize an Ion, call as:

        >>> ionize.Ion('my_acid', [-1, -2], [1.2, 3.4], [-10e-8, -21e-8])
    """

    _state = ('name',
              'valence',
              'reference_pKa',
              'reference_mobility',
              'enthalpy',
              'heat_capacity',
              'nightingale_data')

    # The reference properties of the ion are stored and used to calculate
    # properties at the current temperature.
    reference_pKa = None
    reference_mobility = None
    reference_temperature = reference_temperature

    # The properties of the ions are stored in public variables.
    # These are the properties at the current temperature, or are treated
    # as temperature independant.
    pKa = None
    absolute_mobility = None
    enthalpy = None
    heat_capacity = None
    nightingale_data = None

    def __init__(self, name, valence, reference_pKa, reference_mobility,
                 reference_temperature=None, enthalpy=None, heat_capacity=None,
                 nightingale_data=None):
        """Initialize an Ion object."""

        self.name = name
        self.valence = np.int_(valence)
        self.reference_pKa = np.float_(reference_pKa)
        self.reference_mobility = np.float_(reference_mobility)

        if reference_temperature is not None:
            self.reference_temperature = float(reference_temperature)

        if enthalpy is not None:
            self.enthalpy = np.float_(enthalpy)

        if heat_capacity is not None:
            self.heat_capacity = np.float_(heat_capacity)

        if nightingale_data is not None:
            self.nightingale_data = nightingale_data
            self._nightingale_function = \
                np.poly1d(self.nightingale_data['fit'])

    def absolute_mobility(self, temperature):
        if self._nightingale_function:
            absolute_mobility = \
                [self._nightingale_function(self.T).tolist() *
                 10.35e-11 * z / self._solvent.viscosity(self.T)
                 for z in self.z]
            if (self.T > self.nightingale_data['max']) or \
                    (self.T < self.nightingale_data['min']):
                warnings.warn('Temperature outside range'
                              'for nightingale data.')
        else:
            absolute_mobility =\
                [self._solvent.viscosity(self._T_ref) /
                 self._solvent.viscosity(self.T)*m
                 for m in self._absolute_mobility_ref]
        return absolute_mobility

    def pKa(self, temperature):
        pass

    def temperature_adjust(self):
        """Temperature adjust the ion."""
        # if self.T == self._T_ref:
        #     self.pKa = self._pKa_ref
        #     self.absolute_mobility = self._absolute_mobility_ref
        # else:
        #     self.pKa = self._correct_pKa()
        #     if self._nightingale_function:
        #         self.absolute_mobility = \
        #             [self._nightingale_function(self.T).tolist() *
        #              10.35e-11 * z / self._solvent.viscosity(self.T)
        #              for z in self.z]
        #         if (self.T > self.nightingale_data['max']) or \
        #                 (self.T < self.nightingale_data['min']):
        #             warnings.warn('Temperature outside range'
        #                           'for nightingale data.')
        #     else:
        #         self.absolute_mobility =\
        #             [self._solvent.viscosity(self._T_ref) /
        #              self._solvent.viscosity(self.T)*m
        #              for m in self._absolute_mobility_ref]
        # After storing the ion properties, ensure that the properties are
        # sorted in order of charge. All other ion methods assume that the
        # states will be sorted by charge.
        self._z_sort()
        self._set_z0()

    def _z_sort(self):
        """Sort the charge states from lowest to highest."""
        # Zip the lists together and sort them by z.
        self.z, self.pKa, self.absolute_mobility =\
            zip(*sorted(zip(self.z, self.pKa, self.absolute_mobility)))
        self.z = list(self.z)
        self.pKa = list(self.pKa)
        self.absolute_mobility = list(self.absolute_mobility)

        full = set(range(min(self.z), max(self.z)+1, 1)) - set([0])
        assert set(self.z) ^ full == set(), "Charge states missing."

        return None

    def acidity(self):
        """Return the acidity constant, Ka, based on the pKa."""
        return 10**(-self.pKa())

    def _valence_zero(self):
        """Create a list of charge states with 0 inserted."""
        return np.sort(np.append(self.valence, [0]))

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

    def actual_mobility(self, ionic_strength=None, temperature=None):
        if ionic_strength is None and \
                temperature is None and \
                self.context() is not None:
            try:
                return self.context().actual_mobility(self)
            except:
                warnings.warn('Context failed to return an actual mobility.')
        return self.robinson_stokes_mobility(ionic_strength, temperature)


    from .ionization_fraction import ionization_fraction
    from .activity import activity
    from .effective_mobility import effective_mobility
    # from .Ka_eff import Ka_eff
    from .pKa import pKa, Ka
    from .L import L
    from .molar_conductivity import molar_conductivity
    from .robinson_stokes_mobility import robinson_stokes_mobility
    # from .correct_pKa import pKa, _vant_hoff, _clark_glew
    from .pKa import pKa, Ka, mid_Ka, mid_pKa

if __name__ == '__main__':
    pass

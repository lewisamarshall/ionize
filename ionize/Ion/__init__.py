"""Module containing the Ion class."""
import warnings
from math import copysign
import json
import numpy as np

from ..Aqueous import Aqueous
from ..BaseIon import BaseIon
from ..constants import reference_temperature, boltzmann, kelvin, \
                        elementary_charge
from ..fixed_state import fixed_state


@fixed_state
class Ion(BaseIon):

    """Describe an ion dissolved in aqueous solution."""

    _state = ('name',
              'valence',
              'reference_pKa',
              'reference_mobility',
              'reference_temperature',
              'enthalpy',
              'heat_capacity',
              'nightingale_data')

    # The reference properties of the ion are stored.
    _reference_pKa = None
    _reference_mobility = None
    _reference_temperature = reference_temperature
    _pKa = None
    _absolute_mobility = None
    _enthalpy = None
    _heat_capacity = None
    _nightingale_data = None

    def __init__(self, name, valence, reference_pKa, reference_mobility,
                 reference_temperature=None, enthalpy=None, heat_capacity=None,
                 nightingale_data=None):
        """Initialize an Ion object."""

        self._name = name
        self._valence = np.int_(valence)

        # TODO: Re-impliment sorting?
        assert np.all(np.diff(self.valence) > 0), 'Valences must be sorted.'

        self._reference_pKa = np.float_(reference_pKa)
        self._reference_mobility = np.float_(reference_mobility)

        if reference_temperature is not None:
            self._reference_temperature = float(reference_temperature)

        if enthalpy is not None:
            self._enthalpy = np.float_(enthalpy)
            assert len(enthalpy) == len(self.reference_pKa)

        if heat_capacity is not None:
            self._heat_capacity = np.float_(heat_capacity)

        if nightingale_data is not None:
            self._nightingale_data = nightingale_data
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
    from .L import L
    from .molar_conductivity import molar_conductivity
    from .robinson_stokes_mobility import robinson_stokes_mobility
    from .pKa import pKa, acidity, mid_Ka, mid_pKa

    # from .mobility import absolute_mobility, actual_mobility, mobility, robinson_stokes_mobility
    # from .acidity import pKa, acidity
    # from .ionization import L,ionization_fraction

if __name__ == '__main__':
    pass

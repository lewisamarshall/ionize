"""Module containing the Ion class."""
import warnings
from math import copysign
import json
import numpy as np

from ..Aqueous import Aqueous
from .BaseIon import BaseIon
from ..constants import reference_temperature, boltzmann, kelvin, \
                        elementary_charge
from .fixed_state import fixed_state


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
    _valence = None
    _reference_pKa = None
    _reference_mobility = None
    _reference_temperature = reference_temperature
    _pKa = None
    _absolute_mobility = None
    _enthalpy = None
    _heat_capacity = None
    _nightingale_data = None

    # Polynomial function is a derived parameter.
    _nightingale_function = None

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

        for prop in ('reference_pKa',
                     'reference_mobility',
                     'enthalpy',
                     'heat_capacity'):
            if getattr(self, prop) is not None:
                assert getattr(self, prop).shape == self.valence.shape, \
                    '{} must have the same shape as valence.'.format(prop)

        if nightingale_data is not None:
            self._nightingale_data = nightingale_data
            self._nightingale_function = \
                np.poly1d(self.nightingale_data['fit'])

    def _valence_zero(self):
        """Create a list of charge states with 0 inserted."""
        return np.sort(np.append(self.valence, [0]))

    from .acidity import pKa, acidity, activity, _clark_glew_pKa, \
        _clark_glew_acidity, _vant_hoff_acidity, _vant_hoff_pKa

    from .ionization import acidity_product, ionization_fraction

    from .mobility import absolute_mobility, actual_mobility, \
        mobility, robinson_stokes_mobility

    from .transport import molar_conductivity, diffusivity

if __name__ == '__main__':
    pass

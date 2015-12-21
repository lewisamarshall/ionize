"""Module containing the Ion class."""
from __future__ import division
import warnings
from math import copysign
import json
import numpy as np

from .BaseIon import BaseIon
from ..constants import reference_temperature, boltzmann, kelvin, \
                        elementary_charge
from .fixed_state import fixed_state


@fixed_state
class Ion(BaseIon):
    """An Ion describes a charged species in solution.

    Ion is the most commonly used subclass of BaseIon. It is used to represent
    small ions that have a set of known valence states with distinct pKas. Ions
    are immutable.

    Example::

        acid = ionize.Ion('acid', [-1], [3], [30e-9], molecular_weight=19)
        acid.mobility(pH=7)
        acid.diffusivity(pH=7)

    :param name: The name of the ion.

    :param valence: An iterable of the integer valence states.

    :param reference_pKa: An iterable of the pKas associated with the valence
        states at the reference temperature.

    :param reference_mobility: An iterable of the fully ionized mobility of
        each valence state at infinite dilution in m^2/V/s at the reference
        temperature.

    :param reference_temperature: The temperature at which other parameters are
        measured. Defaults to 25 degrees C.

    :param enthalpy: The enthalpy change on dissociation for each valence
        state. This is optional, but allows more accurate models for calculating
        temperature-dependance of properties.

    :param heat_capacity: The change in heat capacity on dissociation for
        each valence state. This optional parameter further improves the accuracy
        of temperature-dependant property calculation.

    :param nightingale_data: Mobiliity correction data for small ions where
        hydration shell dynamics are important.

    :param molecular_weight: The molecular weight of uncharged species, in
        Daltons.

    :param alias: An iterable of alias strings that can be used to refer to
        the ion.
    """

    _state = {'name': 'The ion name.',
              'valence': 'The valence of each state.',
              'reference_pKa': 'The pKas of each ionization state at the reference state.',
              'reference_mobility': 'The mobility of each ionization at the reference state.',
              'reference_temperature': 'The temperature at whcih properties were measured.',
              'enthalpy': 'The change in enthalpy on ionization.',
              'heat_capacity': 'The change in heat capacity on ionization.',
              'nightingale_data': 'Temperature dependance data for small ions.',
              'molecular_weight': 'The ion molecular weight',
              'alias': 'Alternative chemical names.'}

    # The reference properties of the ion are stored in private variables.
    _valence = None
    _reference_pKa = None
    _reference_mobility = None
    _reference_temperature = reference_temperature
    _pKa = None
    _absolute_mobility = None
    _enthalpy = None
    _heat_capacity = None
    _nightingale_data = None
    _molecular_weight = None
    _alias = None

    # Polynomial function is a derived parameter.
    _nightingale_function = None

    def __init__(self, name, valence, reference_pKa, reference_mobility,
                 reference_temperature=None, enthalpy=None, heat_capacity=None,
                 nightingale_data=None, molecular_weight=None, alias=None):
        """Initialize an Ion object."""

        self._name = str(name)
        if alias is not None:
            self._alias = tuple(alias)
        self._valence = np.int_(valence)

        if len(self.valence) > 1:
            assert np.all(np.diff(self.valence) > 0), \
                'Valences must be sorted.'

        self._reference_pKa = np.float_(reference_pKa)
        self._reference_mobility = np.float_(reference_mobility)

        if reference_temperature is not None:
            self._reference_temperature = float(reference_temperature)

        if enthalpy is not None:
            self._enthalpy = np.float_(enthalpy)
            assert len(enthalpy) == len(self.reference_pKa)

        if heat_capacity is not None:
            self._heat_capacity = np.float_(heat_capacity)

        self._molecular_weight = molecular_weight

        for prop in ('reference_pKa',
                     'reference_mobility',
                     'enthalpy',
                     'heat_capacity'):
            if getattr(self, prop) is not None:
                assert getattr(self, prop).shape == self.valence.shape, \
                    '{} must have the same shape as valence.'.format(prop)

        assert np.all((self.reference_mobility / self.valence) > 0.), \
            'Mobilities must be signed. {}, {}'.format(self.reference_mobility,
                                                       self.valence)

        if nightingale_data is not None:
            self._nightingale_data = nightingale_data
            self._nightingale_function = \
                np.poly1d(self.nightingale_data['fit'])

    def _valence_zero(self):
        """Create a list of charge states with 0 inserted."""
        return np.sort(np.append(self.valence, [0]))

    from .acidity import pKa, acidity, _clark_glew_pKa, \
        _clark_glew_acidity, _vant_hoff_acidity, _vant_hoff_pKa

    from .ionization import acidity_product, ionization_fraction, charge

    from .mobility import absolute_mobility, actual_mobility, \
        mobility, robinson_stokes_mobility, onsager_fuoss_mobility

    from .transport import molar_conductivity, diffusivity

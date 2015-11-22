"""Module containing the BaseIon class."""
import json
import numpy as np
import contextlib
import operator
from math import sqrt
import warnings

from .fixed_state import fixed_state
from ..Solvent import Aqueous
from ..constants import reference_temperature
from ..serialize import _serialize


@fixed_state
class BaseIon(object):

    """BaseIon class describing basic ion properties.

    All ions on ionize should be subclassed from BaseIon.
    """

    _solvent = Aqueous

    _state = ('name', 'reference_temperature')
    _name = 'BaseIon'
    _reference_temperature = reference_temperature

    _context = None

    def __repr__(self):
        inner = []
        for prop in self._state:
            prop = str(prop)  # convert unicode to string
            value = getattr(self, prop)
            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif isinstance(value, str):
                value = str(value)
            inner.append('{}={}'.format(str(prop), repr(value)))
        return '{}({})'.format(type(self).__name__, ', '.join(inner))

    def __str__(self):
        return "{}('{}')".format(type(self).__name__, self.name)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        try:
            assert self._state == other._state
            for prop in self._state:
                if isinstance(getattr(self, prop), np.ndarray):
                    assert(np.all(getattr(self, prop) == getattr(other, prop)))
                else:
                    assert getattr(self, prop) == getattr(other, prop)
            return True
        except Exception as e:
            return False

    def serialize(self, nested=False, compact=False):
        """Return a serialized representation of the ion.

        If "nested" is true, returns a dictionary. Otherwise, returns
        a json-formatted string. This string can be returned to an Ion object
        using deserialize.
        """
        serial = {key: getattr(self, key) for key in self._state}
        serial['__ion__'] = type(self).__name__

        return _serialize(serial, nested, compact)

    # TODO: figure out why this dumps text
    def save(self, filename):
        """Save a serialized Ion to a file."""
        with open(filename, 'w') as file:
            file.write(self.serialize())

    def mobility(self):
        """The effective mobility of the ion.

        Mobility must be overridden by subclasses."""
        raise NotImplementedError

    def diffusivity(self):
        """The diffusivity of the ion.

        Diffusivity must be overridden by subclasses."""
        raise NotImplementedError

    def molar_conductivity(self):
        """The molar conductivity of the ion, in S/m/M

        molar_conductivity must be overridden by subclasses."""
        raise NotImplementedError

    def charge(self, moment=1):
        """The charge of the ion.

        molar_conductivity must be overridden by subclasses."""
        raise NotImplementedError

    def context(self, context=False):
        """Control the context of the ion.

        The context is used as a convenience to set pH, temperature,
        ionic_strength, and other information used in calculations about
        the ion. If no

        """
        if context is False:
            return self._context

        # Update to new context, save old_context
        old_context, self._context = self._context, context

        # Return a context manager to revert to old_context
        @contextlib.contextmanager
        def manager():
            yield
            self._context = old_context
        return manager()

    def _resolve_context(self, pH, ionic_strength, temperature):
        if pH is None:
            try:
                pH = self.context().pH
            except AttributeError:
                pH = None

        try:
            temperature = temperature or \
                          self.context().temperature() or \
                          self.reference_temperature
        except AttributeError:
            temperature = self.reference_temperature

        if pH is not None:
            lower_limit = self._solvent.ionic_strength(pH, temperature)
        else:
            lower_limit = sqrt(self._solvent.dissociation(0., temperature))

        if ionic_strength is None:
            try:
                ionic_strength = self.context().ionic_strength
            except AttributeError:
                ionic_strength = lower_limit

        return pH, ionic_strength, temperature

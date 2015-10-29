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
        """Return a representation of the ion."""
        inner = []
        for prop in self._state:
            prop = str(prop)  # convert unicode to string
            value = getattr(self, prop)
            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif isinstance(value, basestring):
                value = str(value)
            inner.append('{}={}'.format(str(prop), repr(value)))
        return '{}({})'.format(type(self).__name__, ', '.join(inner))

    def __str__(self):
        return "{}('{}')".format(type(self).__name__, self.name)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        """Equality compares type and state."""
        try:
            return self.serialize() == other.serialize()
        except:
            return False

    def serialize(self, nested=False, compact=False):
        """Return a serialized representation of the ion.

        If "nested" is true, returns a dictionary. Otherwise, returns
        a json-formatted string.
        """
        serial = {key: getattr(self, key) for key in self._state}
        serial['__ion__'] = type(self).__name__

        return _serialize(serial, nested, compact)

    # TODO: figure out why this dumps text
    def save(self, filename):
        """Save a serialized version of the ion to a file."""
        with open(filename, 'w') as file:
            json.dump(self.serialize(), file)

    def mobility(self):
        """Mobility must be overridden by subclasses."""
        raise NotImplementedError

    def diffusivity(self):
        """Diffusivity must be overridden by subclasses."""
        raise NotImplementedError

    def molar_conductivity(self):
        raise NotImplementedError

    def charge(self, moment=1):
        raise NotImplementedError

    def context(self, context=False):
        """Set or get the context."""
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

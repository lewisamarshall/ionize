"""Module containing the BaseIon class."""
import json
import numpy as np
import contextlib
import operator

from .fixed_state import fixed_state
from ..Solvent import Aqueous
from ..constants import reference_temperature


def _encode(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    # Let the base class default method raise the TypeError
    return json.JSONEncoder().default(obj)


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

    # TODO: should list all parameters in state
    def __repr__(self):
        """Return a representation of the ion."""
        return "{}('{}')".format(type(self).__name__,
                                             self.name)

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

        if compact:
            sort_keys, indent, separators = True, None, (',', ':')
        else:
            sort_keys, indent, separators = True, 4, (', ', ': ')

        if nested:
            return serial
        else:
            return json.dumps(serial, default=_encode, sort_keys=sort_keys,
                              indent=indent, separators=separators)

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
        return manager

    def _resolve_context(self, pH, ionic_strength, temperature):
        if pH is None:
            try:
                pH = self.context().pH
            except AttributeError:
                pH = None
        try:
            ionic_strength = ionic_strength or self.context().ionic_strength or 0.
        except AttributeError:
            ionic_strength = 0.

        try:
            temperature = temperature or \
                          self.context().temperature() or \
                          self.reference_temperature
        except AttributeError:
            temperature = self.reference_temperature
        # print "Context(pH: {}, I: {}, T: {})".format(pH,
        #                                              ionic_strength,
        #                                              temperature)
        return pH, ionic_strength, temperature

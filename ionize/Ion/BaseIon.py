"""Module containing the BaseIon class."""
import numpy as np
import contextlib
from math import sqrt
import json

from .fixed_state import fixed_state
from ..Solvent import Aqueous
from ..constants import reference_temperature, elementary_charge, lpm3, \
    boltzmann, faraday, kelvin
from ..serialize import _serialize


@fixed_state
class BaseIon(object):
    """BaseIon is the basic implementation of a ion.

    It represents a charged species in solution. All ion implementations in
    ionize are subclassed from BaseIon. The most common version is the
    :class:`Ion`.
    """

    _solvent = Aqueous

    _state = {'name': 'The name of the ion.',
              'reference_temperature': 'The temperature at which ion '
                                       'properties were measured.'}
    _name = 'BaseIon'
    _alias = None
    _reference_temperature = reference_temperature

    _context = None

    def __repr__(self):
        """Return an unambiguous string representation."""
        inner = []
        for prop in sorted(self._state):
            prop = str(prop)  # convert unicode to string
            value = getattr(self, prop)
            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif isinstance(value, str):
                value = str(value)
            inner.append('{}={}'.format(str(prop), repr(value)))
        return '{}({})'.format(type(self).__name__, ', '.join(inner))

    def __str__(self):
        """Return a readable string representation."""
        return "{}('{}')".format(type(self).__name__, self.name)

    def __hash__(self):
        """Return the hash value for the object."""
        return hash(self.serialize())

    def __eq__(self, other):
        """Test equality between two ions."""
        try:
            assert self._state == other._state
            for prop in self._state:
                if isinstance(getattr(self, prop), np.ndarray):
                    assert getattr(other, prop) is not None
                    assert(np.all(getattr(self, prop) == getattr(other, prop)))
                else:
                    if isinstance(getattr(other, prop), np.ndarray):
                        assert getattr(self, prop) is not None
                    assert getattr(self, prop) == getattr(other, prop)
            return True
        except:
            return False

    def serialize(self, nested=False, compact=False):
        """Returns a JSON serialized string representation of the ion.
        This representation can be returned to an Ion object using deserialize.

        :param nested: If True, skips JSON serialization and returns a python
        dictionary.

        :param compact: If True, returns a JSON representation with minimal
        characters. Otherwise, whitespace is introduced for readability.
        """
        serial = {key: getattr(self, key) for key in self._state}
        serial['__ion__'] = type(self).__name__

        return _serialize(serial, nested, compact)

    def save(self, filename):
        """Save a serialized Ion to a file."""
        with open(filename, 'w') as file:
            file.write(self.serialize())

    def mobility(self):
        """The effective mobility of the ion, in meter^2/Volt/second.

        Mobility must be overridden by subclasses."""
        raise NotImplementedError

    def diffusivity(self, pH=None, ionic_strength=None, temperature=None):
        """The diffusivity of the ion, in meter^2/second."""
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        return (self.mobility(pH, ionic_strength, temperature) /
                self.charge(pH, ionic_strength, temperature) /
                elementary_charge * boltzmann * kelvin(temperature)
                )

    def molar_conductivity(self, pH=None, ionic_strength=None,
                           temperature=None):
        """The molar conductivity of the ion, in Seimens/meter/Molar."""
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        return (lpm3 * faraday *
                self.mobility(pH, ionic_strength, temperature) *
                self.charge(pH, ionic_strength, temperature)
                )

    def charge(self, pH=None, ionic_strength=None, temperature=None,
	       moment=1):
        """The average charge of the ion divided by the charge of an electron.

        Charge must be overridden by subclasses."""
        raise NotImplementedError

    def context(self, context=False):
        """Control the context in which ion properties are calculated.

        The context is used as a convenience to specify what values of pH,
        ionic strength, and temperature are used to compute ion properties.

        Example::
            ``water = ionize.Solution()

            my_ion.mobility()  # Raises without a pH for calculation

            with my_ion.context(water):
                my_ion.mobility()  # Uses water pH

            my_ion.context(water)
            my_ion.mobility()      # Uses water pH
            my_ion.context()       # Returns water
            my_ion.context(None)``

        There are multiple sources from which pH, ionic strength, and
        temperature can be drawn. They are, in order of priority::

            Method parameters. Manually setting a value when calling an ion
            method always has the highest priority.

            Context.

            Default values. The default temperature is the package reference
            temperature, 25 degrees C. This is assumed throughout the package
            if temperature is not specified. The default ionic strength is the
            solvent ionic strength. This ionic strength is calculated from the
            pH if available. Otherwise, it is calculated as the minimum ionic
            strength associated with the solvent dissociation constant. There
            is no default value for the pH.


        :param context: The new context to use for ion calculations. If context
        is False, the current context is returned. Context can be None (to
        erase the current context) or a Solution object.
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

        if temperature is None:
            try:
                temperature = self.context().temperature() or \
                              self.reference_temperature
            except AttributeError:
                temperature = self.reference_temperature

        if ionic_strength is None:
            try:
                ionic_strength = self.context().ionic_strength
            except AttributeError:
                if pH is not None:
                    ionic_strength = self._solvent.ionic_strength(pH,
                                                                  temperature)
                else:
                    ionic_strength = \
                        sqrt(self._solvent.dissociation(0., temperature))

        return pH, ionic_strength, temperature

    def separability(self, other, pH=None,
                     ionic_strength=None, temperature=None):
        """Return the separability between this ion and the other ion.

        Separabillity is a normalized difference between this ion and the other
        ion, defined as:
            abs((self.mobility - other.mobility)/ self.mobility)

        :param other: The other ion for comparison. Must have valid context
        and mobility methods.

        The context from this ion is used to override the context of the
        other ion.
        """
        with other.context(self.context()):
            mobilities = [ion.mobility(pH, ionic_strength, temperature)
                          for ion in (self, other)]

            # Self.mobility in denominator
            return abs((mobilities[0] - mobilities[1]) / mobilities[0])

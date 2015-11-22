from __future__ import division
import json
import copy
from collections import OrderedDict
import numbers
import contextlib
import operator
import numpy as np

from ..Ion import Ion
from ..Solvent import Aqueous
from ..Database import Database
from ..serialize import _serialize
from ..constants import permittivity, avogadro, boltzmann, \
                        elementary_charge, lpm3, reference_temperature

database = Database()


class Solution(object):

    """Describe an aqueous solution.

    Args:
        ions (list): A list of valence states for the ion, as integers.

        concentrations (list): The pKa of each valence at the refernce
        temperature, as floats.

        T (float): The temperature to use to calculate the properties of the
        ions, in degrees C.

    Attributes:
        ions (list): A list of valence states for the ion, as integers.

        concentrations (list): The pKa of each valence at the reference
        temperature, as floats.

        T (float): The temperature to use to calculate the properties of the
        ions, in degrees C.

    Raises:
        None

    Example:
        To to initialize an Soltuion, call as:

        >>> ionize.Solution([ion1, ion2], [c1, c2], T=30)

        When a new Solution is initialized, it will immediately calculate the
        equilibrium state, including the pH and the ionic strength (I) of the
        solution. These values willl be stored as permenant attributes of the
        object. Other solution properties can be calculated by invoking the
        appropriate method.
    """

    _solvent = Aqueous
    _hydronium = None
    _hydroxide = None

    _pH = 7.                # Normal pH units.
    _ionic_strength = 0.    # Expected in molar.
    _interaction_matrix = None
    _temperature = reference_temperature  # Temperature in C

    # Ions and concentrations fully represent state.
    # class uses _contents internally to enforce consistency
    # and help with lookups
    _state = ('ions', 'concentrations')
    _contents = OrderedDict()

    @property
    def _name_lookup(self):
        return {ion.name: ion for ion in self.ions}

    @property
    def ions(self):
        """Return a tuple of the ions in the solution."""
        return tuple(self._contents.keys())

    @property
    def concentrations(self):
        """Return a numpy array of the ion concentrations in the solution."""
        return np.array(list(self._contents.values()))

    pH = property(operator.attrgetter("_pH"))
    ionic_strength = property(operator.attrgetter("_ionic_strength"))

    def __init__(self, ions=[], concentrations=[], temperature=None):
        """Initialize a solution object."""

        try:
            len(ions)
        except:
            ions = (ions,)

        try:
            len(concentrations)
        except:
            concentrations = (concentrations,)

        assert len(ions) == len(concentrations), \
            'There must be an ion for each concentration.'

        self._contents = OrderedDict()
        for ion, concentration in zip(ions, concentrations):
            if isinstance(ion, str):
                ion = database.load(ion)
            else:
                ion = copy.copy(ion)
            ion.context(self)
            assert concentration >= 0, 'Concentrations must be positive.'
            self._contents[ion] = concentration

        # TODO: move these into database.
        self._hydronium = Ion('H+', [1], [100], [362E-9])
        self._hydroxide = Ion('OH-', [-1], [-100], [-205E-9])
        self._hydronium = copy.copy(self._hydronium)
        self._hydroxide = copy.copy(self._hydroxide)
        self._hydronium.context(self)
        self._hydroxide.context(self)

        if temperature is not None:
            self.temperature(temperature)
        self._equilibrate()

    def temperature(self, temperature=None):
        """Set or get the temperature of the solution.

        If no argument is supplied, returns the current temperature of the
        solution.

        If a numerical temperature is supplied, returns a context manager that reverts
        to the original temperature.
        """
        if temperature is None:
            return self._temperature
        else:
            old_temperature = self._temperature
            if temperature != old_temperature:
                self._temperature = float(temperature)
                self._equilibrate()
        return self._manage_temperature(old_temperature)

    @contextlib.contextmanager
    def _manage_temperature(self, old_temperature):
        yield
        self.temperature(old_temperature)

    def _cH(self):
        """Return the concentration of protons in solution."""
        cH = 10**(-self.pH)/self._solvent.activity(1, self.ionic_strength,
                                                   self.temperature())
        return cH

    def _cOH(self):
        """Return the concentration of hydroxyls in solution."""
        cOH = (self._solvent.dissociation(self.ionic_strength,
                                          self.temperature()) /
               self._cH() /
               self._solvent.activity(1, self.ionic_strength,
                                      self.temperature()) ** 2.)
        return cOH

    def concentration(self, ion):
        """Return the concentration of the input ion.

        The input may be an Ion or an ion name as a string.
        """
        if ion in ('H+', self._hydronium):
            return self._cH()
        elif ion in ('OH-', self._hydroxide):
            return self._cOH()
        else:
            ion = self._name_lookup.get(ion, ion)
            return self._contents.get(ion, 0)

    def __add__(self, other):
        new_i = list(self.ions)
        new_c = self.concentrations.tolist()
        if isinstance(other, Solution):
            for ion, c in zip(other.ions, other.concentrations):
                if ion in self.ions:
                    new_c[self.ions.index(ion)] += c
                else:
                    new_i.append(ion)
                    new_c.append(c)
            return Solution(new_i, new_c)
        elif isinstance(other, (list, tuple)) and\
                len(other) == 2 and\
                isinstance(other[0], Ion) and\
                isinstance(other[1], (int, float)):
            ion, c = other
            if ion in self.ions:
                new_c[self.ions.index(ion)] += c
            else:
                new_i.append(ion)
                new_c.append(c)
            return Solution(new_i, new_c)
        else:
            raise TypeError

    __radd__ = __add__

    def __mul__(self, other):
        if other >= 0:
            return Solution(self.ions,
                            [c * other for c in self.concentrations])
        else:
            raise TypeError

    __rmul__ = __mul__

    def __str__(self):
        """Return a string representing the Solution."""
        return "Solution(pH={:.3g}, I={:.3g} M)".format(self.pH,
                                                        self.ionic_strength)

    def __repr__(self):
        """Return a representation of the Solution."""
        template = "Solution(ions={}, concentrations={}, temperature={})"
        return template.format([ion.name for ion in self.ions],
                               self.concentrations.tolist(),
                               self.temperature())

    def __len__(self):
        return len(self.ions)

    def __eq__(self, other):
        try:
            return self.serialize() == other.serialize()
        except:
            return False

    def __contains__(self, other):
        return other in self.ions or other in self._name_lookup.keys()

    def __iter__(self):
        return (ion for ion in self.ions)

    def __getitem__(self, item):
        try:
            if hasattr(item, 'name'):
                return self._name_lookup[item.name]
            else:
                return self._name_lookup[item]
        except:
            raise KeyError(item)

    def serialize(self, nested=False, compact=False):
        """Return a JSON-formatted serialization of the object."""
        serial = {'__solution__': True}
        serial['concentrations'] = self.concentrations
        serial['ions'] = self.ions
        return _serialize(serial, nested, compact)

    # TODO: Figure out why this dumps text.
    def save(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.serialize(), file)

    from .equilibrium import _equilibrate
    from .conductivity import conductivity, hydroxide_conductivity, \
        hydronium_conductivity
    from .titrate import titrate, buffering_capacity, \
        equilibrate_CO2, displace
    from .debye import debye
    from .transference import transference, zone_transfer
    from .conservation import kohlrausch, alberty, jovin, gas

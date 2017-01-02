from __future__ import division
import copy
from collections import OrderedDict
import contextlib
import numpy as np

from ..Ion import BaseIon
from ..Solvent import Aqueous
from ..Database import Database
from ..serialize import _serialize
from ..constants import reference_temperature

database = Database()


class Solution(object):

    """Represent an aqueous Solution.

    Solutions automatically compute their equilibrium state, and store it
    as pH and ionic_strength.

    :param ions: An iterable of ion objects, or strings. Each string is
    converted to an ion object using the ionize database. If ion objects are
    used, new copies are stored inside the Solution. Each ion inside the
    Solution has the Solution as their context by default. However, the
    original version of the ion retains its original context.

    :param concentrations: An iterable of concentrations, in moles per liter.

    Example:
        ``sol = ionize.Solution(['chloride', 'tris'], [0.02, 0.05])
        sol.pH, sol.ionic_strength
        'tris' in sol  # True
        sol.concentration('chloride') # 0.02
        with sol.temperature(35):
            sol.conductivity()
        sol.serialize()
        sol['tris']  # Ion('tris')``
    """

    _solvent = Aqueous
    # TODO: Create a _solvent_components member.
    _hydronium = None
    _hydroxide = None

    _pH = 7.                # Normal pH units.
    _ionic_strength = 0.    # Expected in molar.
    _temperature = reference_temperature  # Temperature in C

    # Ions and concentrations fully represent state.
    # class uses _contents internally to enforce consistency
    # and help with lookups
    _state = ('ions', 'concentrations')
    _contents = OrderedDict()

    @property
    def _name_lookup(self):
        name_lookup = dict()
        for ion in self.ions:
            name_lookup[ion.name] = ion
            try:
                for alias in ion.alias:
                    name_lookup[alias] = ion
            except (AttributeError, TypeError):
                pass
        name_lookup['H+'] = self._hydronium
        name_lookup['OH-'] = self._hydroxide
        name_lookup['hydronium'] = self._hydronium
        name_lookup['hydroxide'] = self._hydroxide
        return name_lookup

    @property
    def ions(self):
        """Return a tuple of the ions in the solution."""
        return tuple(self._contents.keys())

    @property
    def concentrations(self):
        """Return a numpy array of the ion concentrations in the solution."""
        return np.array(list(self._contents.values()))

    @property
    def pH(self):
        """The pH of the solution."""
        return self._pH

    @property
    def ionic_strength(self):
        """The ionic strength of the solution."""
        return self._ionic_strength

    def __init__(self, ions=[], concentrations=[]):
        """Initialize a solution object."""

        if isinstance(ions, (str, BaseIon)):
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
            if concentration == 0:
                continue
            elif concentration <0:
                raise ValueError('Concentrations must be positive.')
            else:
                self._contents[ion] = concentration

        self._hydronium = database['hydronium']
        self._hydroxide = database['hydroxide']
        self._hydronium.context(self)
        self._hydroxide.context(self)

        for solvent_component in self._hydroxide, self._hydronium:
            assert solvent_component not in self.ions, \
                "Solvent components cannot be manually added to Solution."

        self._equilibrate()

    def temperature(self, temperature=None):
        """Set or get the temperature of the solution.

        :param temperature: The new temperature. If temperature is None, return
        the current temperature. Otherwise, sets the temperature and returns
        a context manager that reverts to the original temperature.
        """
        if temperature is None:
            return self._temperature
        else:
            old_temperature = self._temperature
            if temperature != old_temperature:
                self._temperature = float(temperature)
                self._equilibrate()

            @contextlib.contextmanager
            def manage_temperature():
                yield
                self.temperature(old_temperature)

            return manage_temperature()

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
        """Return the concentration of the ion.

        The input may be an ion object or an ion name as a string.
        """
        if ion in ('H+', self._hydronium):
            return self._cH()
        elif ion in ('OH-', self._hydroxide):
            return self._cOH()
        else:
            ion = self._name_lookup.get(ion, ion)
            return self._contents.get(ion, 0)

    def __add__(self, other):
        if isinstance(other, Solution):
            ions = list(set(self.ions + other.ions))
            return Solution(ions, [self.concentration(ion) +
                                   other.concentration(ion)
                                   for ion in ions]
                            )
        else:
            try:
                ion, concentration = other
                new_contents = dict(self._contents)
                new_contents[ion] = self.concentration(ion) + concentration
                return Solution(new_contents.keys(), new_contents.values())
            except:
                raise TypeError('Solutions add to other Solutions or to an'
                                 '(Ion, concentration) iterable pair.')

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, Solution):
            ions = list(set(self.ions + other.ions))
            return Solution(ions, [self.concentration(ion) -
                                   other.concentration(ion)
                                   for ion in ions]
                            )
        else:
            try:
                ion, concentration = other
                new_contents = dict(self._contents)
                new_contents[ion] = self.concentration(ion) - concentration
                if new_contents[ion]==0:
                    del new_contents[ion]
                return Solution(new_contents.keys(), new_contents.values())
            except:
                raise TypeError('Solutions add to other Solutions or to an'
                                 '(Ion, concentration) iterable pair.')

    def __mul__(self, other):
        if other >= 0:
            return Solution(self.ions,
                            [c * other for c in self.concentrations])
        else:
            raise TypeError

    __rmul__ = __mul__

    def __truediv__(self, other):
        if other > 0:
            return Solution(self.ions,
                            [c / other for c in self.concentrations])
        else:
            raise TypeError

    def __str__(self):
        """Return a string representing the Solution."""
        return "Solution(pH={:.3g}, I={:.3g} M)".format(self.pH,
                                                        self.ionic_strength)

    def __repr__(self):
        """Return a representation of the Solution."""
        template = "Solution(ions={}, concentrations={})"
        return template.format([ion.name for ion in self.ions],
                               self.concentrations.tolist())

    def __len__(self):
        return len(self.ions)

    def __eq__(self, other):
        try:
            return self.serialize() == other.serialize()
        except:
            return False

    def __hash__(self):
        return hash(self.serialize())

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

    def safe(self):
        """Return True if the solution has a safe pH.

        We define moderate pH as the regime in which the
        conductivity contributed by H+ and OH- is less
        than 10% the sum of charge contributed by other ions.
        """
        diss = sum([self.concentration(ion) * self[ion].molar_conductivity()
                    for ion in ('H+', 'OH-')])
        other = sum([self.concentration(ion) * ion.molar_conductivity()
                     for ion in self.ions])
        return 10. * diss < other

    def moderate(self):
        """Return True if the solution has a moderate pH.

        We define moderate pH as the regime in which the
        sum of charge contributed by H+ and OH- is less
        than 10% the sum of charge contributed by other ions.
        """
        diss = sum([abs(self.concentration(ion) * self[ion].charge())
                    for ion in ('H+', 'OH-')])
        other = sum([abs(self.concentration(ion) * ion.charge())
                     for ion in self.ions])
        return 10. * diss < other

    def serialize(self, nested=False, compact=False):
        """Return a JSON-formatted serialization of the object."""
        serial = {'__solution__': True}
        serial['concentrations'] = self.concentrations
        serial['ions'] = self.ions
        return _serialize(serial, nested, compact)

    def save(self, filename):
        """Save the solution to file."""
        with open(filename, 'w') as file:
            file.write(self.serialize())

    from .equilibrium import _equilibrate
    from .conductivity import conductivity, hydroxide_conductivity, \
        hydronium_conductivity
    from .titrate import titrate, buffering_capacity, \
        equilibrate_CO2, displace
    from .debye import debye
    from .transference import transference, zone_transfer
    from .conservation import kohlrausch, alberty, jovin, gas

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
from ..constants import permittivity, avagadro, boltzmann, \
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
    # TODO: H+ and OH- context may be incorrect.
    _hydronium = Ion('H+', [1], [100], [362E-9])
    _hydroxide = Ion('OH-', [-1], [-100], [-205E-9])

    _contents = OrderedDict()
    _pH = 7.                # Normal pH units.
    _ionic_strength = 0.    # Expected in molar.
    _temperature = reference_temperature  # Temperature in C

    @property
    def ions(self):
        return self._contents.keys()

    @property
    def concentrations(self):
        return np.array(self._contents.values())

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
            if isinstance(ion, basestring):
                ion = database.load(ion)
            else:
                ion = copy.copy(ion)
            ion.context(self)
            assert concentration >= 0, 'Concentrations must be positive.'
            self._contents[ion] = concentration

        if temperature is not None:
            self.temperature(temperature)
        self._equilibrate()

    def temperature(self, temperature=None):
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

    def cH(self, pH=None, ionic_strength=None):
        """Return the concentration of protons in solution."""
        if pH is None:
            pH = self.pH

        if ionic_strength is None:
            ionic_strength = self.ionic_strength

        cH = 10**(-pH)/self._hydronium.activity([1], ionic_strength)[0]
        return cH

    def cOH(self, pH=None, ionic_strength=None):
        """Return the concentration of hydroxyls in solution."""
        if pH is None:
            pH = self.pH

        if ionic_strength is None:
            ionic_strength = self.ionic_strength

        cOH = self.effective_dissociation(ionic_strength)/self.cH(pH)
        return cOH

    def concentration(self, ion):
        if ion in ('H+', self._hydronium):
            raise NotImplementedError
        elif ion in ('OH-', self._hydroxide):
            raise NotImplementedError
        else:
            return self._contents.get(ion, 0)

    def __add__(self, other):
        new_i = self.ions[:]
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
            raise NotImplementedError

    __radd__ = __add__

    def __mul__(self, other):
        if other >= 0:
            return Solution(self.ions,
                            [c * other for c in self.concentrations])
        else:
            raise NotImplementedError

    __rmul__ = __mul__

    def __str__(self):
        """Return a string representing the Solution."""
        return "Solution(pH={:.3g}, I={:.3g} M)".format(self.pH,
                                                        self.ionic_strength)

    def __repr__(self):
        """Return a representation of the Solution."""
        return self.__str__()

    def __len__(self):
        return len(self.ions)

    def serialize(self, nested=False):
        serial = {'__solution__': True}
        serial['concentrations'] = self.concentrations
        serial['ions'] = [ion.serialize(nested=True) for ion in self.ions]

        if nested:
            return serial
        else:
            return json.dumps(serial)

    def save(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.serialize(), file)

    from .equilibrium import _equilibrate, equilibrate_CO2, \
                            effective_dissociation
    from .conductivity import conductivity, hydroxide_conductivity, \
                              hydronium_conductivity
    from .titrate import titrate, buffering_capacity
    from .debye import debye
    from .transference import transference, zone_transfer
    from .conservation import kohlrausch, alberty, jovin, gas

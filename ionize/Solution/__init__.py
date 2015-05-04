from ..Ion import Ion
from math import log, log10, sqrt
from ..Aqueous import Aqueous
from ..load_ion import load_ion
import warnings
import json


class Solution(Aqueous):

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

    _F = 96485.3415        # Faraday's const.           [C/mol]
    _R = 8.31              # Universal gas const.       [J/mol*K]
    _Kw_ref = 1E-14        # Water equilibrium constant [mol^2]
    _Lpm3 = 1000.0         # Liters per meter^3         []
    _visc = 1E-3           # Dynamic viscosity (water)  [Pa s]
    _Adh = 0.512           # L^1/2 / mol^1/2, approximate for room temperature
    _aD = 1.5              # mol^-1/2 mol^-3/2, approximation
    _permittivity = 8.85e-12   # permittivity of free space
    _k = 1.38e-23              # Boltzman constant
    _Na = 6.02e23          # Avagadro's number
    _e = 1.6e-19               # elementary charge

    _H = Ion('H+', [1], [100], [362E-9])
    _OH = Ion('OH-', [-1], [-100], [-205E-9])

    ions = []              # Should be a list of ion objects.
    concentrations = []    # A list of concentrations in molar.
    pH = 7.0               # Normal pH units.
    I = 0.0                # Expected in molar.
    T = 25                 # Temperature in C
    _T_ref = 25            # reference temperature

    _dHw = 55.815e3
    _dCpw = -224

    def __init__(self, ions=[], concentrations=[], T=25):
        """Initialize a solution object."""

        self.T = float(T)
        if self.T == self._T_ref:
            self._Kw = self._Kw_ref
        else:
            self._Kw = self.adjust_Kw()

        if isinstance(ions, basestring):
            ions = [ions]
        try:
            self.ions = [i for i in ions]
        except:
            self.ions = [ions]

        for idx, ion in enumerate(self.ions):
            if isinstance(ion, basestring):
                self.ions[idx] = load_ion(ion)

        self.ions = [i.set_T(self.T) for i in self.ions]

        try:
            self.concentrations = [c for c in concentrations]
        except:
            self.concentrations = [concentrations]

        assert len(self.ions) == len(self.concentrations),\
            """Must be initialized with the same number of ions and concentrations.
        """

        assert all([c >= 0 for c in self.concentrations]),\
            """Concentrations must be positive."""

        if self.ions:
            (self.pH, self.I) = self._find_equilibrium()
        else:
            self.pH = -log10(sqrt(self._Kw))
            self.I = self._calc_I(self.pH)

        if self.I > 1:
            warnings.warn('Ionic strength > 1M. Ionic stregth correction may be inaccurate.')

        for ion in self.ions:
            ion._I = self.I
            ion._pH = self.pH

        actual_mobilities = self._onsager_fuoss()

        for i in range(len(self.ions)):
            self.ions[i].actual_mobility = actual_mobilities[i]

        self._H.actual_mobility = [actual_mobilities[-1][0]]
        self._OH.actual_mobility = [actual_mobilities[-1][1]]

    def set_T(self, T):
        return Solution(self.ions, self.concentrations, T=T)

    def adjust_Kw(self):
        pKw_ref = -log10(self._Kw_ref)
        T_ref = self._T_ref + 273.15
        T = self.T + 273.15
        pKw = pKw_ref -\
            (self._dHw/2.303/self._R)*(1.0/T_ref - 1.0/T) -\
            (self._dCpw/2.303/self._R)*(T_ref/T-1.0+log(T/T_ref))
        Kw = 10.0**(-pKw)
        return Kw

    def buffering_capacity(self):
        """Return the buffering capacity of the solution.

        This function generates an approximate solution to the buffering
        capacity by finding the derivative of the pH with respect to
        the addition of an acid insult at small concentration.
        """
        # Remove any ions at concentration 0.
        c = 0.001*min([cp for cp in self.concentrations if cp > 0])
        Cb = 0

        # Add an acid insult at 0.1% the lowest concentration in the solution.
        # If the buffering capacity is measured as above the insult c,
        # make the insult c lower.
        while Cb < c:
            new_sol = self + (Ion('Acid Insult', -1, -2, -1), c)
            # Find the slope of the pH.
            Cb = abs(c/(self.pH-new_sol.pH))
            c = 0.01 * Cb
        return Cb

    def cH(self, pH=None, I=None):
        """Return the concentration of protons in solution."""
        if not pH:
            pH = self.pH

        if not I:
            I = self.I

        cH = 10**(-pH)/self._H.activity_coefficient(I, [1])[0]
        return cH

    def cOH(self, pH=None, I=None):
        """Return the concentration of hydroxyls in solution."""
        if not pH:
            pH = self.pH

        if not I:
            I = self.I

        cOH = self.Kw_eff(I)/self.cH(pH)
        return cOH

    def H_conductivity(self):
        """Return the conductivity of protons in solution.

        Corrects for the mobility of the ion using the
        ion objects's actual mobility.
        """
        H_conductivity = self.cH()*self._H.molar_conductivity(self.pH, self.I)
        return H_conductivity

    def OH_conductivity(self):
        """Return the conductivity of hydroxyls in solution.

        Corrects for the mobility of the ion using the
        ion object's actual mobility.
        """
        OH_conductivity = self.cOH() *\
            self._OH.molar_conductivity(self.pH, self.I)

        return OH_conductivity

    def debye(self):
        """Return the Debye length of the solution.

        Uses the Debye-Huckel approximation for the calculation
        """
        dielectric = self._dielectric(self.T)
        viscosity = self._viscosity(self.T)
        epsilon = dielectric * self._permittivity
        lamda = (epsilon*self._k*self.T/self._e**2/self.I/self._Na)**.5
        return lamda

    def get_concentration(self, ion):
        if not ion in self.ions:
            return 0
        else:
            idx = self.ions.index(ion)
            return self.concentrations[idx]

    def __add__(self, other):
        new_i = self.ions[:]
        new_c = self.concentrations[:]
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
        return "Solution(pH={:.3g}, I={:.3g} M)".format(self.pH, self.I)

    def __repr__(self):
        """Return a representation of the Solution."""
        return self.__str__()

    def __len__(self):
        return len(self.ions)

    def serialize(self):
        serial = {'type': 'ionize solution',
                  'ions': [ion.serialize() for ion in self.ions],
                  'concentrations': self.concentrations,
                  'T': self.T
                  }

        return serial

    def save(self, filename):
        with open(filename, 'w') as file:
            json.dump(self.serialize(), file)

    from .calc_I import calc_I as _calc_I
    from .calc_pH import calc_pH as _calc_pH
    from .conductivity import conductivity
    from .equil_offset import equil_offset as _equil_offset
    from .find_equilibrium import find_equilibrium as _find_equilibrium
    from .Kw_eff import Kw_eff
    from .onsager_fuoss import onsager_fuoss as _onsager_fuoss
    from .transference import transference
    from .zone_transfer import zone_transfer
    from .conservation import kohlrausch, alberty, jovin, gas
    from .titrate import titrate

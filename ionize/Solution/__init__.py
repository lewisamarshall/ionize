import warnings
from ..Ion import Ion
import sys
from math import log, log10, sqrt


class Solution(object):

    """Represent a solution containing a set of ions.

    Initialize with solution(ions, concentrations), where ions is a list
    containing AspPy ion objects, and concentrations is a list containing
    positive scalars. Can also take T as an optional temperature argument, in
    Celcius.

    When a new solution is initialized, it will immediately calculate the
    equilibrium state, including the pH and the ionic strength (I) of the
    solution. These values willl be stored as permenant attributes of the
    object. Other solution properties can be calculated by invoking the
    appropriate method.

    See also Ion.
    """

    _F = 96485.3415        # Faraday's const.           [C/mol]
    _R = 8.31              # Universal gas const.       [J/mol*K]
    _Kw_ref = 1E-14        # Water equilibrium constant [mol^2]
    _Lpm3 = 1000.0         # Liters per meter^3         []
    _visc = 1E-3           # Dynamic viscosity (water)  [Pa s]
    _Adh = 0.512           # L^1/2 / mol^1/2, approximate for room temperature
    _aD = 1.5              # mol^-1/2 mol^-3/2, approximation

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

        try:
            self.ions = [i.set_T(self.T) for i in ions]
        except:
            self.ions = [ions.set_T(self.T)]

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
            (self.pH, self.I) = self.find_equilibrium()
        else:
            self.pH = -log10(sqrt(self._Kw))
            self.I = self.calc_I(self.pH)

        for ion in self.ions:
            ion._I = self.I
            ion._pH = self.pH

        actual_mobilities = self.onsager_fuoss()

        for i in range(len(self.ions)):
            self.ions[i].actual_mobility = actual_mobilities[i]

        self._H.actual_mobility = [actual_mobilities[-1][0]]
        self._OH.actual_mobility = [actual_mobilities[-1][1]]

    def adjust_Kw(obj):
        pKw_ref = -log10(obj._Kw_ref)
        T_ref = obj._T_ref + 273.15
        T = obj.T + 273.15
        pKw = pKw_ref -\
            (obj._dHw/2.303/obj._R)*(1.0/T_ref - 1.0/T) -\
            (obj._dCpw/2.303/obj._R)*(T_ref/T-1.0+log(T/T_ref))
        Kw = 10.0**(-pKw)
        return Kw

    def buffering_capacity(obj):
        """Return the buffering capacity of the solution.

        This function generates an approximate solution to the buffering
        capacity by finding the derivative of the pH with respect to
        the addition of an acid insult at small concentration.
        """
        # Remove any ions at concentration 0.
        c = 0.001*min([cp for cp in obj.concentrations if cp > 0])
        Cb = 0

        # Add an acid insult at 0.1% the lowest concentration in the solution.
        # If the buffering capacity is measured as above the insult c,
        # make the insult c lower.
        while Cb < c:
            new_sol = obj + (Ion('Acid Insult', -1, -2, -1), c)
            # Find the slope of the pH.
            Cb = abs(c/(obj.pH-new_sol.pH))
            c = 0.01 * Cb
        return Cb

    def cH(obj, pH=None, I=None):
        """Return the concentration of protons in solution."""
        if not pH:
            pH = obj.pH

        if not I:
            I = obj.I

        cH = 10**(-pH)/obj._H.activity_coefficient(I, [1])[0]
        return cH

    def cOH(obj, pH=None, I=None):
        """Return the concentration of hydroxyls in solution."""
        if not pH:
            pH = obj.pH

        if not I:
            I = obj.I

        cOH = obj.Kw_eff(I)/obj.cH(pH)
        return cOH

    def H_conductivity(obj):
        """Return the conductivity of protons in solution.

        Corrects for the mobility of the ion using the
        ion object's actual mobility.
        """
        H_conductivity = obj.cH()*obj._H.molar_conductivity(obj.pH, obj.I)
        return H_conductivity

    def OH_conductivity(obj):
        """Return the conductivity of hydroxyls in solution.

        Corrects for the mobility of the ion using the
        ion object's actual mobility.
        """
        OH_conductivity = obj.cOH()*obj._OH.molar_conductivity(obj.pH, obj.I)
        return OH_conductivity

    def __add__(obj, other):
        new_i = obj.ions[:]
        new_c = obj.concentrations[:]
        if isinstance(other, Solution):
            for ion, c in zip(other.ions, other.concentrations):
                if ion in obj.ions:
                    new_c[obj.ions.index(ion)] += c
                else:
                    new_i.append(ion)
                    new_c.append(c)
            return Solution(new_i, new_c)
        elif isinstance(other, (list, tuple)) and len(other) == 2 and\
                isinstance(other[0], Ion) and isinstance(other[1], (int, float)):
            ion, c = other
            if ion in obj.ions:
                new_c[obj.ions.index(ion)] += c
            else:
                new_i.append(ion)
                new_c.append(c)
            return Solution(new_i, new_c)
        else:
            raise NotImplementedError

    __radd__ = __add__

    def __mul__(obj, other):
        if other >= 0:
            return Solution(obj.ions,
                            [c * other for c in obj.concentrations])
        else:
            raise NotImplementedError

    __rmul__ = __mul__

    def __str__(obj):
        """Return a string representing the Solution."""
        return "Solution object -- pH = " + str(obj.pH) + \
            ", I = " + str(obj.I) + ' M'

    def __repr__(obj):
        """Return a representation of the Solution."""
        return obj.__str__()

    def __len__(obj):
        return len(obj.ions)

    from .calc_I import calc_I
    from .calc_pH import calc_pH
    from .conductivity import conductivity
    from .equil_offset import equil_offset
    from .find_equilibrium import find_equilibrium
    from .Kw_eff import Kw_eff
    from .onsager_fuoss import onsager_fuoss
    from .transference import transference
    from .zone_transfer import zone_transfer
    from ..dielectric import dielectric
    from ..viscosity import viscosity
    from .conservation import kohlrausch, alberty, jovin

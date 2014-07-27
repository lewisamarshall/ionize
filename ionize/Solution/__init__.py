import warnings
from ..Ion import Ion
import sys


class Solution(object):

    """Represent a solution containing a set of ions.

    Initialize with solution(ions, concentrations), where ions is a list
    containing AspPy ion objects, and concentrations is a list containing
    positive scalars.

    When a new solution is initialized, it will immediately calculate the
    equilibrium state, including the pH and the ionic strength (I) of the
    solution. These values willl be stored as permenant attributes of the
    object. Other solution properties can be calculated by invoking the
    appropriate method.

    See also Ion.
    """

    _F = 96485.3415        # Faraday's const.           [C/mol]
    _Rmu = 8.31            # Universal gas const.       [J/mol*K]
    _Kw = 1E-14            # Water equilibrium constant [mol^2]
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

    def __init__(self, ions=[], concentrations=[], T=25):
        """Initialize a solution object."""
        self.T = T
        try:
            self.ions = [i for i in ions]
        except:
            self.ions = [ions]

        try:
            self.concentrations = [c for c in concentrations]
        except:
            self.concentrations = [concentrations]

        assert len(self.ions) == len(self.concentrations),\
            """Must be initialized with the same number of ions and concentrations.
        """

        assert all([c >= 0 for c in concentrations]),\
            """Concentrations must be positive."""

        if self.ions:
            (self.pH, self.I) = self.find_equilibrium()
        else:
            self.I = self.calc_I(self.pH)

        actual_mobilities = self.onsager_fuoss()

        for i in range(len(self.ions)):
            self.ions[i].actual_mobility = actual_mobilities[i]

        self._H.actual_mobility = [actual_mobilities[-1][0]]
        self._OH.actual_mobility = [actual_mobilities[-1][1]]

    def add_ion(obj, new_ions, new_concentrations):
        """add_ion initializes a new solution with additional ions.

        The returned solution will contain all of the ions in the current
        solution plus a new set of ions from new_ions at a new set of
        concentrations from new_concentrations.
        """
        new_solution = Solution(obj.ions + new_ions,
                                obj.concentrations + new_concentrations)
        return new_solution

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
        if isinstance(other, Solution):
            return Solution(obj.ions + other.ions,
                            obj.concentrations + other.concentrations)
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

    from buffering_capacity import buffering_capacity
    from calc_I import calc_I
    from calc_pH import calc_pH
    from conductivity import conductivity
    from equil_offset import equil_offset
    from find_equilibrium import find_equilibrium
    from Kw_eff import Kw_eff
    from onsager_fuoss import onsager_fuoss
    from transference import transference
    from zone_transfer import zone_transfer
    from ..dielectric import dielectric
    from ..viscosity import viscosity

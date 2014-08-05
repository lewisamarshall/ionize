import warnings
from math import sqrt, copysign


class Ion(object):

    """Describe an ion dissolved in aqueous solution.

    Initialize with Ion(name, z, pKa, absolute_mobility).
    """
    # Weakly private variables
    # These are constants and should not change.
    # Eventually, _T may be  removed from the constants list.
    _F = 96485.34         # Faraday's const.[C/mol]
    _Lpm3 = 1000.0        # Conversion from liters to m^3

    # The following are constants in eqtn 6 of Bahga 2010.
    _Adh = 0.5102         # L^1/2 / mol^1/2, approximate for RT
    _aD = 1.5             # mol^-1/2 mol^-3/2, approximation
    R = 8.314             # J/mol-K
    # The reference properties of the ion are stored in private variables.
    _pKa_ref = []
    _absolute_mobility_ref = []  # m^2/V/s.
    dH = None
    dCp = None
    _pH = None
    _I = None

    def __init__(self, name, z, pKa_ref, absolute_mobility_ref,
                 dH=None, dCp=None, nightingale_function=None,
                 T=25.0, T_ref=25.0):
        """Initialize an ion object."""
        self.name = name
        self.T = T
        self._T_ref = T_ref
        self.dH = dH
        self.dCp = dCp
        self.nightingale_function = nightingale_function
        try:
            self.z = [zp for zp in z]
        except:
            self.z = [z]

        try:
            self._pKa_ref = [p for p in pKa_ref]
        except:
            self._pKa_ref = [pKa_ref]

        try:
            self._absolute_mobility_ref = [m for m in absolute_mobility_ref]
        except:
            self._absolute_mobility_ref = [absolute_mobility_ref]

        if T == T_ref:
            self.pKa = self._pKa_ref
            self.absolute_mobility = self._absolute_mobility_ref
        else:
            _Adh = self.get_Adh()
            self.pKa = self.correct_pKa()
            if self.nightingale_function:
                self.absolute_mobility = [self.nightingale_function(self.T).tolist()] *\
                    len(self.z)
            else:
                self.absolute_mobility =\
                    [self._viscosity(self._T_ref)/self._viscosity(self.T)*m
                     for m in self._absolute_mobility_ref]

        self.actual_mobility = None                 # Fill by solution

        # Check that z is a vector of integers
        assert all([isinstance(zp, int) for zp in self.z]), \
            "z contains non-integer"

        # Check that the pKa is a vector of numbers of the same length as z.
        assert len(self.pKa) == len(self.z), "pKa is not the same length as z"

        assert len(self.absolute_mobility) == len(self.z), '''absolute_mobility is not
                                                    the same length as z'''

        # Force the sign of the fully ionized mobilities to match the sign of
        # the charge. This command provides a warning.
        if not all([copysign(z, m) == z for z, m in zip(self.z,
                    self.absolute_mobility)]):
            self.absolute_mobility = [copysign(m, z) for z, m in zip(self.z,
                                      self.absolute_mobility)]
            warnings.warn("Mobility signs and charge signs don't match.")

        # After storing the ion properties, ensure that the properties are
        # sorted in order of charge. All other ion methods assume that the
        # states will be sorted by charge.
        self = self.z_sort()
        self.Ka = self.get_Ka()
        self.z0 = self.get_z0()

    def z_sort(obj):
        """Sort the charge states from lowest to highest."""
        # Zip the lists together and sort them by z.
        obj.z, obj.pKa, obj.absolute_mobility = zip(*sorted(zip(obj.z, obj.pKa,
                                                    obj.absolute_mobility)))
        obj.z = list(obj.z)
        obj.pKa = list(obj.pKa)
        obj.absolute_mobility = list(obj.absolute_mobility)

        full = set(range(min(obj.z), max(obj.z)+1, 1)) - {0}
        assert set(obj.z) ^ full == set(), "Charge states missing."

        return obj

    def get_Ka(obj):
        """Return the Kas based on the pKas.

        These values are not corrected for ionic strength.
        """
        Ka = [10.**-p for p in obj.pKa]
        return Ka

    def get_z0(obj):
        """Return the list of charge states with 0 inserted."""
        z0 = [0]+obj.z
        z0 = sorted(z0)
        return z0

    from ..dielectric import dielectric

    def get_Adh(obj, T=None):
        """Account for the temperature dependance of Adh."""
        if not T:
            T = obj.T
        T_ref = 25
        Adh_ref = 0.5102
        d = obj.dielectric(T)
        d_ref = obj.dielectric(T)
        Adh = Adh_ref * ((T_ref+273.15)*d_ref/(T+273.15)/d)**(-1.5)
        return Adh

    def set_T(obj, T):
        """Return a new ion at the specified temperature."""
        return Ion(obj.name, obj.z, obj._pKa_ref, obj._absolute_mobility_ref,
                   obj.dH, obj.dCp, obj.nightingale_function,
                   T=T, T_ref=obj._T_ref)

    def __str__(obj):
        """Return a string representing the ion."""
        return ("Ion object -- " + obj.name + ": " +
                str(dict(zip(obj.z, zip(obj.pKa, obj.absolute_mobility)))))

    def __repr__(obj):
        """Return a representation of the ion."""
        return obj.__str__()

    def __eq__(self, other):
        if self.name == other.name and\
                self.z == other.z and \
                self._pKa_ref == other._pKa_ref and\
                self._absolute_mobility_ref == other._absolute_mobility_ref and\
                self.dH == other.dH and\
                self.dCp == other.dCp and\
                self._T_ref == other._T_ref:
                # self.nightingale_function == other.nightingale_function and\
            return True
        else:
            return False

    from .ionization_fraction import ionization_fraction
    from .activity_coefficient import activity_coefficient
    from .effective_mobility import effective_mobility
    from .Ka_eff import Ka_eff
    from .L import L
    from .molar_conductivity import molar_conductivity
    from .robinson_stokes_mobility import robinson_stokes_mobility
    from .correct_pKa import correct_pKa, vant_hoff, clark_glew
    from ..viscosity import viscosity as _viscosity

if __name__ == '__main__':
    pass

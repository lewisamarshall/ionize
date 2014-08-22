import warnings
from math import copysign


class Ion(object):

    """Describe an ion dissolved in aqueous solution.

    Args:
        name (str): The chemical name of the ion.

        z (list): A list of valence states for the ion, as integers.

        pKa_ref (list): The pKa of each valence at the refernce temperature,\
        as floats.

        absolute_mobility_ref (list): The signed absolute mobility of each\
        valence at the reference temperature, as floats, in units\
        of m^2/V/s. Expect O(10^-8).

        dH (list): The enthalpy of dissociation of each valence, at the\
        reference temperature, as floats.

        dCp (list): The change in heat capacity of dissociation of each\
        valence, at the reference temperature, as floats.

        nightingale_function (function): A function describing absolute\
        mobility as a function of temperature, for special ions.

        T (float): The temperature to use to calculate the properties of the
        ions, in degrees C.

        T_ref (float): The reference temperature for the reference properties,
        in degrees C.

    Attributes:
        z (list): A list of valence states for the ion, as integers.

        pKa (list): The pKa of each valence at the refernce temperature,\
            as floats.

        absolute_mobility (list): The signed absolute mobility of each\
            valence at the reference temperature, as floats, in units\
            of m^2/V/s. Expect O(10^-8).

        T (float): The temperature to use to calculate the properties of the\
            ions, in degrees C.

    Raises:
        None

    Example:
        To to initialize an Ion, call as:

        >>> ionize.Ion('my_acid', [-1, -2], [1.2, 3.4], [-10e-8, -21e-8])
    """
    # These are constants and should not change.
    _F = 96485.34         # Faraday's const.[C/mol]
    _Lpm3 = 1000.0        # Conversion from liters to m^3
    _R = 8.314             # J/mol-K

    # The following are constants in eqtn 6 of Bahga 2010.
    # _Adh is updated for temperature on initialize.
    _Adh = 0.5102         # L^1/2 / mol^1/2, approximate for RT
    # _aD is treated as a constant, though it does vary slightly with temp.
    _aD = 1.5             # L^3/2 mol^-1

    # The reference properties of the ion are stored and used to calculate
    # properties at the current temperature.
    _pKa_ref = []
    _absolute_mobility_ref = []  # m^2/V/s.
    _T_ref = 25

    # The properties of the ions are stored in public variables.
    # These are the properties at the current temperature, or are treated
    # as temperature independant.
    pKa = None
    Ka = None
    absolute_mobility = None
    dH = None
    dCp = None
    z0 = None
    T = 25

    # If the Ion is in a solution object, copy the pH and I of the Solution
    # locally for reference, in a private variable. Also store the Onsager-
    # Fouss mobility in actual_mobility.
    _pH = None
    _I = None
    actual_mobility = None

    def __init__(self, name, z, pKa_ref, absolute_mobility_ref,
                 dH=None, dCp=None, nightingale_function=None,
                 T=25.0, T_ref=25.0):
        """Initialize an Ion object."""
        # Copy properties into the ion.
        self.name = name
        self.T = T
        self._T_ref = T_ref
        self.dH = dH
        self.dCp = dCp
        self.nightingale_function = nightingale_function

        # Copy in the properties that should be lists, as long as they are
        # single values or iterables.
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

        # Temperature adjust the ion.
        self._set_Adh()
        if T == T_ref:
            self.pKa = self._pKa_ref
            self.absolute_mobility = self._absolute_mobility_ref
        else:
            self.pKa = self._correct_pKa()
            if self.nightingale_function:
                self.absolute_mobility = [self.nightingale_function(self.T).tolist()] *\
                    len(self.z)
            else:
                self.absolute_mobility =\
                    [self._viscosity(self._T_ref)/self._viscosity(self.T)*m
                     for m in self._absolute_mobility_ref]

        # Force the sign of the fully ionized mobilities to match the sign of
        # the charge. This command provides a warning.
        if not all([copysign(z, m) == z for z, m in zip(self.z,
                    self.absolute_mobility)]):
            self.absolute_mobility = [copysign(m, z) for z, m in zip(self.z,
                                      self.absolute_mobility)]
            warnings.warn("Mobility signs and charge signs don't match.")

        # Check that z is a vector of integers
        assert all([isinstance(zp, int) for zp in self.z]), \
            "z contains non-integer"

        # Check that the pKa is a vector of numbers of the same length as z.
        assert len(self.pKa) == len(self.z), "pKa is not the same length as z"

        assert len(self.absolute_mobility) == len(self.z), '''absolute_mobility is not
                                                    the same length as z'''

        # After storing the ion properties, ensure that the properties are
        # sorted in order of charge. All other ion methods assume that the
        # states will be sorted by charge.
        self._z_sort()
        self._set_Ka()
        self._set_z0()

    def _z_sort(self):
        """Sort the charge states from lowest to highest."""
        # Zip the lists together and sort them by z.
        self.z, self.pKa, self.absolute_mobility =\
            zip(*sorted(zip(self.z, self.pKa, self.absolute_mobility)))
        self.z = list(self.z)
        self.pKa = list(self.pKa)
        self.absolute_mobility = list(self.absolute_mobility)

        full = set(range(min(self.z), max(self.z)+1, 1)) - {0}
        assert set(self.z) ^ full == set(), "Charge states missing."

        return None

    def _set_Ka(self):
        """Set the Kas based on the pKas.

        These values are not corrected for ionic strength.
        """
        self.Ka = [10.**-p for p in self.pKa]
        return None

    def _set_z0(self):
        """Set the list of charge states with 0 inserted."""
        self.z0 = sorted([0]+self.z)
        return None

    def _set_Adh(self, T=None):
        """Account for the temperature dependance of Adh."""
        if not T:
            T = self.T
        T_ref = 25
        Adh_ref = 0.5102
        d = self._dielectric(T)
        d_ref = self._dielectric(T)
        self._Adh = Adh_ref * ((T_ref+273.15)*d_ref/(T+273.15)/d)**(-1.5)
        return None

    def set_T(self, T):
        """Return a new ion at the specified temperature."""
        return Ion(self.name, self.z, self._pKa_ref,
                   self._absolute_mobility_ref, self.dH, self.dCp,
                   self.nightingale_function,
                   T=T, T_ref=self._T_ref)

    def __str__(self):
        """Return a string representing the ion."""
        obj_str = "Ion('{}', z={})".format(self.name, self.z)
        return obj_str

    def __repr__(self):
        """Return a representation of the ion."""
        return self.__str__()

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

    from ..dielectric import dielectric as _dielectric
    from .ionization_fraction import ionization_fraction
    from .activity_coefficient import activity_coefficient
    from .effective_mobility import effective_mobility
    from .Ka_eff import Ka_eff
    from .L import L
    from .molar_conductivity import molar_conductivity
    from .robinson_stokes_mobility import robinson_stokes_mobility
    from .correct_pKa import _correct_pKa, _vant_hoff, _clark_glew
    from ..viscosity import viscosity as _viscosity

if __name__ == '__main__':
    pass

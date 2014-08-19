def equil_offset(self, I_i):
    """Return the error in ionic strength.

    Takes an ionic strength, then uses it to calculate a new ionic strength
    using the new equilibrum coefficents. find_equilibrium finds the root of
    this function.
    """
    pH = self._calc_pH(I_i)
    I_f = self._calc_I(pH, I_i)

    res = (I_f-I_i)
    return res

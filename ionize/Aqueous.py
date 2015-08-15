"""Create the Aqueous class to hold the properties of water."""


class Aqueous(object):

    """Describe the properties of water."""

    _Kw_ref = 1E-14        # Water equilibrium constant [mol^2]
    _dHw = 55.815e3
    _dCpw = -224
    _T_ref = 25.
    _R = 8.31              # Universal gas const.       [J/mol*K]

    def dielectric(self, T=None):
        """Return the dielectric constant of water at a specified temperature.

        The temperature should be specified in Celcius. Correlation is based on
        the CRC handbook.
        """
        Tk = T+273.15  # Convert to Kalvin.
        d = 249.21 - 0.79069*Tk + 0.72997e-3*Tk**2
        return d

    def viscosity(self, T=None):
        """Return the viscosity of water at the specified temperature.

        The temperature should be specified in celcius. The function internally
        converts the temperature to Kelvin. Correlation is based on Fox and
        McDonald's Intro to FLuid Mechanics, as implimented in STEEP.
        """
        Tk = T+273.15
        v = (2.414e-5)*10**(247.8/(Tk-140))
        return v

    def dissociation(self, T=None):
        """Return the dissociation constant of water."""
        if not T:
            return self._Kw_ref

        pKw_ref = -log10(self._Kw_ref)
        T_ref = self._T_ref + 273.15
        T += 273.15
        pKw = pKw_ref -\
            (self._dHw/2.303/self._R)*(1.0/T_ref - 1.0/T) -\
            (self._dCpw/2.303/self._R)*(T_ref/T-1.0+log(T/T_ref))
        Kw = 10.0**(-pKw)
        return Kw

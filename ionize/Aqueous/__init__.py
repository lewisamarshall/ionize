"""Create the Aqueous class to hold the properties of water."""


class Aqueous(object):

    """Describe the properties of water."""

    def __init__(self):
        """Initialize."""
        pass

    def _dielectric(self, T=None):
        """Return the dielectric constant of water at a specified temperature.

        The temperature should be specified in Celcius. Correlation is based on
        the CRC handbook.
        """
        Tk = T+273.15  # Convert to Kalvin.
        d = 249.21 - 0.79069*Tk + 0.72997e-3*Tk**2
        return d

    def _viscosity(self, T=None):
        """Return the viscosity of water at the specified temperature.

        The temperature should be specified in celcius. The function internally
        converts the temperature to Kelvin. Correlation is based on Fox and
        McDonald's Intro to FLuid Mechanics, as implimented in STEEP.
        """
        Tk = T+273.15
        v = (2.414e-5)*10**(247.8/(Tk-140))
        return v

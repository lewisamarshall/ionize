from math import sqrt


def activity_coefficient(self, I=None, z=None):
    """Return activity coefficients of a charge state at ionic strength I."""
    if I is None:
        if self._I:
            I = self._I
        else:
            I = 0.0

    if not z:
        z = self.z0
    else:
        try:
            z = [zp for zp in z]
        except:
            z = [z]

    # There are two coefficients that are used repeatedly.
    # Specified in Bahga.
    A = self._Adh*sqrt(I)/(1.0+self._aD*sqrt(I))
    B = 0.1*I  # Matching STEEP implementation.

    # Use them to calculate the activity coefficients.
    # These coefficients are for z=+-1, for H+ and OH-
    gamma = [10**(zp**2*(B-A)) for zp in z]

    return gamma

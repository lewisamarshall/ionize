from math import sqrt


def activity_coefficient(obj, I, z=None):
    """Return activity coefficients of a charge state at ionic strength I."""
    if not z:
        z = obj.z0()
    else:
        try:
            z = [zp for zp in z]
        except:
            z = [z]

    # There are two coefficients that are used repeatedly.
    # Specified in Bahga.
    A = obj._Adh*sqrt(I)/(1+obj._aD*sqrt(I))
    B = 0.1*I  # Matching STEEP implementation.
    # Use them to calculate the activity coefficients.
    # These coefficients are for z=+-1, for H+ and OH-
    gamma = [10**(zp**2*(B-A)) for zp in z]
    return gamma

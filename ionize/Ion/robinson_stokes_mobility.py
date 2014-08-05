from math import copysign, sqrt
import warnings


def robinson_stokes_mobility(obj, I=None, T=25):
    """Return the Robinson-Stokes correction to fully ionized mobility.

    This correction is appropriate if a generic ionic strength is known,
    but the specific ions in solution are unknown.
    """
    # Currently using the ionic strength where Bahga 2010
    # uses twice the ionic strength. This appears to work, and follows the
    # SPRESSO implimentation.
    # Likely typo in paper.
    if I is None:
        if obj._I:
            I = obj._I
        else:
            I = 0

    if not T:
        T = obj.T
    T_ref = 25
    d = obj._dielectric(T)
    d_ref = obj._dielectric(T)

    A = 0.2297*((T_ref+273.15)*d_ref/(T+273.15)/d)**(-1.5)
    B = 31.410e-9 * ((T_ref+273.15)*d_ref/(T+273.15)/d)**(-0.5) *\
        obj._viscosity(T_ref)/obj._viscosity(T)
    actual_mobility = []
    for abs_mob, z in zip(obj.absolute_mobility, obj.z):
        actual_mobility.append(abs_mob -
                               (A * abs_mob +
                                B * copysign(1, z)) * sqrt(I) /
                               (1 + obj._aD * sqrt(I)))

    return actual_mobility

from math import copysign, sqrt


def robinson_stokes_mobility(obj, I):
    """Return the robinson stokes correction to fully ionized mobility.

    If only an ionic strength is specified, use the Robinson-Stokes
    correction to calculate a new fully ionized mobility.

    If a solution object is supplied, use the full onsager fouss correction.
    """
    if I >= 0:
        # Currently using the ionic strength where Bahga 2010
        # uses twice the ionic strength. This appears to work, and follows the
        # SPRESSO implimentation.
        # Likely typo in paper.
        A = 0.2297
        B = 31.410e-9
        actual_mobility = []
        for abs_mob, z in zip(obj.absolute_mobility, obj.z):
            actual_mobility.append(abs_mob -
                                   (A * abs_mob +
                                    B * copysign(1, z)) * sqrt(I) /
                                   (1 + obj._aD * sqrt(I)))
        else:
            error('''Ionic strength must be specified as
                    a scalar positive value.''')

    return actual_mobility

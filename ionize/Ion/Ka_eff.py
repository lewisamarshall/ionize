def Ka_eff(obj, I=None):
    """Return the effective Ka values for the ion.

    Args:
        I (float): The ambiant ionic strength.

    This function correct the Ka for ionic strength, using the Dubye-Huckle
    theory to calculate activity coefficients. If no ionic strength is supplied,
    and the Ion is nested in a Solution, the solution ionic strength will be
    used. Otherwise, the ionic strength is assumed to be 0.
    """
    if I is None:
        if obj._I:
            I = obj._I
        else:
            I = 0.0

    # If the ionic strength is zero, simply return the Ka's.
    if I is 0:
        return obj.Ka[:]

    # Make the effective Ka vector the same size as the Ka vector.
    Ka_eff = []

    gam_i = obj.activity_coefficient(I)
    gam_h = obj.activity_coefficient(I, [1])

    # For each acidity coefficient, get the effective
    # coefficient by multiplying by activities.
    for i, Kp in enumerate(obj.Ka):
        Ka_eff.append(Kp*gam_i[i+1]/gam_i[i]/gam_h[0])

    return Ka_eff

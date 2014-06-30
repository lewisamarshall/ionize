def Ka_eff(obj, I=0):
    """Return the effective Ka values for the ion.

    This function uses the ionic strength correction function from
    Dubye-Huckle theory to calculate the activity coefficients, and uses
    these to correct Ka.
    """
    # If the ionic strength is zero, simply return the Ka's.
    if I is 0:
        return obj.Ka()

    # Make the effective Ka vector the same size as the Ka vector.
    Ka_eff = obj.Ka()

    gam_i = obj.activity_coefficient(I)
    gam_h = obj.activity_coefficient(I, 1)

    # For each acidity coefficient, get the effective
    # coefficient by multiplying by activities.
    for i, Kp in enumerate(obj.Ka):
        Ka_eff[i] = Kp*gam_i[i+1]/gam_i[i]/gam_h
    return Ka_eff

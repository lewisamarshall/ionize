def effective_mobility(obj, pH, I=0):
    """Return the effective mobility of the ion at a given pH and I.

    If an actual mobility from Onsager-Fouss is available, it is used,
    otherwise, the Robinson-Stokes correction is used.
    """
    if obj.actual_mobility:
        actual_mobility = obj.actual_mobility
    else:
        actual_mobility = obj.robinson_stokes_mobility(I)

    i_frac = obj.ionization_fraction(pH, I)
    effective_mobility = sum([f*m for (f, m) in zip(i_frac, actual_mobility)])

    return effective_mobility

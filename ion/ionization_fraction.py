def ionization_fraction(obj, pH, I=0):
    """Return the ionization fractions of an ion.

    The ionization fraction is based on a given pH and ionic strength.
    """
    # Get the vector of products of acidity constants.
    L = obj.L(I)
    # Compute the concentration of H+ from the pH.
    cH = 10**(-pH)/obj.activity_coefficient(I, [1])[0]

    # Calculate the numerator of the function for ionization fraction.
    i_frac_vector = [Lp * cH ** z for (Lp, z) in zip(L, obj.z0())]

    # Calculate the vector of ionization fractions
    i_frac = [i/sum(i_frac_vector) for i in i_frac_vector]

    # filter out the neutral fraction.
    i_frac = [i for (i, z) in zip(i_frac, obj.z0()) if z]

    return i_frac

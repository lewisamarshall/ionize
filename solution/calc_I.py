def calc_I(obj, pH, I_guess=0):
    """Return an estimate of the ionic strength of the solution.

    calc_I should only be used during solution initialization.
    After that, the equilibrium ionic strength is stored in obj.I.

    pH must be supplied. If a guess for ionic strength is not supplied,
    ionic strength corrections will not be used.
    """
    # For each ion, add the contribution to ionic strength to the sum.
    I = sum([c * sum(
        [z**2*f/2 for (z, f) in
         zip(ionp.z, ionp.ionization_fraction(pH, I_guess))])
        for c, ionp in zip(obj.concentrations, obj.ions)])

    # Add the ionic strength due to water dissociation.
    I = I+(obj.cH(pH) + obj.cOH(pH, I_guess))/2
    return I

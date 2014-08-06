def molar_conductivity(obj, pH=None, I=None):
    """Retun the molar conductivity of the ion based on the pH and I.

    Provides conducitivity in Siemens per meter per mole.

    Args:
        pH (float): The ambiant pH.
        I (float): The ambiant ionic strength.

    If an actual mobility from Onsager-Fouss is available, it is used,
    otherwise, the Robinson-Stokes mobility estimate is used.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].molar_conductivity()
    
    Otherwise, always call with a pH argument.
    """
    if pH is None:
        assert obj._pH, 'requires an input pH'
        pH = obj._pH

    if I is None:
        if obj._I:
            I = obj._I
        else:
            I = 0

    if obj.actual_mobility:
        actual_mobility = obj.actual_mobility
    else:
        actual_mobility = obj.robinson_stokes_mobility(I)

    i_frac = obj.ionization_fraction(pH, I)

    m_conductivity = (obj._Lpm3 * obj._F *
                      sum(z * f * m for (z, f, m)
                          in zip(obj.z, i_frac, actual_mobility)))

    return m_conductivity

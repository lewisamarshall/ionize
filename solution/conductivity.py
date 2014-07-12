def conductivity(obj):
    """Return the electrical conductivity of the solution.

    Relies on molar conductivity calculations from ion and total conductivity
    of H+ and OH-.
    """
    con = 0
    for c, i in zip(obj.concentrations, obj.ions):
        con += c * i.molar_conductivity(obj.pH, obj.I)
        print c *i.molar_conductivity(obj.pH, obj.I)

    con += obj.OH_conductivity()    # Add OH- contribution.
    con += obj.H_conductivity()     # Add H+ contribution.
    return con

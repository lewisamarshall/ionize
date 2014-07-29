def conductivity(obj):
    """Return the electrical conductivity of the solution.

    Relies on molar conductivity calculations from ion and total conductivity
    of H+ and OH-.

    The conductivity is returned in Seimens per meter.
    """
    con = 0
    for c, i in zip(obj.concentrations, obj.ions):
        con += c * i.molar_conductivity()

    con += obj.OH_conductivity()    # Add OH- contribution.
    con += obj.H_conductivity()     # Add H+ contribution.
    return con

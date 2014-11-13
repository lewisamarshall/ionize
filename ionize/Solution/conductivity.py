def conductivity(self):
    """Return the electrical conductivity of the solution.

    Relies on molar conductivity calculations from ion and total conductivity
    of H+ and OH-.

    The conductivity is returned in Seimens per meter.
    """
    conduct = 0
    for c, i in zip(self.concentrations, self.ions):
        conduct += c * i.molar_conductivity()

    conduct += self.OH_conductivity()    # Add OH- contribution.
    conduct += self.H_conductivity()     # Add H+ contribution.
    return conduct

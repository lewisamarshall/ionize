def conductivity(self):
    """Return the electrical conductivity of the solution.

    Relies on molar conductivity calculations from ion and total conductivity
    of H+ and OH-.

    The conductivity is returned in Seimens per meter.
    """
    conduct = 0
    for c, i in zip(self.concentrations, self.ions):
        conduct += c * i.molar_conductivity()

    conduct += self.hydronium_conductivity() + self.hydroxide_conductivity()
    return conduct


def hydronium_conductivity(self):
    """Return the conductivity of protons in solution."""
    H_conductivity = self.cH() * \
        self._hydronium.molar_conductivity(self.pH, self.ionic_strength)
    return H_conductivity


def hydroxide_conductivity(self):
    """Return the conductivity of hydroxyls in solution."""
    OH_conductivity = self.cOH() *\
        self._hydroxide.molar_conductivity(self.pH, self.ionic_strength)

    return OH_conductivity

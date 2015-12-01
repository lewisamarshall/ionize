def conductivity(self):
    """Return the electrical conductivity of the solution, in Seimens/meter.
    """
    conduct = 0
    for c, i in zip(self.concentrations, self.ions):
        conduct += c * i.molar_conductivity()

    conduct += self.hydronium_conductivity() + self.hydroxide_conductivity()
    return conduct


def hydronium_conductivity(self):
    """Return the conductivity of protons in solution, in Seimens/meter."""
    H_conductivity = self.concentration('H+') * \
        self._hydronium.molar_conductivity(self.pH,
                                           self.ionic_strength,
                                           self.temperature())
    return H_conductivity


def hydroxide_conductivity(self):
    """Return the conductivity of hydroxyls in solution, in Seimens/meter."""
    OH_conductivity = self.concentration('OH-') *\
        self._hydroxide.molar_conductivity(self.pH,
                                           self.ionic_strength,
                                           self.temperature())
    return OH_conductivity

def transference(self):
    """Return the fraction of charge carried by each of the ions as a list.

    Should not precisely add to 1, because some charge is carried by protons
    and hydroxyls.
    """
    T = [0] * len(self.ions)
    for i in range(len(T)):
        T[i] = (self.ions[i].molar_conductivity(self.pH, self.I) *
                self.concentrations[i])

    T = [tp / self.conductivity() for tp in T]
    return T

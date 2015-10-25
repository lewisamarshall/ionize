from ..constants import lpm3

def transference(self):
    """Return the fraction of charge carried by each of the ions as a list.

    Should not precisely add to 1, because some charge is carried by protons
    and hydroxyls.
    """
    T = [0] * len(self.ions)
    for i in range(len(T)):
        T[i] = (self.ions[i].molar_conductivity(self.pH, self.ionic_strength) *
                self.concentrations[i])

    T = [tp / self.conductivity() for tp in T]
    return T

def zone_transfer(self):
    """Return the zone transfer charge of the solution per liter."""
    Qi = [0]*len(self.ions)
    transference = self.transference()
    for i in range(len(Qi)):
        if self.concentrations[i] == 0:
            Qi[i] = 0
        else:
            Qi[i] = (self.ions[i].molar_conductivity(self.pH,
                                                     self.ionic_strength) *
                     self.concentrations[i]/transference[i] /
                     abs(self.ions[i].mobility(self.pH,
                                               self.ionic_strength)))

    Qi = [Qp/lpm3 for Qp in Qi]
    return Qi

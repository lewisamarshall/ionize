def zone_transfer(self):
    """Return the zone transfer charge of the solution per liter."""
    Qi = [0]*len(self.ions)
    transference = self.transference()
    for i in range(len(Qi)):
        if self.concentrations[i] == 0:
            Q[i] = 0
        else:
            Qi[i] = (self.ions[i].molar_conductivity(self.pH, self.I) *
                     self.concentrations[i]/transference[i] /
                     abs(self.ions[i].effective_mobility(self.pH, self.I)))

    Qi = [Qp/self._Lpm3 for Qp in Qi]
    return Qi

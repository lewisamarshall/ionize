def zone_transfer(self, vol):
    """Return the zone transfer charge of the solution at a given volume.

    The volume of the solution is specified in liters.
    """
    Qi = [0]*len(self.ions)
    transference = self.transference()
    for i in range(len(Qi)):
        Qi[i] = (self.ions[i].molar_conductivity(self.pH, self.I) *
                 self.concentrations[i]/transference[i] /
                 abs(self.ions[i].effective_mobility(self.pH, self.I)))
    Qi = [Qp*vol/self._Lpm3 for Qp in Qi]
    return Qi

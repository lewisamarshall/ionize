from ..constants import permittivity, boltzmann, kelvin, elementary_charge, \
                        lpm3, avagadro

def debye(self, temperature=None):
    """Return the Debye length of the solution.

    Uses the Debye-Huckel approximation for the calculation
    """
    with self.temperature(temperature):
        dielectric = self._solvent.dielectric(temperature)
        viscosity = self._solvent.viscosity(temperature)
        lamda = (dielectric * permittivity * boltzmann * kelvin(temperature) /
                 elementary_charge**2 /
                 (self.ionic_strength * lpm3) / avagadro) ** .5
        return lamda

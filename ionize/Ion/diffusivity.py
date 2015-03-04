def diffusivity(self, pH):
    """Return the diffusivity of the species at a specified pH.

    The diffusivity is returned in units of m^2/s.
    """
    diffusivity = sum([m * f / float(z) for
                      m, f, z in zip(self.absolute_mobility,
                                     self.ionization_fraction(pH),
                                     self.z
                                     )]) * self._kB * (self.T + 273.15)
    diffusivity /= sum(self.ionization_fraction(pH))
    return diffusivity

from __future__ import division


def separability(self, other, pH=None, ionic_strength=None, temperature=None):
    """Return the separability between this ion and the other ion.

    Separabillity is a normalized difference between this ion and the other
    on. Context is taken from this ion if available. 
    """
    with other.context(self.context()):
        mobilities = [ion.mobility(pH, ionic_strength, temperature)
                      for ion in (self, other)]

        # Self.mobility in denominator
        return abs((mobilities[0] - mobilities[1]) / mobilities[0])

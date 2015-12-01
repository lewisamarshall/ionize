def debye(self):
    """Return the Debye length of the solution.

    Uses the Debye-Huckel approximation for the calculation.
    """
    return self._solvent.debye(self.ionic_strength, self.temperature())

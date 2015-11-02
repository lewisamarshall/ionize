def debye(self, temperature=None):
    """Return the Debye length of the solution.

    Uses the Debye-Huckel approximation for the calculation
    """
    if temperature is None:
        temperature = self.temperature()

    with self.temperature(temperature):
        return self._solvent.debye(self.ionic_strength, self.temperature())

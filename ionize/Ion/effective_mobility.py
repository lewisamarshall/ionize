def effective_mobility(self, pH=None, ionic_strength=None, temperature=None):
    """Return the effective mobility of the ion at a given pH and I.

    Args:
        pH (float): The ambiant pH.

        I (float): The ambiant ionic strength.

    If an actual mobility from Onsager-Fouss is available, it is used,
    otherwise, the Robinson-Stokes mobility estimate is used.

    If the Ion is nested in a Solution, ok to call without a pH.

    >>> Solution(myIon, .1).ions[0].effective_mobility()

    Otherwise, always call with a pH argument.
    """
    pH, ionic_strength, temperature = \
        self._resolve_context(pH, ionic_strength, temperature)

    ionization_fraction = self.ionization_fraction(pH, ionic_strength,
                                                   temperature)
    actual_mobility = self.actual_mobility(ionic_strength, temperature)

    effective_mobility = sum(ionization_fraction * actual_mobility)

    return effective_mobility

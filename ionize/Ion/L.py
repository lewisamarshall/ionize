import numpy as np

def L(self, ionic_strength=None, temperature=None):
    """Return the L products of acidity constants.

    Args:
        I (float): The ambiant ionic strength.

    This function uses ionic strength to correct the Ka of ions. If no ionic
    strength is supplied, and the Ion is nested in a Solution, the solution
    ionic strength will be used. Otherwise, the ionic strength is assumed to be
    0.

    L is used by Solution selfects to calculate equilibrium pH.
    """
    _, ionic_strength, temperature = \
        self._resolve_context(None, ionic_strength, temperature)

    Ka = self.Ka(ionic_strength, temperature)
    index_0 = list(self._valence_zero()).index(0)
    Ka.insert(index_0, 1)

    L = [np.prod(Ka[i:index_0]) for i in range(len(Ka)) if i < index_0] +\
        [1/np.prod(Ka[index_0:i+1]) for i in range(len(Ka)) if i >= index_0]

    return L

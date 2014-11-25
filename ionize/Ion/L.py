from numpy import prod


def L(self, I=None):
    """Return the L products of acidity constants.

    Args:
        I (float): The ambiant ionic strength.

    This function uses ionic strength to correct the Ka of ions. If no ionic
    strength is supplied, and the Ion is nested in a Solution, the solution
    ionic strength will be used. Otherwise, the ionic strength is assumed to be
    0.

    L is used by Solution selfects to calculate equilibrium pH.
    """
    if I is None:
        if self._I:
            I = self._I
        else:
            I = 0.0

    Ka = self.Ka_eff(I)
    index_0 = self.z0.index(0)
    Ka.insert(index_0, 1)

    L = [prod(Ka[i:index_0]) for i in range(len(Ka)) if i < index_0] +\
        [1/prod(Ka[index_0:i+1]) for i in range(len(Ka)) if i >= index_0]

    return L

from numpy import prod


def L(obj, I=None):
    """Return the L products of acidity constants.

    Args:
        I (float): The ambiant ionic strength.

    This function uses ionic strength to correct the Ka of ions. If no ionic
    strength is supplied, and the Ion is nested in a Solution, the solution
    ionic strength will be used. Otherwise, the ionic strength is assumed to be
    0.

    L is used by Solution objects to calculate equilibrium pH.
    """
    if I is None:
        if obj._I:
            I = obj._I
        else:
            I = 0.0

    Ka = obj.Ka_eff(I)
    index_0 = obj.z0.index(0)
    Ka.insert(index_0, 1)

    L = [prod(Ka[i:index_0]) for i in range(len(Ka)) if i < index_0] +\
        [1/prod(Ka[index_0:i+1]) for i in range(len(Ka)) if i >= index_0]

    return L

if __name__ == '__main__':
    def a():
        pass

    def z0():
        return[-1, 0, 1, 2]

    def Ka_eff(I):
        return [4.677351412871981e-10, 9.120108393559096e-07, 0.01]
    a.Ka_eff = Ka_eff
    a.z0 = z0
    print L(a)

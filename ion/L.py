from numpy import cumprod


def L(obj, I=0):
    """Return the L products of acidity constants.

    These products are used in the pH calculating routine.
    It can use ionic strength correction if an ionic strength is specified.
    Otherwise, it uses uncorrected acidity coefficients.
    """
    L = obj.z0()
    Ka = obj.Ka_eff(I)
    index_0 = obj.z0().index(0)
    L[index_0] = 1

    if index_0 is not 0:
        for i in range(index_0-1, -1, -1):
            L[i] = L[i+1]*Ka[i]

    if index_0 is not len(L)-1:
        for i in range((index_0+1), len(L)):
            L[i] = L[i-1]/Ka[i-1]
    return L

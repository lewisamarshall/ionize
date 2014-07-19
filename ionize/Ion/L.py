from numpy import cumprod, prod


def L(obj, I=0):
    """Return the L products of acidity constants.

    These products are used in the pH calculating routine.
    It can use ionic strength correction if an ionic strength is specified.
    Otherwise, it uses uncorrected acidity coefficients.
    """
    Ka = obj.Ka_eff(I)
    index_0 = obj.z0().index(0)
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

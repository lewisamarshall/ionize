import numpy
from math import log10


def calc_pH(obj, I=0):
    """Return the pH of the object.

    If an ionic strength is specified, uses the corrected acidity constants.
    This function should be used only when finding the equilibrium state.
    After that, the value should be pulled from obj.pH.

    If ionic strength does not exist, assume it is zero.
    This function is used to find the equilibrium state,
    so it cannot pull the ionic strength from the object.
    """
    # Find the order of the polynomial. This is the maximum
    # size of the list of charge states in an ion.
    MaxCol = max([max(i.z)-min(i.z)+2 for i in obj.ions])

    # Set up the matrix of Ls, the multiplication
    # of acidity coefficients for each ion.
    LMat = numpy.zeros([len(obj.ions), MaxCol])

    for i in range(len(obj.ions)):
        LMat[i, 0:len(obj.ions[i].z)+1] = obj.ions[i].L(I)

    # Construct Q vector.
    Q = 1.0
    for j in range(LMat.shape[0]):
        Q = numpy.convolve(Q, LMat[j, :])

    # Convolve with water dissociation.
    Q = numpy.convolve(Q, [-obj.Kw_eff(I), 0.0, 1.0])

    # Construct P matrix
    PMat = []
    for i in range(len(obj.concentrations)):
        z_list = obj.ions[i].z0

        tmp = numpy.zeros([1, LMat.shape[1]])
        tmp[0, 0:len(z_list)] = z_list
        Mmod = LMat.copy()
        Mmod[i, :] = Mmod[i, :] * tmp

        Pi = 1
        for kl in range(Mmod.shape[0]):
            Pi = numpy.convolve(Pi, Mmod[kl, :])

        Pi = numpy.convolve([0.0, 1.0], Pi)  # Convolve with P2
        PMat.append(Pi)

    PMat = numpy.array(PMat, ndmin=2)

    # Multiply P matrix by concentrations, and sum.
    C = numpy.tile((obj.concentrations), ((PMat.shape[1], 1))).transpose()
    P = numpy.sum(PMat*C, 0)

    # Construct polynomial.
    poly = numpy.array([0.0] * max(len(P), len(Q)))
    poly[0:len(P)] = numpy.add(poly[0:len(P)], P)
    poly[0:len(Q)] = numpy.add(poly[0:len(Q)], Q)  # from QMat

    # format for the poly function.
    poly = list(poly)
    poly.reverse()

    # Solve Polynomial for concentration
    roo = numpy.roots(poly)

    roo_reduced = [r for r in roo if r.real > 0 and r.imag == 0]
    if roo_reduced:
        cH = float(roo_reduced[-1].real)
    else:
        print 'Failed to find pH.'

    # Convert to pH. Use the activity to correct the calculation.
    pH = -log10(cH*obj._H.activity_coefficient(I, [1])[0])
    return pH

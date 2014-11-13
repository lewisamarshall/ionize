import numpy
from math import log10
# pylint: disable=W0212


def calc_pH(self, I=0):
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
    max_columns = max([max(i.z)-min(i.z)+2 for i in self.ions])
    n_ions = len(self.ions)

    # Set up the matrix of Ls, the multiplication
    # of acidity coefficients for each ion.
    l_matrix = numpy.zeros([n_ions, max_columns])

    for i in range(n_ions):
        l_matrix[i, 0:len(self.ions[i].z)+1] = self.ions[i].L(I)

    # Construct Q vector.
    Q = 1.0
    for j in range(l_matrix.shape[0]):
        Q = numpy.convolve(Q, l_matrix[j, :])

    # Convolve with water dissociation.
    Q = numpy.convolve(Q, [-self.Kw_eff(I), 0.0, 1.0])

    # Construct P matrix
    PMat = []
    for i in range(n_ions):
        z_list = self.ions[i].z0

        tmp = numpy.zeros([1, l_matrix.shape[1]])
        tmp[0, 0:len(z_list)] = z_list
        Mmod = l_matrix.copy()
        Mmod[i, :] = Mmod[i, :] * tmp

        Pi = 1
        for kl in range(Mmod.shape[0]):
            Pi = numpy.convolve(Pi, Mmod[kl, :])

        Pi = numpy.convolve([0.0, 1.0], Pi)  # Convolve with P2
        PMat.append(Pi)

    PMat = numpy.array(PMat, ndmin=2)

    # Multiply P matrix by concentrations, and sum.
    C = numpy.tile((self.concentrations), ((PMat.shape[1], 1))).transpose()
    P = numpy.sum(PMat*C, 0)

    # Construct polynomial. Change the shapes as needed.
    if len(P) < len(Q):
        P.resize(Q.shape)
    elif len(P) > len(Q):
        Q.resize(P.shape)

    poly = numpy.add(P, Q)

    # Solve Polynomial for concentration
    # reverse order for poly function
    roo = numpy.roots(poly[::-1])

    roo_reduced = [r for r in roo if r.real > 0 and r.imag == 0]
    if roo_reduced:
        cH = float(roo_reduced[-1].real)
    else:
        print 'Failed to find pH.'

    # Convert to pH. Use the activity to correct the calculation.
    pH = -log10(cH * self._H.activity_coefficient(I, [1])[0])
    return pH

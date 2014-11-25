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
    l_matrix = numpy.array([numpy.resize(i.L(I),
                           [max_columns]) for i in self.ions])

    # Construct Q vector.
    Q = 1.0
    for j in range(n_ions):
        Q = numpy.convolve(Q, l_matrix[j, :])

    # Convolve with water dissociation.
    Q = numpy.convolve(Q, [-self.Kw_eff(I), 0.0, 1.0])

    # Construct P matrix
    PMat = []
    for i in range(n_ions):
        z_list = numpy.resize(self.ions[i].z0, [max_columns])

        Mmod = l_matrix.copy()
        Mmod[i, :] *= numpy.array(z_list)

        Pi = 1.
        for kl in range(n_ions):
            Pi = numpy.convolve(Pi, Mmod[kl, :])

        Pi = numpy.convolve([0.0, 1.0], Pi)  # Convolve with P2
        PMat.append(Pi)

    # Multiply P matrix by concentrations, and sum.
    P = numpy.sum(numpy.array(PMat, ndmin=2) *
                  numpy.array(self.concentrations)[:, numpy.newaxis], 0)

    # Construct polynomial. Change the shapes as needed, then reverse the order
    if len(P) < len(Q):
        P.resize(Q.shape)
    elif len(P) > len(Q):
        Q.resize(P.shape)
    poly = (P+Q)[::-1]

    # Solve Polynomial for concentration
    # reverse order for poly function
    cH = numpy.roots(poly)

    # Parse the real roots of the root finding algorithms
    cH = [c for c in cH if c.real > 0 and c.imag == 0]
    if len(cH) == 1:
        cH = float(cH[0])
    elif len(cH) > 1:
        print 'Found multiple possible pH solutions. Choosing one.'
        cH = float(cH[0])
    else:
        print 'Failed to find pH.'

    # Convert to pH. Use the activity to correct the calculation.
    pH = -log10(cH * self._H.activity_coefficient(I, [1])[0])
    return pH

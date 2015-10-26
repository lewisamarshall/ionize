from scipy.optimize import newton, brentq
import numbers

from ..Ion import Ion
from ..Database import Database

database = Database()


def buffering_capacity(self):
    """Return the buffering capacity of the solution.

    This function generates an approximate solution to the buffering
    capacity by finding the derivative of the pH with respect to
    the addition of an acid insult at small concentration.
    """
    # Remove any ions at concentration 0.
    c = 0.001*min([cp for cp in self.concentrations if cp > 0])
    Cb = 0

    # Add an acid insult at 0.1% the lowest concentration in the solution.
    # If the buffering capacity is measured as above the insult c,
    # make the insult c lower.
    while Cb < c:
        new_sol = self + (Ion('Acid Insult', [-1], [-2], [-1]), c)
        # Find the slope of the pH.
        Cb = abs(c/(self.pH-new_sol.pH))
        c = 0.01 * Cb
    return Cb


def titrate(self, titrant, target, titration_property='pH', return_c=False):
    """Return a Solution titrated to the target pH using the titrant.

    Titrate uses root finding to determine the concentration of titrant
    to add to achieve the desired pH. If the desired pH cannot be achieved
    by addition of the titrant, titrate raises an error.

    To titrate to a target property other than pH, simply set the property
    to a property of the Solution class.
    """
    if isinstance(titrant, basestring):
        titrant = database.load(titrant)

    att = getattr(self, titration_property)
    if isinstance(att, numbers.Number):
        min_func = lambda c:\
            getattr(self + (titrant, c), titration_property)-target
    else:
        min_func = lambda c: \
            getattr(self + (titrant, c), titration_property)()-target

    c, r = brentq(min_func, 0, 55, full_output=True)

    if r.converged:
        if return_c:
            return c
        else:
            return (self + (titrant, c))
    else:
        raise RuntimeError('Solver did not converge.')

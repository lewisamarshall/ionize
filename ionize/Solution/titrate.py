from scipy.optimize import newton, brentq
from ..Ion import Ion


def titrate(self, titrant, target, property='pH', return_c=False):
    """Return a Solution titrated to the target pH using the titrant.

    Titrate uses root finding to determine the concentration of titrant
    to add to achieve the desired pH. If the desired pH cannot be achieved
    by addition of the titrant, titrate raises an error.

    To titrate to a target property other than pH, simply set the property
    to a property of the Solution class.
    """

    min_func = lambda c: (self + (titrant, c)).__dict__[property]-target

    c, r = brentq(min_func, 0, 55, full_output=True)
    if r.converged:
        if return_c:
            return c
        else:
            return (self + (titrant, c))
    else:
        raise

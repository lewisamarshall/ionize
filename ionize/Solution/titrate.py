from scipy.optimize import newton, brentq
from ..Ion import Ion


def titrate(self, titrant, target, property='pH', return_c=False):
    """Return a Solution titrated to the target pH using the titrant.

    Titrate uses root finding to determine the concentration of titrant
    to add to achieve the desired pH. If the desired pH cannot be achieved
    by addition of the titrant, titrate raises an error.

    Property can be sent to titrate to a chosen conductivity or ionic strength
    instead of pH, by setting property to "I" or "Conductivity".

    Normally, titrate returns the titrated solution. However, to instead return
    the titrant concentration, set return_c="True". 
    """
    if property == 'pH':
        min_func = lambda c: (self + (titrant, c)).pH - target
    elif property == 'conductivity':
        min_func = lambda c: (self + (titrant, c)).conductivity - target
    elif property == 'I':
        min_func = lambda c: (self + (titrant, c)).I - target
    else:
        raise NotImplementedError

    c, r = brentq(min_func, 0, 55, full_output=True)
    if r.converged:
        if return_c:
            return c
        else:
            return (self + (titrant, c))
    else:
        raise

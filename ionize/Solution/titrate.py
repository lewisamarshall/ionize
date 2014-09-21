from scipy.optimize import newton, brentq
from ..Ion import Ion


def titrate(self, titrant, target, property='pH'):
    if property == 'pH':
        min_func = lambda c: (self + (titrant, c)).pH - target
    else:
        raise NotImplementedError
    c, r = brentq(min_func, 0, 55, full_output=True)
    if r.converged:
        return c
    else:
        raise

from __future__ import division, print_function
from scipy.optimize import newton, brentq, root
import numpy as np
import numbers
import warnings
from copy import deepcopy

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
    if isinstance(titrant, str):
        titrant = database.load(titrant)

    att = getattr(self, titration_property)
    if isinstance(att, numbers.Number):
        def min_func(c):
            getattr(self + (titrant, c), titration_property)-target
    else:
        def min_func(c):
            getattr(self + (titrant, c), titration_property)()-target

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, r = brentq(min_func, 0, 55, full_output=True)

    if r.converged:
        if return_c:
            return c
        else:
            return (self + (titrant, c))
    else:
        raise RuntimeError('Solver did not converge.')


def equilibrate_CO2(self):
    CO2 = database['carbonic acid']
    CO2.context(self)
    # TODO: move to constants
    atmospheric_CO2 = 0.0004
    eq = atmospheric_CO2 * self._solvent.henry_CO2(self.temperature())

    def min_func(concentration):
        self._contents[CO2] = concentration
        self._equilibrate()
        deionized = concentration * (1-sum(CO2.ionization_fraction()))
        return deionized - eq

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, r = brentq(min_func, 0, 1., full_output=True)

    if r.converged:
        self._contents[CO2] = c
    else:
        raise RuntimeError('Solver did not converge')


def displace(self, receding, advancing):
    """Electrophoretically displace an ion."""
    # Convert ion names to ions
    if isinstance(advancing, str):
        advancing = database[advancing]
    if isinstance(receding, str):
        receding = database[receding]

    # Check that displacement is possible.
    assert advancing not in self, \
        'The advancing ion is already present.'
    assert receding in self, \
        'The receding ion is not present.'

    # Make a new solution to perform optimization.
    new_solution = deepcopy(self)
    new_solution._contents[advancing] = \
        new_solution._contents.pop(receding)

    # Give the advancing and receding ions the right context.
    advancing.context(new_solution), receding.context(self)
    new_solution._equilibrate()

    def min_func(concentrations):
        for ion, concentration in zip(new_solution.ions,
                                      concentrations):
            new_solution._contents[ion] = concentration
        new_solution._equilibrate()

        field_ratio = receding.mobility()/advancing.mobility()

        delta = [new_solution.concentration(ion) -
                 self.concentration(ion) / field_ratio *
                 (self[ion].mobility() - receding.mobility()) /
                 (new_solution[ion].mobility() - advancing.mobility())
                 for ion in new_solution.ions if ion is not advancing]

        delta.append(field_ratio -
                     self.conductivity() / new_solution.conductivity())

        return np.array(delta)

    sol = root(min_func, new_solution.concentrations)

    if sol['success']:
        return new_solution
    else:
        msg = 'Solver failed to converge. {}'.format(sol['message'])
        raise RuntimeError(msg)

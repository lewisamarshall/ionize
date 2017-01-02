from __future__ import division, print_function
from scipy.optimize import brentq, root
import numpy as np
import numbers
import warnings
from copy import deepcopy

from ..Ion import Ion
from ..Database import Database
from ..constants import atmospheric_CO2


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


def titrate(self, titrant, target, titration_property='pH'):
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
            return getattr(self + (titrant, c), titration_property)-target
    else:
        def min_func(c):
            return getattr(self + (titrant, c), titration_property)()-target

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, r = brentq(min_func, 0, 55, full_output=True)

    if r.converged:
        return (self + (titrant, c))
    else:
        raise RuntimeError('Solver did not converge.')


def equilibrate_CO2(self, partial_pressure=atmospheric_CO2):
    """Titrate the CO2 in solution to equilibrium with the atmosphere.

    :param partial_pressure: The partial pressure of CO2 in the atmosphere,
    in bar. Defaults to the typical value of Earth's atmosphere.
    """
    CO2 = database['carbonic acid']
    eq = partial_pressure * self._solvent.henry_CO2(self.temperature())

    def min_func(concentration):
        new_sol = self + (CO2, concentration)
        with CO2.context(new_sol):
            deionized = concentration * (1-sum(CO2.ionization_fraction()))
        return deionized - eq

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, r = brentq(min_func, 0, 1., full_output=True)

    if r.converged:
        return self + (CO2, c)
    else:
        raise RuntimeError('Solver did not converge')


def displace(self, receding, advancing=None, guess=None):
    """Electrophoretically displace an ion."""
    # Convert ion names to ions
    if isinstance(advancing, str):
        advancing = database[advancing]
    if isinstance(receding, str):
        receding = database[receding]

    receding.context(self)

    # Check that displacement is possible.
    assert advancing not in self, \
        'The advancing ion is already present.'
    assert receding in self, \
        'The receding ion is not present.'

    # Make a new solution to perform optimization.
    new_solution = deepcopy(self)
    if advancing is not None:
        new_solution._contents[advancing] = \
            new_solution._contents.pop(receding)
        advancing.context(new_solution)
    else:
        new_solution._contents.pop(receding)

    # If there is a guess, use it to update the concentrations
    if guess is not None:
        assert len(new_solution.ions) == len(guess), 'Incorrect guess length.'
        for ion, c_guess in zip(new_solution._contents.keys(), guess):
            new_solution._contents[ion] = c_guess

    new_solution._equilibrate()

    def min_func(concentrations):
        for ion, concentration in zip(new_solution.ions,
                                      concentrations):
            new_solution._contents[ion] = abs(concentration)

        new_solution._equilibrate()

        velocity = 1./self.zone_transfer(receding)

        err = []

        for ion in new_solution.ions:
            if ion is advancing:
                err.append(1./new_solution.zone_transfer(ion) - velocity)
            else:
                err.append(new_solution.concentration(ion) -
                           self.concentration(ion) *
                           (1./self.zone_transfer(ion)-velocity) /
                           (1./new_solution.zone_transfer(ion)-velocity)
                           )

        return np.array(err)

    try:
        sol = root(min_func, new_solution.concentrations, method='lm')
    except Exception as e:
        msg = 'Solver failed on iteration. \nSolution: {}\nError: {}'
        raise RuntimeError(msg.format(repr(new_solution), e))

    if sol['success']:
        return new_solution
    else:
        msg = 'Solver failed to converge. \nSolution: {}\nError: {}'
        raise RuntimeError(msg.format(repr(new_solution), sol['message']))

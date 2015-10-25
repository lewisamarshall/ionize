import numpy as np
from ..constants import lpm3


def transference(self):
    """Return the fraction of charge carried by each of the ions as a list.

    Should not precisely add to 1, because some charge is carried by protons
    and hydroxyls.
    """
    T = self.concentrations * _molar_conductivities(self) / self.conductivity()
    return T


def zone_transfer(self):
    """Return the zone transfer charge of the solution per liter."""

    return self.conductivity() / _mobilities(self) / lpm3


def _molar_conductivities(solution):
    return np.array([ion.molar_conductivity() for ion in solution.ions])


def _mobilities(solution):
    return np.array([ion.mobility() for ion in solution.ions])

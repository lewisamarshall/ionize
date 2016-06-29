from ..constants import lpm3
from ..Database import Database

database = Database()


def transference(self, ion):
    """Return the fraction of charge carried by each of the ions as a list.

    Should not precisely add to 1, because some charge is carried by
    dissociated water.
    """
    if ion in self:
        if isinstance(ion, str):
            ion = self._name_lookup[ion]
        with ion.context(self):
            return (self.concentration(ion) * ion.molar_conductivity() /
                    self.conductivity())
    else:
        return 0


def zone_transfer(self, ion):
    """Return the zone transfer ion in the solution, in Coulombs/liter."""
    if isinstance(ion, str):
        ion = database[ion]

    with ion.context(self):
        return self.conductivity() / ion.mobility() / lpm3

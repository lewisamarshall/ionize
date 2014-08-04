"""Represent ions in aqueous solution.

Class Ion represents ion species.

Class Solution represents an aqueous solution containing one or more ions.

Function load_ion loads these ions from a database housed in ions_shelve.db.

Function search_ion searches the database.

Function get_db returns the database as a dictionary.
"""

from .Ion import Ion
from .Solution import Solution
from .get_db import get_db
from .load_ion import load_ion
from .search_ion import search_ion
from .viscosity import viscosity as _viscosity
from .dielectric import dielectric as _dielectric


def viscosity(T):
    """Return the viscosity of water at temperature T.

    T should be entered in Celcius.
    """
    return _viscosity(None, T)


def dielectric(T):
    """Return the dielectric constant of water at temperature T.

    T should be entered in Celcius.
    """
    return _dielectric(None, T)

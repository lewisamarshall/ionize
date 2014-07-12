"""Represent ions in aqueous solution.

Class Ion represents ion species.

Function load_ion loads these ions from a database housed in ions_shelve.db.
The ions in this database are largely based on the Hirokawa database.

Class Solution represents an aqueous solution containing one or more ions.
"""

from Ion import Ion
from Solution import Solution
from load_ion import load_ion

if __name__ == "__main__":
    hcl = load_ion('hydrochloric acid')
    print hcl
    tris = load_ion('histidine')
    print tris
    buf = Solution([hcl, tris], [0.03, 0.06])
    print buf.I
    print buf.pH
    print buf.conductivity()
    print buf.ions[1].actual_mobility
    print buf.ions[1].absolute_mobility

"""Represent ions in aqueous solution.

Class Ion represents ion species.

Function load_ion loads these ions from a database housed in ions_shelve.db.
The ions in this database are largely based on the Hirokawa database.

Class Solution represents an aqueous solution containing one or more ions.
"""

from Ion import Ion
from Solution import Solution
from load_ion import load_ion
from search_ion import search_ion

if __name__ == "__main__":
    # water = Solution()
    hcl = load_ion('hydrochloric acid')
    print hcl
    tris = load_ion('tris')
    print tris
    buf = Solution([hcl, tris], [0.03, 0.06])
    print "I", buf.I
    print "pH", buf.pH
    print "conductivity", buf.conductivity()
    print "tris actual mobility", buf.ions[1].actual_mobility
    print "tris absolute mobility", buf.ions[1].absolute_mobility
    print "buffering capacity", buf.buffering_capacity()
    print "Kw_eff", buf.Kw_eff()
    print "new solution pH", buf.add_ion([load_ion('histidine')], [0.01]).pH
    print "transference", buf.transference()
    print "zone transfer", buf.zone_transfer(0.001)
    print buf.__dict__
    help(Solution)

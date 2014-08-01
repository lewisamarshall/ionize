import warnings
from numpy import mean


def alberty(self):
    """Return the Alberty conservation function value of a solution.

    Throw a warning if the ions are not univalent, or if water dissociation
    cannot be neglected.
    """
    al = 0

    for ion, c in zip(self.ions, self.concentrations):
        if len(ion.actual_mobility) == 1:
            al += c * ion._Lpm3 / abs(ion.actual_mobility[0])
        else:
            f = ion.ionization_fraction()
            if max(f)/sum(f) < .9:
                warnings.warn('Ion not in single valance. Alberty invalid.')
            elif abs(ion.z[f.index(max[f])]) != 1:
                warnings.warn('Ion valance is not 1. Alberty invalid.')
            al += c * ion._Lpm3 / abs(ion.actual_mobility[f.index(max(f))])

    return al

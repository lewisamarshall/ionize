import warnings
from numpy import mean


def kohlrausch(obj):
    """Return the Kohlrausch regulating function (KRF) value of a solution.

    Throw a warning if the ions in the solution are not near fully ionized.
    """
    KRF = 0

    for ion, c in zip(obj.ions, obj.concentrations):
        z_eff = (mean([z*f for z, f in
                 zip(ion.z, ion.ionization_fraction())]))
        KRF += abs(z_eff) * c * ion._Lpm3 / ion.effective_mobility()
        if max(ion.ionization_fraction()) < .9:
            warnings.warn('ions are not fully ionized. KRF is a poor approx.')

    return KRF


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


def jovin(self):
    """Return the Jovin conservation function value of a solution.

    Throw a warning if the ions are not univalent, or if water dissociation
    cannot be neglected.
    """
    jov = 0

    for ion, c in zip(self.ions, self.concentrations):
        if len(ion.actual_mobility) == 1:
            jov += c * ion.z[0]
        else:
            f = ion.ionization_fraction()
            if max(f)/sum(f) < .9:
                warnings.warn('Ion not in single valance. Jovin invalid.')
            elif abs(ion.z[f.index(max[f])]) != 1:
                warnings.warn('Ion valance is not 1. Jovin invalid.')
            jov += c * ion.z[f.index(max(f))]

    return jov


def gas(self):
    pass

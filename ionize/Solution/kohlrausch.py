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

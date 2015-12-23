from __future__ import division
import warnings
import numpy as np
from ..constants import lpm3


def kohlrausch(self):
    """Return the Kohlrausch regulating function (KRF) value of a solution.

    The Kohlrausch regulating function is only valid if ions are near full
    ionization. This function will deliver a warning where this is not the
    case.
    """
    KRF = 0

    for ion in self.ions:

        z_eff = (np.mean([z*f for z, f in
                 zip(ion.valence, ion.ionization_fraction())]))
        KRF += abs(z_eff) * self.concentration(ion) * lpm3 / ion.mobility()
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
        if len(ion.reference_mobility) == 1:
            al += c * lpm3 / abs(ion.actual_mobility()[0])
        else:
            f = ion.ionization_fraction()
            if max(f)/sum(f) < .9:
                warnings.warn('Ion not in single valance. Alberty invalid.')
            elif abs(ion.valence[f.argmax()]) != 1:
                warnings.warn('Ion valance is not 1. Alberty invalid.')
            al += c * lpm3 / abs(ion.actual_mobility()[f.argmax()])

    return al


def jovin(self):
    """Return the Jovin conservation function value of a solution.

    Throw a warning if the ions are not univalent, or if water dissociation
    cannot be neglected.
    """
    jov = 0

    for ion, c in zip(self.ions, self.concentrations):
        if len(ion.reference_mobility) == 1:
            jov += c * ion.valence[0]
        else:
            f = ion.ionization_fraction()
            if max(f)/sum(f) < .9:
                warnings.warn('Ion not in single valance. Jovin invalid.')
            # elif abs(ion.valence[f.index(max(f))]) != 1:
            elif abs(ion.valence[f.argmax()]) != 1:
                warnings.warn('Ion valance is not 1. Jovin invalid.')
            jov += c * ion.valence[f.argmax()]

    return jov


def gas(self):
    """Return the Gas conservation function value of the solution.

    There are two Gas conservation function values, one that depends
    on the mobility of hydronium, the other depending on the mobility
    of hydroxide. Each of these is only valid when the concentration
    of the other ion is negligable. When this assumption is broken,
    this function returns NaN for one value and a float for the other.
    """
    alberty = self.alberty()
    jovin = self.jovin()
    gas = [alberty - jovin / self._hydronium.mobility(self.pH),
           alberty - jovin / self._hydroxide.mobility(self.pH)]
    if not self.safe():
        if self.concentration('H+') > self.concentration('OH-'):
            gas[1] = float('NaN')
        else:
            gas[0] = float('NaN')
    return np.array(gas)

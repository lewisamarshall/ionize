from __future__ import division

from .PolyIon import PolyIon
from ..Ion import fixed_state
from ..constants import boltzmann, kelvin, reference_temperature, \
    elementary_charge, avogadro, lpm3, gpkg

from math import pi, exp
import numpy as np

from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint, positive_pKs, \
    negative_pKs, pKcterminal, pKnterminal


@fixed_state
class Peptide(PolyIon):
    """Peptide represents single protein chains in solution.

    Peptides properties are based entirely on analysis of the sequence of the
    peptide.
    """

    _state = {'name': 'Name of the peptide.',
              'sequence': 'Amino acid sequence of the peptide.'
              }

    _sequence = None
    _analysis = None

    # TODO: move h to function or constants. Unify with pitts?
    _h_max = 1
    _h_min = 2./3.
    _h = 5./6.

    def __init__(self, name=None, sequence=None):
        self._name = name
        self._sequence = sequence
        self._analysis = ProteinAnalysis(str(self.sequence))

    @property
    def molecular_weight(self):
        return SeqUtils.molecular_weight(self.sequence, 'protein')

    def charge(self, pH=None, ionic_strength=None, temperature=None,
               moment=1):
        """Return the time-averaged charge of the peptide.

        :param pH
        :param ionic_strength
        :param temperature
        """
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        amino_acid_count = self._analysis.count_amino_acids()

        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)

        nterm = self.sequence[0]
        cterm = self.sequence[-1]

        if nterm in pKnterminal:
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs['Cterm'] = pKcterminal[cterm]

        charge = IsoelectricPoint(self.sequence,
                                  amino_acid_count)._chargeR(pH,
                                                             pos_pKs,
                                                             neg_pKs)
        return charge**moment

    def isoelectric_point(self, ionic_strength=None, temperature=None):
        """Return the isoelectric point of the peptide."""
        # _, ionic_strength, temperature = \
        #     self._resolve_context(None, ionic_strength, temperature)
        return self._analysis.isoelectric_point()

    def volume(self):
        """Return the approximate volume of the folded peptide in m^3."""
        v = self.molecular_weight / avogadro / self.density() / lpm3 / gpkg
        return v

    def radius(self):
        """Return the approximate radius of the folded peptide in m."""
        return (self.volume() * 3. / 4. / pi) ** (1. / 3.)

    def density(self):
        """Return the approximate density of the folded peptide in kg/L."""
        return 1.410 + 0.145 * exp(-self.molecular_weight / 13.)

    def mobility(self, pH=None, ionic_strength=None, temperature=None):
        """Return the effective mobility of the ion in m^2/V/s.

        If a context solution is available, mobility uses the full Onsager-Fuoss
        correction to mobility. Otherwise, the Robinson-Stokes model is used.

        :param pH
        :param ionic_strength
        :param temperature
        """
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        mobility = self.charge(pH) * elementary_charge /\
            (6 * pi * self._solvent.viscosity(temperature) * self.radius() *
             (1 + self.radius() /
              self._solvent.debye(ionic_strength, temperature)
              )
             ) * self._h
        return mobility

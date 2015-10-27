from .PolyIon import PolyIon
from ..constants import boltzmann, kelvin, reference_temperature, \
                        elementary_charge, avagadro
from ..Ion import fixed_state

from math import pi, exp
import numpy as np

from Bio import SeqUtils
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint, \
                                          positive_pKs, negative_pKs, \
                                          pKcterminal, pKnterminal
from Bio.SeqUtils.ProtParam import ProteinAnalysis


@fixed_state
class Peptide(PolyIon):

    _state = ('name',
              'sequence')

    _sequence = None
    _analysis = None

    # TODO: move h to function or constants
    _h_max = 1
    _h_min = 2./3.
    _h = 5./6.

    def __init__(self, name=None, sequence=None):
        self._name = name
        self._sequence = sequence
        self._analysis = ProteinAnalysis(str(self.sequence))

    def molecular_weight(self):
        return SeqUtils.molecular_weight(self.sequence, 'protein')

    def charge(self, pH=None, ionic_strength=None, temperature=None):
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

        return IsoelectricPoint(self.sequence,
                                amino_acid_count)._chargeR(pH,
                                                           pos_pKs,
                                                           neg_pKs)

    def isoelectric_point(self, ionic_strength=None, temperature=None):
        # _, ionic_strength, temperature = \
        #     self._resolve_context(None, ionic_strength, temperature)
        return self._analysis.isoelectric_point()

    def volume(self):
        v = self.molecular_weight() / avagadro / self.density() / (100**3)
        return v

    def radius(self):
        return (self.volume() * 3. / 4. / pi) ** (1. / 3.)

    def density(self):
        return 1.410 + 0.145 * exp(-self.molecular_weight() / 13.)

    def mobility(self, pH=None, ionic_strength=None, temperature=None):
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        mobility = self.charge(pH) * elementary_charge /\
            (6 * pi * self._solvent.viscosity(temperature) * self.radius() *
             (1 + self.radius() /
              self._solvent.debye(ionic_strength, temperature)
              )
             ) * self._h
        return mobility

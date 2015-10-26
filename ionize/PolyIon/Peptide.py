from .PolyIon import PolyIon
from ..constants import boltzmann, kelvin, reference_temperature, \
                        elementary_charge, avagadro
from ..Ion import fixed_state

from Bio import PDB, SeqUtils
from math import pi, exp
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint, \
                                          positive_pKs, negative_pKs, \
                                          pKcterminal, pKnterminal
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import tempfile
import numpy as np

lister = PDB.PDBList()
parser = PDB.PDBParser()
builder = PDB.PPBuilder()


@fixed_state
class Peptide(PolyIon):

    _state = ('name',
              'sequence')

    _sequence = None
    _analysis = None

    h_max = 1
    h_min = 2./3.
    h = 5./6.


    def __init__(self, name=None, sequence=None):
        self.name = name
        if sequence is None:
            self.get_file()
            self.set_structure()
            self.set_sequence()
        else:
            self._sequence = sequence

        # TODO: Fix that this only looks at first sequence
        self._analysis = ProteinAnalysis(str(self.sequence[0]))

    def _get_sequence(self):
        self.file = lister.retrieve_pdb_file(self.name)
        self.structure = parser.get_structure(self.name, self.file)
        self._sequence = [pp.get_sequence()
                          for pp in
                          builder.build_peptides(self.structure)]

    def molecular_weight(self):
        return np.array([SeqUtils.molecular_weight(sequence, 'protein')
                         for sequence in self.sequence])

    def charge(self, pH=None, ionic_strength=None, temperature=None):
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        amino_acid_count = self._analysis.count_amino_acids

        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)

        nterm = self.sequence[0]
        cterm = self.sequence[-1]

        if nterm in pKnterminal:
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs['Cterm'] = pKcterminal[cterm]
        return IsoelectricPoint(self.sequence,
                                self.amino_acid_count)._chargeR(pH,
                                                                pos_pKs,
                                                                neg_pKs)

    def isoelectric_point(self, ionic_strength=None, temperature=None):
        pH, ionic_strength, temperature = \
            self._resolve_context(pH, ionic_strength, temperature)

        return self.analysis.isoelectric_point()

    def volume(self):
        return self.mw[0]/self.density()

    def radius(self):
        self.radius = self.volume() * 3. / 4. / pi ** (1. / 3.) / 100.

    def density(self):
        return 1.410 + 0.145 * exp(-self.molecular_weight[0] / 13.)

    def mobility(self, pH, ionic_strength, temperature):
        mobility = self.charge(pH) /\
            (6 * pi * self._solvent.viscosity(temperature) * self.radius() *
             (1 + self.kappa * self.radius()))*self.h
        return mobility

    # def set_debye(self):
    #     I = 0.001 * 1000
    #     self.lamda = (self.epsilon*self.k*self.T/self.e**2/I/self.Na)**.5
    #     self.kappa = 1./self.lamda

from .IonComplex import IonComplex
from ..PolyIon import Peptide
from ..Ion import fixed_state

from Bio import PDB
import tempfile
lister = PDB.PDBList()
parser = PDB.PDBParser()
builder = PDB.PPBuilder()


@fixed_state
class Protein(IonComplex):

    _state = ('name', 'members', 'sequences', )

    sequences = tuple()

    def __init__(self, name=None, sequences=None):
        self._name = name

        if sequences is None:
            self._get_sequences()
        else:
            self._sequences = sequences

        self._members = tuple([Peptide(name='{}:{}'.format(self.name, idx),
                                       sequence=sequence)
                               for idx, sequence
                               in enumerate(self._sequences)])

    def _get_sequences(self):
        temploc = tempfile.mkdtemp()
        file_ = lister.retrieve_pdb_file(self.name, pdir=temploc)
        structure = parser.get_structure(self.name, file_)
        sequences = [str(peptide.get_sequence()) for
                     peptide in builder.build_peptides(structure)]
        self._sequences = tuple(sequences)

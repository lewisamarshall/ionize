from .IonComplex import IonComplex
from ..PolyIon import Peptide
from ..Ion import fixed_state

import tempfile
from string import ascii_uppercase
from Bio import PDB
lister = PDB.PDBList(obsolete_pdb='override')
parser = PDB.PDBParser()
builder = PDB.PPBuilder()


@fixed_state
class Protein(IonComplex):

    _state = ('name', 'members')

    sequences = tuple()

    def __init__(self, name=None, ids=None, sequences=None, members=None):
        self._name = name

        if members is not None:
            self._members = tuple(members)
            return

        if sequences is None:
            ids, sequences = self._from_pdb()
        elif ids is None:
            ids = tuple(['{}:{}'.format(self.name, ascii_uppercase[idx])
                         for idx in range(len(sequences))])

        self._members = tuple([Peptide(name=id,
                                       sequence=sequence)
                               for id, sequence
                               in zip(ids, sequences)])

    def _from_pdb(self):
        temploc = tempfile.mkdtemp()
        try:
            file_ = lister.retrieve_pdb_file(self.name, pdir=temploc)
        except:
            raise RuntimeError(
                'Could not download {} from the PDB.'.format(self.name))
        structure = parser.get_structure(self.name, file_)

        ids = []
        sequences = []
        for chain in structure.get_chains():
            ids.append('{}:{}'.format(self.name, chain.id))
            sequences.append(
                str(builder.build_peptides(chain)[0].get_sequence()))

        return tuple(ids), tuple(sequences)

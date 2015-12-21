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
    """Protein represents an ion composed of a complex of peptides.

    :param name: Name of the protein.
    :param ids: Names of the peptide members.
    :param sequences: Sequences of the peptide members.
    :param members: An iterable of the peptide members.

    If members and sequences are not provided, the name will be searched in the
    Protein DataBase (PDB). If a protein of the same name is available, the
    sequences of the peptides will be gathered from the PDB.
    """

    _state = {'name': 'Protein name.',
             'members': 'Name of the peptide members.'
             }

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

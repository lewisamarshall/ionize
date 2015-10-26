from .PolyIon import PolyIon
from ..Ion import fixed_state

@fixed_state
class NucleicAcid(PolyIon):

    _state = ('name',
              'size',
              'species',
              'sequence')

    _name = 'Nucleic Acid'
    _size = float('inf')
    _species = 'DNA'
    _sequence = None

    __species_options = ('DNA', 'RNA', 'dsDNA', 'ssDNA', 'ssRNA', '')


    def __init__(self, name=None,
                 size=None, sequence=None,
                 species='DNA'):

        self._name = name or self.name
        self._size = size or self.size
        self._sequence = sequence or self.sequence

        assert species in self.__species_options
        self._species = species or self.species

    def mobility(self):
        mu = (3.75 - 1.8 * (self.size**-.6)) * 1e-8
        return mu

    def charge(self):
        raise NotImplementedError

    def diffusivity(self):
        raise NotImplementedError

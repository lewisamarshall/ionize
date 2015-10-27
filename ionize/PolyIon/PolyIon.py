from ..Ion import BaseIon


class PolyIon(BaseIon):

    def molecular_weight(self):
        raise NotImplementedError

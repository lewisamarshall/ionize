from .Aqueous import Aqueous


class BaseIon(object):
    _solvent = Aqueous()
    name = 'BaseIon'
    charge = None

    def __repr__(self):
        """Return a representation of the ion."""
        return "Ion('{}', z={})".format(self.name, self.charge)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def mobility(self):
        raise NotImplementedError

    def diffusivity(self):
        raise NotImplementedError

    def ionization_fraction(self):
        raise NotImplementedError

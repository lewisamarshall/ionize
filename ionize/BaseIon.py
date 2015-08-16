from .Aqueous import Aqueous


class BaseIon(object):
    _solvent = Aqueous()
    name = 'BaseIon'
    charge = None

    def __repr__(self):
        """Return a representation of the ion."""
        return "Ion('{}', z={})".format(self.name, self.charge)

    def __eq__(self, other):
        selfstate, otherstate = [{k: v for k, v in obj.__dict__.items()
                                  if not k.startswith('_')}
                                 for obj in self, other]

        return selfstate == otherstate

    def mobility(self):
        raise NotImplementedError

    def diffusivity(self):
        raise NotImplementedError

    def ionization_fraction(self):
        raise NotImplementedError

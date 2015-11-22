from ..Ion import BaseIon, fixed_state


@fixed_state
class IonComplex(BaseIon):

    _state = ('name', 'members')

    _members = tuple()

    def __init__(self, name, members):
        assert all(isinstance(content, Ion) for content in contents), \
            'All inputs to a complex must be Ions.'

        self._members = tuple([member for member in members])

    def context(self, context):
        [member.context(context) for member in self.members]
        self._context = context

    def charge(self, pH, ionic_strength, temperature):
        return sum([member.charge(pH, ionic_strength, temperature) for
                    member in self.members])

    def mobility(self, pH, ionic_strength, temperature):
        raise NotImplementedError

    def diffusivity(self, pH, ionic_strength, temperature):
        raise NotImplementedError

    def molar_conductivity(self, pH, ionic_strength, temperature):
        raise NotImplementedError

    def molecular_weight(self):
        return sum([member.molecular_weight() for member in self.members])

    def __getitem__(self, idx):
        return self.members[idx]

    def __iter__(self):
        return (member for member in self.members)

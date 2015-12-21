from ..Ion import BaseIon, fixed_state

import contextlib


@fixed_state
class IonComplex(BaseIon):
    """IonComplex is a collection of ions that moves together.

    :param name: The name for the complex.
    :param members: An iterable of ion objects that make up the complex.
    """

    _state = {'name': 'Name of the ion.',
              'members': 'Components of the ion complex.'}

    _members = tuple()

    def __init__(self, name, members):
        assert all(isinstance(member, BaseIon) for member in members), \
            'All inputs to a complex must be Ions.'

        self._name = name
        self._members = tuple([member for member in members])

    def context(self, context=False):
        """Control the context of the ion.

        Complex contexts works by setting the context of the members of the
        complex. The Complex has on context of its own.

        See the context method for the members for additional information.
        """
        old_context = [member.context() for member in self.members]

        if context is False:
            return tuple(old_context)

        # Update to new context, save old_context
        [member.context(context) for member in self.members]

        # Return a context manager to revert to old_context
        @contextlib.contextmanager
        def manager():
            yield
            [member.context(context)
             for member, old in zip(self.members, old_context)]
        return manager()

    def charge(self, pH=None, ionic_strength=None, temperature=None,
               moment=1):
        """Return the net charge for the complex."""
        return sum([member.charge(pH, ionic_strength, temperature) for
                    member in self.members])**moment

    def mobility(self, pH=None, ionic_strength=None, temperature=None):
        """Return the mobility of the complex, in m^2/V/s."""
        # NOTE: relies on the molecular-weight averaged mobility
        return sum([member.mobility(pH, ionic_strength, temperature) *
                    member.molecular_weight / self.molecular_weight
                    for member in self])

    @property
    def molecular_weight(self):
        try:
            return sum([member.molecular_weight for member in self.members])
        except TypeError:
            raise TypeError('Each member must have a numeric molecular '
                            'weight to compute the complex weight.')

    def __getitem__(self, idx):
        return self.members[idx]

    def __iter__(self):
        return (member for member in self.members)

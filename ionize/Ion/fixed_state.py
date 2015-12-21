import operator


def fixed_state(cls):
    for state_property in cls._state:
        setattr(cls,
                state_property,
                property(operator.attrgetter("_{}".format(state_property)),
                         doc=cls._state[state_property])
                )
    return cls

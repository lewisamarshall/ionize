from ..Ion import Ion


def buffering_capacity(obj):
    """Return the buffering capacity of the solution.

    This function generates an approximate solution to the buffering
    capacity by finding the derivative of the pH with respect to
    the addition of an acid insult at small concentration.
    """
    # Remove any ions at concentration 0.
    c = [cp for cp in obj.concentrations if cp > 0]
    # Choose a concentration 1% the lowest ion concentration.
    c = 0.01*min(c)

    # Add an acid insult at 1% the lowest concentration in the solution.
    new_sol = obj.add_ion([Ion('Acid Insult', -1, -2, -1)], [c])

    # Find the slope of the pH.
    Cb = abs(c/(obj.pH-new_sol.pH))
    return Cb

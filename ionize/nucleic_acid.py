import warnings
if not __name__ == '__main__':
    from Ion import Ion


def nucleic_acid(bp=float('inf')):
    """Return an ion representing nucleic acid.

    Takes bp as an argument, representing the length of the nucleic acid in
    base pairs. If no bp argument is given, assumes that the nucleic acid is
    long.

    This function returns the free solution mobility of nucleic acid. If this
    mobility is modified by a seiving matrix, it will differ significantly
    from these estimates. The correlation between size and mobility is drawn
    from Stellwagen 1997, 'The Free Solution Mobility of DNA'.
    """
    mu = (3.75 - 1.8*(bp**-.6)) * 1e-8
    return Ion('DNA: {} bp'.format(bp), -1, -1, -mu)

if __name__ == '__main__':
    DNA()
    DNA(1000)

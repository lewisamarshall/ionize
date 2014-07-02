import warnings


def load_ion(ion_name):
    """Return an ion by name from the database.
    ion=load_ion('ion_name') pulls the named ion from the database.
    Database derived from Peakmaster."""

    load database.mat
    if ion_name in NAMES:
        N = find(strcmpi(ion_name, NAME))
        indices =~ isnan(PKA(N, :))
        loaded_ion = ion(NAME{N},
                         Z(indices),
                         PKA(N, indices),
                         MOBILITY(N, indices).*1e-9.*sign(Z(indices)))
        return loaded_ion
    else:
        warnings.warn('Ion not found in database. Returning None.')
        return None

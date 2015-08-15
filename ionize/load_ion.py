import warnings
from .Ion import Ion
from get_db import get_db


def load_ion(ion_name, solvation=True):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.

    By default load solvation data for small ions from the Nightingale
    files. Pass in solvation=False to disable.
    """
    ion_name = ion_name.lower()  # Convert to lower case.
    ion_list = get_db()  # Open the database.
    if ion_name in ion_list.keys():
        ion_entry = ion_list[ion_name]
        ion_entry.pop('__ion__')
        return Ion(**ion_entry)
    else:
        # warnings.warn('Ion not found in database.)
        raise NameError('Ion {} not found in database.'.format(ion_name))

if __name__ == "__main__":
    ions = get_db()
    for name in ions.keys():
        print load_ion(name)

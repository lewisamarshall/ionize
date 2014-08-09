import warnings
from Ion import Ion
import sys
from get_db import get_db
from load_nightingale import load_nightingale


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

        if solvation:
            nightingale_function = load_nightingale(ion_name)
        else:
            nightingale_function = None

        if len(ion_entry) == 3:
            return Ion(ion_name.lower(),
                       ion_entry[0], ion_entry[1], ion_entry[2],
                       nightingale_function)
        elif len(ion_entry) == 5:
            return Ion(ion_name.lower(),
                       ion_entry[0], ion_entry[1], ion_entry[2], ion_entry[3],
                       ion_entry[4], nightingale_function)
    else:
        warnings.warn('Ion not found in database. Returning None.')
        return None

if __name__ == "__main__":
    ion_list = get_db()
    for name in ion_list.keys():
        print load_ion(name)

import warnings
from Ion import Ion
import shelve
import os
import sys
from get_db import get_db
from load_nightingale import load_nightingale


def load_ion(ion_name):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.
    Database derived from Peakmaster, with additions from literature.
    """
    ion_list = get_db()
    if ion_name.lower() in ion_list.keys():
        ion_entry = ion_list[ion_name.lower()]
        ion_list.close()

        nightingale = load_nightingale(ion_name.lower())
        if len(ion_entry) == 3:
            return Ion(ion_name.lower(),
                       ion_entry[0], ion_entry[1], ion_entry[2])
        elif len(ion_entry) == 5:
            return Ion(ion_name.lower(),
                       ion_entry[0], ion_entry[1], ion_entry[2], ion_entry[3],
                       ion_entry[4], nightingale)
    else:
        warnings.warn('Ion not found in database. Returning None.')
        ion_list.close()
        return None

if __name__ == "__main__":
    ion_list = get_db()
    for name in ion_list.keys():
        print load_ion(name)
    ion_list.close()

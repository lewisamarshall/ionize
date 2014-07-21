import warnings
from Ion import Ion
import shelve
import os
import sys


def load_ion(ion_name):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.
    Database derived from Peakmaster, with additions from literature.
    """
    if sys.platform=='win32':
        path = os.path.join(os.getcwd(), os.path.dirname(__file__), 'ions_shelve.bin')
    else:
        path = os.path.join(os.getcwd(), os.path.dirname(__file__), 'ions_shelve.db')

    ion_list = shelve.open(path, flag='r')

    if ion_name.lower() in ion_list.keys():
        ion_entry = ion_list[ion_name.lower()]
        return Ion(ion_name.lower(), ion_entry[0], ion_entry[1], ion_entry[2])
    else:
        warnings.warn('Ion not found in database. Returning None.')
        return None

if __name__ == "__main__":
    path = os.path.join(os.getcwd(), os.path.dirname(__file__), 'ions_shelve.bin')
    for name in sorted(shelve.open(path, 'r').keys()):
        print load_ion(name)

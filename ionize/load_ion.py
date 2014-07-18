import warnings
from Ion import Ion
import shelve


def load_ion(ion_name):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.
    Database derived from Peakmaster, with additions from literature.
    """
    ion_list = shelve.open('ionize/ions_shelve', flag='r')

    if ion_name.lower() in ion_list.keys():
        ion_entry = ion_list[ion_name.lower()]
        return Ion(ion_name.lower(), ion_entry[0], ion_entry[1], ion_entry[2])
    else:
        warnings.warn('Ion not found in database. Returning None.')
        return None

if __name__ == "__main__":
    for name in sorted(shelve.open('ions_shelve').keys()):
        print load_ion(name)

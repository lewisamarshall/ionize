import warnings
from Ion import Ion
import shelve


def load_ion(ion_name):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.
    Database derived from Peakmaster, with additions.
    """
    ion_list = shelve.open('ions_shelve')

    if ion_name.lower() in ion_list.keys():
        ion_entry = ion_list[ion_name.lower()]
        return Ion(ion_name.lower(), ion_entry[0], ion_entry[1], ion_entry[2])
    else:
        warnings.warn('Ion not found in database. Returning None.')
        return None

if __name__ == "__main__":
    print load_ion('hydrochloric acid')
    print load_ion('tris')
    print load_ion('sodium')
    print load_ion('acetic acid')
    print load_ion('e-aminocaproic acid')
    print load_ion('caproic acid')
    print load_ion('magnesium')
    print load_ion('alexa fluor 488')
    print load_ion('fluorescein')
    print load_ion('vanillic acid')
    print load_ion('hydrofluoric acid')
    print load_ion('potassium')
    print load_ion('ammonium')
    print load_ion('histidine')

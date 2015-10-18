import re
from get_db import get_db
import collections


def search_ion(searchstring=None, valence_search=None):
    """Print the names of ions in the database that contain the search string.
    """
    ion_list = get_db()

    ions = []

    if searchstring:
        for name in ion_list.keys():
            if not re.search(searchstring, name):
                ion_list.pop(name)

    if valence_search is not None:
        if not isinstance(valence_search, collections.Iterable):
            valence_search = [valence_search]

        for name in ion_list.keys():
            z = ion_list[name]['valence']
            if not any([zp in valence_search for zp in z]):
                ion_list.pop(name)

    ions = sorted([str(ion) for ion in ion_list.keys()])
    return [str(i) for i in ions]

if __name__ == '__main__':
    for ion in search_ion('fluor'):
        print ion

    print '\n'
    print search_ion('fluor', -1)
    print search_ion(z_search=(3))

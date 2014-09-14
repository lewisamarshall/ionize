import re
from get_db import get_db
import collections


def search_ion(searchstring=None, z_search=None):
    """Print the names of ions in the database that contain the search string.
    """
    ion_list = get_db()

    ions = []

    if searchstring:
        for name in ion_list.keys():
            if not re.search(searchstring, name):
                ion_list.pop(name)

    if z_search is not None:
        if isinstance(z_search, collections.Iterable):
            if len(z_search) > 1:
                z_range = range(min(z_search), max(z_search)+1)
            else:
                z_range = z_search
        else:
            z_range = [z_search]

        for name in ion_list.keys():
            z = ion_list[name][0]
            if not any([zp in z_range for zp in z]):
                ion_list.pop(name)

    ions = sorted([str(ion) for ion in ion_list.keys()])
    return [str(i) for i in ions]

if __name__ == '__main__':
    for ion in search_ion('fluor'):
        print ion

    print '\n'
    print search_ion('fluor', -1)
    print search_ion(z_search=(3))

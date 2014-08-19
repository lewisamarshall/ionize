import re
from get_db import get_db


def search_ion(searchstring):
    """Print the names of ions in the database that contain the search string.

    This function always returns None. It only prints the search results.
    """
    ion_list = get_db()

    ions = []

    for name in sorted(ion_list.keys()):
        if re.search(searchstring, name):
            ions.append(name)

    return [str(i) for i in ions]

if __name__ == '__main__':
    for ion in search_ion('fluor'):
        print ion

    print '\n'
    print search_ion('fluor')

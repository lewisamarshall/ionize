import re
import os
from get_db import get_db


def search_ion(searchstring):
    """Print the names of ions in the database that contain the search string.

    This function always returns None. It only prints the search results.
    """
    ion_list = get_db()

    if searchstring.lower() in ion_list.keys():
        print searchstring
    else:
        for name in sorted(ion_list.keys()):
            if re.search(searchstring, name):
                print name
    return None

if __name__ == '__main__':
    search_ion('fluor')

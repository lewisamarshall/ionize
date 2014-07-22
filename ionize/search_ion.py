import re
from load_ion import load_ion
import os
import shelve
from get_db import get_db


def search_ion(namepart):
    ion_list = get_db()

    if namepart.lower() in ion_list.keys():
        print namepart
    else:
        for name in sorted(ion_list.keys()):
            if re.search(namepart, name):
                print name
    ion_list.close()
    return None

if __name__ == '__main__':
    search_ion('fluor')

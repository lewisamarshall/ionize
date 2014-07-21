import re
from load_ion import load_ion
import os
import shelve


def search_ion(namepart):
    path = os.path.dirname(__file__)
    ion_list = shelve.open(os.path.join(path, 'ions_shelve'), flag='r')

    if namepart.lower() in ion_list.keys():
        print namepart
    else:
        for name in sorted(ion_list.keys()):
            if re.search(namepart, name):
                print name
    return None

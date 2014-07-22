import sys
import os
import shelve


def get_db(flag='r'):
    """Opens the ion database and returns it as a shelve.

    By default opens in write-only state. Pass in a different flag for write
    access.
    """
    if sys.platform == 'win32':
        db_name = 'ions_shelve.bin'
    else:
        db_name = 'ions_shelve.db'
    path = os.path.join(os.getcwd(), os.path.dirname(__file__), db_name)

    ion_list = shelve.open(path, flag=flag)
    return ion_list

if __name__ == '__main__':
    ion_list = get_db()
    print len(ion_list), 'ions in database.'
    ion_list.close()

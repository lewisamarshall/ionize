import warnings
from Ion import Ion
from get_db import get_db
from load_nightingale import load_nightingale


def deserialize_ion(serial):
        return Ion(serial['name'],
                   serial['z'],
                   serial['pKa_ref'],
                   serial['absolute_mobility_ref'],
                   serial['dH'],
                   serial['dCp'],
                   serial['nightingale_function'],
                   T=25.0, T_ref=25.0)

def load_ion(ion_name, solvation=True):
    """Return an ion by name from the database.

    load_ion('ion_name') pulls the named ion from the database.

    By default load solvation data for small ions from the Nightingale
    files. Pass in solvation=False to disable.
    """
    ion_name = ion_name.lower()  # Convert to lower case.
    ion_list = get_db()  # Open the database.
    if ion_name in ion_list.keys():
        ion_entry = ion_list[ion_name]

        return deserialize_ion(ion_entry)
    else:
        # warnings.warn('Ion not found in database.)
        raise NameError(ion_name)

if __name__ == "__main__":
    ions = get_db()
    for name in ions.keys():
        print load_ion(name)

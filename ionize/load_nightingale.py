from scipy import interpolate
from viscosity import viscosity
import sys
import os


def load_nightingale(name):
    """Return the Nightingale mobility function if available for the ion named.

    The function includes data for small ions, including hydrochloric acid,
    sodium, potassium, lithium, magnesium, perchloric acid, rubidium, cesium,
    calcium silver, and sulfuric acid. This data takes the place of viscosity
    correction for these ions, and also includes emperical effects from changes
    in solvation.

    If no data is available, returns None.
    """
    namedict = {
        'hydrochloric acid': 'Nightingale_Cl_data.txt',
        'sodium': 'Nightingale_Na.txt',
        'potassium': 'Nightingale_K_data.txt',
        'lithium': 'Nightingale_Li_data.txt',
        'magnesium': 'Nightingale_Mg_data.txt',
        'perchloric acid': 'Nightingale_Perchlorate_data.txt',
        'rubidium': 'Nightingale_Rb and Cs _data.txt',
        'cesium': 'Nightingale_Rb and Cs _data.txt',
        'calcium': 'Nightingale_Rb and Cs _data.txt',
        'silver': 'Nightingale_Ag_data.txt',
        'sulfuric acid': 'Nightingale_Sulfate_data.txt'
    }

    z_dict = {
        'hydrochloric acid': -1.0,
        'sodium': 1.0,
        'potassium': 1.0,
        'lithium': 1.0,
        'magnesium': 2.0,
        'perchloric acid': -1.0,
        'rubidium': 1.0,
        'cesium': 1.0,
        'calcium': 2.0,
        'silver': 1.0,
        'sulfuric acid': -2.0
    }

    if name in namedict:
        temp = []
        state = []
        datafilename = os.path.join(os.getcwd(), os.path.dirname(__file__),
                                    'nightingale_data', namedict[name])
        datafile = open(datafilename)
        z = z_dict[name]
        datafile.readline()
        for line in datafile:
            entries = [float(i) for i in line.strip().split(',')]
            temp.append(entries[0])
            # Convert from limiting conductivity to mobility
            state.append(entries[1]*10.35e-11/z/viscosity(None, entries[0]))
        statefunc = interpolate.interp1d(temp, state)
        return statefunc
    else:
        return None

if __name__ == "__main__":
    print load_nightingale('silver')(22)

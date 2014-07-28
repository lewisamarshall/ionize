from scipy import interpolate
from viscosity import viscosity
import sys
import os


def load_nightingale(name):
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
                                    'STEEP_files', namedict[name])
        datafile = open(datafilename)
        z = z_dict[name]
        datafile.readline()
        for line in datafile:
            entries = line.strip().split(',')
            entries = map(float, entries)
            temp.append(entries[0])
            state.append(entries[1]*10.35e-21/z/viscosity(None, entries[0])*10**10)
        statefunc = interpolate.interp1d(temp, state)
        return statefunc
    else:
        return None

if __name__ == "__main__":
    print load_nightingale('silver')(22)

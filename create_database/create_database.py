from __future__ import print_function
from math import copysign
import json
import os
import csv
import numpy as np
import pandas as pd

ALIASES = {'chloride': 'hydrochloric acid',
           'dodecylsulfonic acid': 'laurylsulfonic acid'
           }

NIGHTINGALE_FILES = {'silver': 'silver',
                     'calcium': 'calcium',
                     'hydrochloric acid': 'hydrochloric_acid',
                     'potassium': 'potassium',
                     'lithium': 'lithium',
                     'magnesium': 'magnesium',
                     'sodium': 'sodium',
                     'perchlorate acid': 'perchlorate_acid',
                     'rubidium': 'rubidium_cesium',
                     'sulfuric acid': 'sulfuric_acid',
                     'cesium': 'rubidium_cesium'}


class DataCreator(object):

    data = None

    def __init__(self):
        self.load_steep()
        self.load_spresso()
        self.load_nightingale()
        # self.create()
        # self.write

    def load_spresso(self):
        z = pd.read_csv(os.path.join(os.path.dirname(__file__),
                                     'data', 'Z.csv'),
                        engine='python').columns.tolist()
        name = pd.read_csv(os.path.join(os.path.dirname(__file__),
                                         'data', 'name.csv'),
                            engine='python', sep='~', header=None,
                            names=['name'])
        pKa = pd.read_csv(os.path.join(os.path.dirname(__file__),
                                     'data', 'pKa.csv'),
                          engine='python', header=None,
                          names=['pKa {}'.format(i) for i in z])
        mobility = pd.read_csv(os.path.join(os.path.dirname(__file__),
                                            'data', 'mobility.csv'),
                               engine='python', header=None,
                               names=['mobility {}'.format(i) for i in z])
        self.spresso = pd.concat([name, pKa, mobility], axis=1)

    def create(self):
        self.data = dict()
        n_steep = 0
        steep_common = []
        for idx, name in enumerate(namelist):
            # get name
            name = name.rstrip().lower()

            # check to see if the ion is in the steep database.
            if name in steep_db.keys():
                n_steep += 1
                steep_common.append(name)
                dH = steep_db[name][3]
                if steep_db[name][4]:
                    dCp = steep_db[name][4]
                else:
                    dCp = None
            else:
                dH = None
                dCp = None

            # Get pKa lists
            pKa = pKalist.readline()
            pKa = pKa.rstrip().split(',')

            # Get mobility lists
            mobility = mobilitylist.readline()
            mobility = mobility.rstrip().split(',')

            # print name, pKa, mobility
            state = [[], [], [], dH, dCp]
            for i in range(len(z)):
                if not pKa[i] == 'NaN':
                    state[0].append(z[i])
                    state[1].append(float(pKa[i]))
                    state[2].append(copysign(float(mobility[i]), z[i])*1e-9)

            # If steep parameters were added, but the length was different,
            # pad the parameters.
            if state[3] and len(state[3]) < len(state[2]):
                state[3].append(0)
            if state[4] and len(state[4]) < len(state[2]):
                state[4].append(0)

            if name in fit_dict.keys():
                nightingale_function = fit_dict[name]
            else:
                nightingale_function = None
            # Add the result to the ion dictionary.
            serial_ion = {'name': name,
                          'valence': state[0],
                          'reference_pKa': state[1],
                          'reference_mobility': state[2],
                          'enthalpy': dH,
                          'heat_capacity': dCp,
                          'nightingale_data': nightingale_function}

            ion_dict[name] = serial_ion
        ion_dict['taps']['heat_capacity'] = [15.0]
        # print ion_dict['taps']
        # print ion_dict['tartaric acid']
        # boric acid is uniquely in the STEEP database but not the Spresso database
        # add manually.
        steep_db['boric acid'][2][0] *= -1
        boric = steep_db['boric acid']
        ion_dict['boric acid'] = {'name': 'boric acid',
                                  'valence': boric[0],
                                  'reference_pKa': boric[1],
                                  'reference_mobility': boric[2],
                                  'enthalpy': boric[3],
                                  'heat_capacity': boric[4],
                                  'nightingale_data': None}

        # Remove valences to fix ion search
        # TODO: Undo this.
        ion_dict['uranyl'] = ion_dict.pop('uranyl(vi)')
        ion_dict['iron'] = ion_dict.pop('iron(ii)')

        n_steep += 1
        steep_common.append('boric acid')

        for key in ion_dict.keys():
            ion_dict[key]['__ion__'] = True





    def load_nightingale(self):
        nightingale_dict = dict()
        self.nightingale_fits = nightingale_fits = dict()

        for ion in NIGHTINGALE_FILES.keys():
            fullpath = os.path.join(os.path.dirname(__file__),
                                    'data', 'nightingale',
                                    '{}.txt'.format(NIGHTINGALE_FILES[ion]))
            data = pd.read_csv(fullpath, engine='python', header=0,
                               index_col=0, names=['value'])

            nightingale_fits[ion] = {'fit': np.polyfit(data.index,
                                                       data['value'],
                                                       deg=8
                                                       ).tolist(),
                                     'min': min(data.index),
                                     'max': max(data.index)
                                     }

    def load_steep(self):
        fullpath = os.path.join(os.path.dirname(__file__),
                                'data',
                                'STEEP_Database.csv')
        self.steep = pd.read_csv(fullpath)

    def create_aliases(self):
        pass

    def write(self):
        path = os.path.join(os.path.dirname(__file__), '..',
                            'ionize', 'Database', 'ion_data.json')

        with open(path, 'wb') as ion_file:
            json.dump(self.data, ion_file,
                      sort_keys=True, indent=4, separators=(',', ': '))

    def check(self):
        # Check to make that all of the entries in the steep db were added
        for name in steep_db.keys():
            if name not in steep_common:
                print(name, 'in Steep database not added')

        print(len(ion_dict), 'ions in database')
        print(n_steep, 'ions from steep')


if __name__ == '__main__':
    DataCreator()

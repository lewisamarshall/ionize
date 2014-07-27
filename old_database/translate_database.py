from math import copysign
import shelve
from load_steep_db import load_steep_db

z = open('old_database/z.csv', 'r')
namelist = open('old_database/name.csv', 'r')
pKalist = open('old_database/pKa.csv', 'r')
mobilitylist = open('old_database/mobility.csv', 'r')

z = map(int, z.readline().split(','))

steep_db = load_steep_db()
# print z
ion_dict = dict()
n_steep = 0
steep_common = []
for idx, name in enumerate(namelist):
    # get name
    name = name.rstrip().lower()
    if name in steep_db.keys():
        n_steep += 1
        steep_common.append(name)
        dH = steep_db[name][3]
        if steep_db[name][4]:
            dCp = steep_db[name][4]
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

    if state[3] and len(state[3])<len(state[2]):
        state[3].append(None)
    if state[4] and len(state[4])<len(state[2]):
        state[4].append(None)

    ion_dict[name] = state

ion_dict['boric acid'] = steep_db['boric acid']
n_steep += 1
steep_common.append('boric acid')

ions = shelve.open('ions_shelve')
ions.update(ion_dict)
# print ions
ions.close()
print n_steep, 'ions from steep'
for name in steep_db.keys():
    if name not in steep_common:
        print name


for ion in sorted(ion_dict.keys()):
    if ion_dict[ion][3]:
        if len(ion_dict[ion][3]) != len(ion_dict[ion][2]):
            print ion, ion_dict[ion]

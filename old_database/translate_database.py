from math import copysign
import shelve

z = open('./z.csv', 'r')
namelist = open('./name.csv', 'r')
pKalist = open('./pKa.csv', 'r')
mobilitylist = open('./mobility.csv', 'r')

z = map(int, z.readline().split(','))

print z
ion_dict = dict()
for idx, name in enumerate(namelist):
    # get name
    name = name.rstrip().lower()

    # Get pKa lists
    pKa = pKalist.readline()
    pKa = pKa.rstrip().split(',')

    # Get mobility lists
    mobility = mobilitylist.readline()
    mobility = mobility.rstrip().split(',')

    # print name, pKa, mobility
    state = [[], [], []]
    for i in range(len(z)):
        if not pKa[i] == 'NaN':
            state[0].append(z[i])
            state[1].append(float(pKa[i]))
            state[2].append(copysign(float(mobility[i]), z[i])*1e-9)

    ion_dict[name] = state

# print ion_dict

ions = shelve.open('ions_shelve')
ions.update(ion_dict)
print ions
ions.close()

import numpy
import os


files = [(Ag, 'Nightingale_Ag_data.txt'),
         (Ca, 'Nightingale_Ca_data.txt'),
         (Cl, 'Nightingale_Cl_data.txt'),
         (K, 'Nightingale_K_data.txt'),
         (Li, 'Nightingale_Li_data.txt'),
         (Mg, 'Nightingale_Mg_data.txt'),
         (Na, 'Nightingale_Na.txt'),
         ('Nightingale_Perchlorate_data.txt'),
         'Nightingale_Rb and Cs _data.txt',
         'Nightingale_Sulfate_data.txt'
         ]

nightingale_dict = dict()

for key, file in files:
    print os.path.dirname(__file__)
    fullfile = os.path.join(os.getcwd(),os.path.dirname(__file__),'../ionize/nightingale_data', file)
    open_file = open(fullfile)
    print open_file.readline()
    nightingale_dict[key]

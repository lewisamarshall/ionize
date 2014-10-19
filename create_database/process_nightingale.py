import numpy
import os


files = [('Ag', 'Nightingale_Ag_data.txt'),
         ('Ca', 'Nightingale_Ca_data.txt'),
         ('Cl', 'Nightingale_Cl_data.txt'),
         ('K', 'Nightingale_K_data.txt'),
         ('Li', 'Nightingale_Li_data.txt'),
         ('Mg', 'Nightingale_Mg_data.txt'),
         ('Na', 'Nightingale_Na.txt'),
         ('Perchlorate Acid', 'Nightingale_Perchlorate_data.txt'),
         ('Rubidium', 'Nightingale_Rb and Cs _data.txt'),
         ('Sulfate', 'Nightingale_Sulfate_data.txt')
         ]

nightingale_dict = dict()
fit_dict = dict()

for key, file in files:
    print os.path.dirname(__file__)
    fullfile = os.path.join(os.getcwd(),
                            os.path.dirname(__file__),
                            '../ionize/nightingale_data', file)
    open_file = open(fullfile)
    print open_file.readline()
    nightingale_dict[key] = ([], [])
    for line in open_file:
        line = line.split(',')
        nightingale_dict[key][0].append(float(line[0].strip()))
        nightingale_dict[key][1].append(float(line[1].strip()))

    fit_dict[key] = numpy.polyfit(nightingale_dict[key][0],
                                  nightingale_dict[key][1],
                                  deg=6
                                  )

if __name__ == '__main__':
    from matplotlib import pyplot as plot
    xp = numpy.linspace(0, 100)
    plot.figure()
    for key in nightingale_dict.keys():
        plot.plot(nightingale_dict[key][0],
                  nightingale_dict[key][1],
                  label=key)

        plot.plot(xp, numpy.poly1d(fit_dict[key])(xp), 'k--')

        plot.ylim([0, 1])
        plot.legend()
    plot.show()

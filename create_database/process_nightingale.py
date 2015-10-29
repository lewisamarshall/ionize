import numpy
import os


files = [('silver', 'Nightingale_Ag_data.txt'),
         ('calcium', 'Nightingale_Ca_data.txt'),
         ('hydrochloric acid', 'Nightingale_Cl_data.txt'),
         ('potassium', 'Nightingale_K_data.txt'),
         ('lithium', 'Nightingale_Li_data.txt'),
         ('magnesium', 'Nightingale_Mg_data.txt'),
         ('sodium', 'Nightingale_Na.txt'),
         ('perchlorate acid', 'Nightingale_Perchlorate_data.txt'),
         ('rubidium', 'Nightingale_Rb and Cs _data.txt'),
         ('sulfuric acid', 'Nightingale_Sulfate_data.txt'),
         ('cesium', 'Nightingale_Rb and Cs _data.txt')
         ]

nightingale_dict = dict()
fit_dict = dict()

for key, file in files:
    print os.path.dirname(__file__)
    fullfile = os.path.join(os.getcwd(),
                            os.path.dirname(__file__),
                            './nightingale_data', file)
    open_file = open(fullfile)
    print open_file.readline()
    nightingale_dict[key] = ([], [])
    for line in open_file:
        line = line.split(',')
        nightingale_dict[key][0].append(float(line[0].strip()))
        nightingale_dict[key][1].append(float(line[1].strip()))

    fit_dict[key] = {'fit': numpy.polyfit(nightingale_dict[key][0],
                                          nightingale_dict[key][1],
                                          deg=8
                                          ).tolist(),
                     'min': min(nightingale_dict[key][0]),
                     'max': max(nightingale_dict[key][0])
                     }
    # print fit_dict

if __name__ == '__main__':
    from matplotlib import pyplot as plot
    xp = numpy.linspace(0, 100).tolist()
    plot.figure()
    for key in nightingale_dict.keys():
        plot.plot(nightingale_dict[key][0],
                  nightingale_dict[key][1],
                  label=key)
        xpl = [x for x in xp if
               fit_dict[key]['min'] < x < fit_dict[key]['max']]

        plot.plot(xpl, numpy.poly1d(fit_dict[key]['fit'])(xpl), 'k--')

        plot.ylim([0, 1])
        plot.legend()
    plot.show()

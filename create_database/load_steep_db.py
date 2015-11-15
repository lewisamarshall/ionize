from __future__ import print_function
import csv


def load_steep_db(filename='create_database/STEEP_Database.csv'):
    """Return a dictionary containing the STEEP database.

    Opens the database CSV file and converts to a dictionary.
    """
    data = dict()
    with open(filename, 'rU') as file:
        filereader = csv.reader(file, delimiter=',')
        for idx, row in enumerate(filereader):
            if idx:
                name = row[1].lower()
                z = map(int, [row[i] for i in [15, 10, 5] if row[i]])
                pKa = map(float, [row[i] for i in [17, 12, 7] if row[i]])
                dH = map(float, [row[i] for i in [18, 13, 8] if row[i]]) or None
                dCp = map(float, [row[i] for i in [19, 14, 9] if row[i]]) or None
                mu = map(float, [row[i] for i in [16, 11, 6] if row[i]])
                data[name] = [z, pKa, mu, dH, dCp]
    return data

if __name__ == '__main__':
    data = load_steep_db()
    for ion in data:
        # if ion == 'taps':
            print(ion+':', data[ion])

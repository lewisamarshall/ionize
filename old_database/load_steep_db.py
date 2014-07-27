import csv


def load_steep_db(filename='old_database/STEEP_Database.csv'):
    data = dict()
    with open(filename, 'rU') as file:
        filereader = csv.reader(file, delimiter=',')
        for idx, row in enumerate(filereader):
            if idx:
                name = row[1].lower()
                z = map(int, [row[i] for i in [15, 10, 5] if row[i]])
                pKa = map(float, [row[i] for i in [17, 12, 7] if row[i]])
                dH = map(float, [row[i] for i in [18, 13, 8] if row[i]])
                dCp = map(float, [row[i] for i in [19, 14, 9] if row[i]])
                mu = map(float, [row[i] for i in [16, 11, 6] if row[i]])
                data[name] = [z, pKa, mu, dH, dCp]
    return data

if __name__ == '__main__':
    data = load_steep_db()
    for ion in data:
        print ion+':', data[ion]

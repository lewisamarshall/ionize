import os
import json


def get_db():
    """Opens the ion database and returns it as a dictionary.
    """

    db_name = 'ions_db.json'
    path = os.path.join(os.getcwd(), os.path.dirname(__file__), db_name)

    with open(path, 'r') as fp:
        ions = json.load(fp)

    return ions

if __name__ == '__main__':
    ion_db = get_db()
    print len(ion_db), 'ions in database.'
    print ion_db['hydrochloric acid']

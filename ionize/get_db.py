import os
import json


def get_db():
    """Opens the ion database and returns it as a dictionary.
    """

    db_name = 'ions_db.json'
    path = os.path.join(os.getcwd(), os.path.dirname(__file__), db_name)

    with open(path, 'r') as fp:
        ion_list = json.load(fp)

    return ion_list

if __name__ == '__main__':
    ion_list = get_db()
    print len(ion_list), 'ions in database.'
    print ion_list['hydrochloric acid']

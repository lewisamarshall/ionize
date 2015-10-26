import json
import re
import os
import operator

from ..Ion import Ion


class Database(object):

    _default_source = os.path.join(os.getcwd(),
                                   os.path.dirname(__file__),
                                   'ion_data.json')

    data = property(operator.attrgetter("_data"))
    _data = None

    source = property(operator.attrgetter("_source"))
    _source = None

    def __init__(self, source=None):
        self._source = source or self._default_source
        self._open()

    def _open(self):
        with open(self.source, 'r') as fp:
            self._data = json.load(fp)

    def load(self, name):
        name = name.lower()
        if name in self.data:
            data = {key: value for
                    key, value in self.data[name].items() if key != '__ion__'}
            return Ion(**data)
        else:
            raise NameError('Ion {} not found in database.'.format(name))

    def search(self, searchstring):
        return tuple(sorted([str(key) for key in self.data.keys()
                             if re.search(searchstring, key)]))

    def keys(self):
        return sorted(self.data.keys())

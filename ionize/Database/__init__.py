"""Database module for ionize."""
import json
import re
import os
import operator

from ..Ion import Ion


class Database(object):
    """A database containing ion information."""

    _default_source = os.path.join(os.getcwd(),
                                   os.path.dirname(__file__),
                                   'ion_data.json')

    data = property(operator.attrgetter("_data"))
    _data = None

    source = property(operator.attrgetter("_source"))
    _source = None

    def __init__(self, source=None):
        """Initialize a Database instance."""
        self._source = source or self._default_source
        self._open()

    def _open(self):
        with open(self.source, 'r') as fp:
            self._data = json.load(fp)

    def load(self, name):
        """Return an ion from the database based on the name."""
        if name in self.data:
            name = self.data[name].get('alias_of', name)
            data = {key: value for
                    key, value in self.data[name].items() if key != '__ion__'}
            return Ion(**data)
        elif name.lower() in self.data:
            return self.load(name.lower())
        else:
            raise NameError('Ion {} not found in database.'.format(name))

    def search(self, searchstring):
        """Return each name in the database that matches the searchstring.

        The searchstring can be a regular expression.
        """
        return tuple(sorted([str(key) for key in self.data.keys()
                             if re.search(searchstring, key)]))

    def keys(self):
        """Return the keys of the database as a list."""
        return sorted(self.data.keys())

    def serialize(self):
        """Return a JSON formatted serialization of the database."""
        return json.dumps(self.data)

    def __getitem__(self, key):
        """Return the ion that matches the key."""
        return self.load(key)

    def __iter__(self):
        """Retern a generator that yields each ion in the Database."""
        for key in self.keys():
            yield self[key]

    def __repr__(self):
        return 'Database("{}")'.format(self.source)

    def __str__(self):
        return 'Database: {} entries'.format(len(self.keys()))

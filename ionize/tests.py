from .Ion import Ion
from .Solution import Solution
from .get_db import get_db
from .load_ion import load_ion
from .search_ion import search_ion
from .nucleic_acid import nucleic_acid
import unittest


class TestIon(unittest.TestCase):

    def setUp(self):
        self.db = get_db()

    def test_import(self):
        print 'Importing ions.'
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            ion.molar_conductivity(7, 0.1)


class TestSolution(unittest.TestCase):

    def setUp(self):
        self.hcl = load_ion('hydrochloric acid')
        self.tris = load_ion('tris')

    def test_titration(self):
        print 'Titrating buffer.'
        c_tris = 1.0
        pH_old = 14
        n = 100
        c_hcl_set = [c_tris * i * 2.0 / n for i in range(n)]
        for c_hcl in c_hcl_set:
            buf = Solution([self.tris, self.hcl], [c_tris, c_hcl])
            self.assertTrue(buf.pH < pH_old)

if __name__ == '__main__':
    unittest.main()

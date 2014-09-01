import ionize
import unittest


class TestIon(unittest.TestCase):

    def setUp(self):
        self.db = ionize.get_db()

    def test_import(self):
        print 'Importing ions.'
        for ion_name in self.db.keys():
            ion = ionize.load_ion(ion_name)
            ion.molar_conductivity(7, 0.1)


class TestSolution(unittest.TestCase):

    def setUp(self):
        self.hcl = ionize.load_ion('hydrochloric acid')
        self.tris = ionize.load_ion('tris')

    def test_titration(self):
        print 'Titrating buffer.'
        c_tris = 1.0
        pH_old = 14
        n = 100
        c_hcl_set = [c_tris * i * 2.0 / n for i in range(n)]
        for c_hcl in c_hcl_set:
            buf = ionize.Solution([self.tris, self.hcl], [c_tris, c_hcl])
            self.assertTrue(buf.pH < pH_old)

if __name__ == '__main__':
    unittest.main()

from .Ion import Ion
from .Solution import Solution
from .get_db import get_db
from .load_ion import load_ion
from .search_ion import search_ion
from .nucleic_acid import nucleic_acid
from .deserialize import deserialize

import unittest
import warnings
import numpy


class TestIon(unittest.TestCase):

    def setUp(self):
        self.db = get_db()
        warnings.filterwarnings('ignore')

    def test_import(self):
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)

    def test_properties(self):
        pH_list = [5, 7, 9]
        I_list = [0, .01, .1]
        T_list = [20., 25., 30.]
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            for T in T_list:
                ion_T = ion.set_T(T)
                for pH in pH_list:
                    for I in I_list:
                        ion_T.molar_conductivity(pH, I)
                        ion_T.effective_mobility(pH, I)
                        ion_T.diffusivity(pH)

    def test_name_search(self):
        for ion_name in self.db.keys():
            if ion_name not in search_ion(ion_name):
                print ion_name, search_ion(ion_name)
            self.assertTrue(ion_name in search_ion(ion_name))

    def test_z_search(self):
        for z in range(-2, 2):
            for name in search_ion(z_search=z):
                self.assertTrue(z in load_ion(name).z)

    def test_equality(self):
        hcl = load_ion('hydrochloric acid')
        sol = Solution([hcl], [0.1])
        self.assertEqual(hcl, sol.ions[0])
        self.assertIn(load_ion('hydrochloric acid'), sol.ions,
                      [ion.__dict__ for ion in sol.ions])

    def test_serialize(self):
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            self.assertEqual(ion, deserialize(ion.serialize()))


class TestSolution(unittest.TestCase):

    def setUp(self):
        self.hcl = load_ion('hydrochloric acid')
        self.tris = load_ion('tris')
        self.solutions = [Solution([self.tris, self.hcl], [k, 0.1-k])
                          for k in numpy.linspace(0, 0.1, 20)]

    def test_titration(self):
        c_tris = 0.3
        pH_old = 14
        n = 100
        c_hcl_set = [c_tris * i * 2.0 / n for i in range(n)]
        for c_hcl in c_hcl_set:
            buf = Solution([self.tris, self.hcl], [c_tris, c_hcl])
            self.assertTrue(buf.pH < pH_old)
            self.solutions.append(buf)

    def test_solution_properties(self):
        for buf in self.solutions:
            buf.zone_transfer()
            buf.conductivity()
            buf.transference()

    def test_conservation_functions(self):
        for buf in self.solutions:
            buf.kohlrausch()
            buf.alberty()
            buf.jovin()
            buf.gas()

    def test_water_properties(self):
        water = Solution([], [])
        self.assertAlmostEqual(water.pH, 7.0, 2)
        self.assertAlmostEqual(water.I, 0, 4)


class TestNucleicAcid(unittest.TestCase):

    def test_nucleic_acid(self):
        nucleic_acid()
        [nucleic_acid(n) for n in [10, 100, 1000, 10000, 100000]]

if __name__ == '__main__':
    unittest.main()

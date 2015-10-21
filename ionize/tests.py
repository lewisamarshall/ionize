"""Test module for ionize."""

from .Aqueous import Aqueous
from .Ion import Ion
from .Solution import Solution
from .get_db import get_db
from .load_ion import load_ion
from .search_ion import search_ion
from .nucleic_acid import nucleic_acid
from .deserialize import deserialize

import unittest
import warnings
import numpy as np


class TestAqueous(unittest.TestCase):

    """Test the Aqueous class."""

    def setUp(self):
        """Create an instance and a temperature range to test over."""
        self.temperature_range = np.linspace(10, 90)
        self.aqueous = Aqueous

    def test_dielectric(self):
        """Test that dielectric constant is monotone decreasing."""
        d = 100
        for t in self.temperature_range:
            d_new = self.aqueous.dielectric(t)
            self.assertLess(d_new, d)
            d = d_new

    def test_viscosity(self):
        """Test that viscosity is monotone decreasing."""
        self.assertAlmostEqual(self.aqueous.viscosity(25.),
                               8.9e-4, places=5)
        v = 1
        for t in self.temperature_range:
            v_new = self.aqueous.viscosity(t)
            self.assertLess(v_new, v)
            v = v_new

    def test_dissociation(self):
        """Test that dissociation constant is monotone increasing."""
        k = 0
        self.assertAlmostEqual(self.aqueous.viscosity(25.),
                               8.9e-4, places=5)
        for t in self.temperature_range:
            k_new = self.aqueous.dissociation(t)
            self.assertGreater(k_new, k)
            k = k_new

    def test_debye_huckel(self):
        """Test that Debye-Huckel constant is monotone increasing."""
        dh = 0
        for t in self.temperature_range:
            dh_new = self.aqueous.debye_huckel(t)
            self.assertGreater(dh_new, dh)
            dh = dh_new


class TestIon(unittest.TestCase):

    def setUp(self):
        self.db = get_db()
        warnings.filterwarnings('ignore')

    def test_import(self):
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)

    def test_acidity(self):
        """Test that all acidities are computable."""
        ionic_strength_list = [0, .01, .1]
        temperature_list = [20., 25., 30.]
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            for T in temperature_list:
                for I in ionic_strength_list:
                    ion.pKa(I, T)
                    ion.acidity(I, T)

    def test_properties(self):
        pH_list = [5, 7, 9]
        I_list = [0, .01, .1]
        T_list = [20., 25., 30.]
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            for T in T_list:
                for pH in pH_list:
                    for I in I_list:
                        ion.molar_conductivity(pH, I, T)
                        ion.effective_mobility(pH, I, T)
                        ion.diffusivity(pH, I, T)


    def test_equality(self):
        hcl = load_ion('hydrochloric acid')
        hcl2 = load_ion('hydrochloric acid')
        hcl2.context({'pH': 8, 'ionic_strength': 0.1, 'temperature': 28})
        # sol = Solution([hcl], [0.1])
        self.assertEqual(hcl, hcl2)


    def test_serialize(self):
        for ion_name in self.db.keys():
            ion = load_ion(ion_name)
            self.assertEqual(ion, deserialize(ion.serialize()))

    def test_reorder(self):
        pass

class TestSearch(unittest.TestCase):
    def setUp(self):
        self.db = get_db()
        warnings.filterwarnings('ignore')

    def test_name_search(self):
        for ion_name in self.db.keys():
            if ion_name not in search_ion(ion_name):
                print ion_name, search_ion(ion_name)
            self.assertTrue(ion_name in search_ion(ion_name))

    def test_z_search(self):
        for z in range(-2, 2):
            for name in search_ion(valence_search=z):
                self.assertTrue(z in load_ion(name).valence)




# class TestSolution(unittest.TestCase):
#
#     def setUp(self):
#         self.hcl = load_ion('hydrochloric acid')
#         self.tris = load_ion('tris')
#         self.solutions = [Solution([self.tris, self.hcl], [k, 0.1-k])
#                           for k in np.linspace(0, 0.1, 20)]
#
#     def test_titration(self):
#         c_tris = 0.3
#         pH_old = 14
#         n = 100
#         c_hcl_set = [c_tris * i * 2.0 / n for i in range(n)]
#         for c_hcl in c_hcl_set:
#             buf = Solution([self.tris, self.hcl], [c_tris, c_hcl])
#             self.assertTrue(buf.pH < pH_old)
#             self.solutions.append(buf)
#
#     def test_solution_properties(self):
#         for buf in self.solutions:
#             buf.zone_transfer()
#             buf.conductivity()
#             buf.transference()
#
#     def test_conservation_functions(self):
#         for buf in self.solutions:
#             buf.kohlrausch()
#             buf.alberty()
#             buf.jovin()
#             buf.gas()
#
#     def test_water_properties(self):
#         water = Solution([], [])
#         self.assertAlmostEqual(water.pH, 7.0, 2)
#         self.assertAlmostEqual(water.I, 0, 4)
#
#
# class TestNucleicAcid(unittest.TestCase):
#
#     def test_nucleic_acid(self):
#         nucleic_acid()
#         [nucleic_acid(n) for n in [10, 100, 1000, 10000, 100000]]

if __name__ == '__main__':
    unittest.main()

"""Test module for ionize."""

from .Solvent import Aqueous
from .Ion import Ion
from .PolyIon import NucleicAcid, Peptide
from .IonComplex import Protein
from .Solution import Solution
from .Database import Database
from .deserialize import deserialize
from .__main__ import cli

import unittest
import warnings
import numpy as np
from copy import copy
from click.testing import CliRunner


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
            k_new = self.aqueous.dissociation(0, t)
            self.assertGreater(k_new, k)
            k = k_new

    def test_debye_huckel(self):
        """Test that Debye-Huckel constant is monotone increasing."""
        dh = 0
        for t in self.temperature_range:
            dh_new = self.aqueous.debye_huckel(t)
            self.assertGreater(dh_new, dh)
            dh = dh_new

    def test_bjerrum(self):
        """Test the Bjerrum length."""
        bj = 0
        for t in self.temperature_range:
            bj_new = self.aqueous.bjerrum(t)
            self.assertGreater(bj_new, bj)
            bj = bj_new

    def test_pKs(self):
        """Test the solvent dissociation."""
        self.aqueous.pKs(0.01, 25)


class BaseTestIon(object):
    """Base class for ion tests."""
    def test_serialize(self):
        for ion in self.ions:
            self.assertEqual(ion, deserialize(ion.serialize()),
                             'Deserializing {} failed.'.format(ion_name))

    def test_immutable(self):
        """Test that parts of ion state are immutable."""
        ion = self.database.load('histidine')
        for prop in ion._state:
            with self.assertRaises(AttributeError):
                setattr(ion, prop, None)

    def test_repr(self):
        for name in self.database.keys():
            ion = self.database.load(name)
            self.assertEqual(ion, eval(repr(ion)),
                             'Evaluating repr({}) was malformed.'.format(name))

    def test_separability(self):
        ref = self.database['tris']
        for name in ['tris', 'bis-tris', 'hydrochloric acid']:
            other = self.database[name]
            if name is not 'tris':
                self.assertGreater(ref.separability(other, pH=8), 0)
            else:
                self.assertAlmostEqual(ref.separability(other, pH=8), 0)

    def test_context(self):
        ion = self.database['tris']
        pH, ionic_strength, temperature = \
            ion._resolve_context(None, None, None)
        self.assertEqual(pH, None)
        self.assertAlmostEqual(ionic_strength, 0, 5)
        self.assertEqual(temperature, ion.reference_temperature)

        pH, ionic_strength, temperature = \
            ion._resolve_context(1, 0.1, 0)
        self.assertEqual(pH, 1)
        self.assertEqual(ionic_strength, 0.1)
        self.assertEqual(temperature, 0)

class TestIon(unittest.TestCase):

    def setUp(self):
        self.database = Database()
        warnings.filterwarnings('ignore')

    def test_malformed(self):
        good_prop = list(range(1, 4))
        bad_prop = list(range(1, 3))
        properties = ('valence',
                      'reference_pKa',
                      'reference_mobility',
                      'enthalpy',
                      'heat_capacity')

        base_initializer = {prop: good_prop for prop in properties}

        for prop in properties:
            initializer = copy(base_initializer)
            initializer['name'] = 'bad_{}'.format(prop)
            initializer[prop] = bad_prop
            with self.assertRaises(AssertionError):
                Ion(**initializer)

    def test_bad_sign(self):
        with self.assertRaises(AssertionError):
            Ion('bad_sign', [-1, 1], [5, 5], [1e-9, -1e-9])

    def test_acidity(self):
        """Test that all acidities are computable."""
        ionic_strength_list = [0, .01, .1]
        temperature_list = [20., 25., 30.]
        for name in self.database.keys():
            ion = self.database.load(name)
            for T in temperature_list:
                for I in ionic_strength_list:
                    ion.pKa(I, T)
                    ion.acidity(I, T)

    def test_properties(self):
        pH_list = [5, 7, 9]
        I_list = [0, .01, .1]
        T_list = [20., 25., 30.]
        for ion_name in self.database.keys():
            ion = self.database.load(ion_name)
            for T in T_list:
                for pH in pH_list:
                    for I in I_list:
                        ion.molar_conductivity(pH, I, T)
                        ion.mobility(pH, I, T)
                        ion.diffusivity(pH, I, T)

    def test_equality(self):
        hcl = self.database.load('hydrochloric acid')
        hcl2 = self.database.load('hydrochloric acid')
        hcl2.context({'pH': 8, 'ionic_strength': 0.1, 'temperature': 28})
        sol = Solution([hcl], [0.1])
        self.assertEqual(hcl, hcl2)

    def test_hash(self):
        for ion in self.database.keys():
            ion1 = self.database[ion]
            ion2 = self.database[ion]
            ion1.context({'pH': 8, 'ionic_strength': 0.1, 'temperature': 28})
            ion2.context({'pH': 9, 'ionic_strength': 0.05, 'temperature': 23})
            self.assertTrue(hash(ion1)==hash(ion2),
                            "Ions from database aren't identitcal.",
                            )

    def test_serialize(self):
        for ion_name in self.database.keys():
            ion = self.database.load(ion_name)
            self.assertEqual(ion, deserialize(ion.serialize()),
                             'Deserializing {} failed.'.format(ion_name))

    def test_order(self):
        """Ensures that valence order is enforced."""
        with self.assertRaises(AssertionError):
            Ion('badIon', [3, 1, 2], [0, 0, 0], [0, 0, 0])

    def test_immutable(self):
        """Test that parts of ion state are immutable."""
        ion = self.database.load('histidine')
        for prop in ion._state:
            with self.assertRaises(AttributeError):
                setattr(ion, prop, None)

    def test_repr(self):
        for name in self.database.keys():
            ion = self.database.load(name)
            self.assertEqual(ion, eval(repr(ion)),
                             'Evaluating repr({}) was malformed.'.format(name))

    def test_separability(self):
        ref = self.database['tris']
        for name in ['tris', 'bis-tris', 'hydrochloric acid']:
            other = self.database[name]
            if name is not 'tris':
                self.assertGreater(ref.separability(other, pH=8), 0)
            else:
                self.assertAlmostEqual(ref.separability(other, pH=8), 0)

    def test_context(self):
        ion = self.database['tris']
        pH, ionic_strength, temperature = \
            ion._resolve_context(None, None, None)
        self.assertEqual(pH, None)
        self.assertAlmostEqual(ionic_strength, 0, 5)
        self.assertEqual(temperature, ion.reference_temperature)

        pH, ionic_strength, temperature = \
            ion._resolve_context(1, 0.1, 0)
        self.assertEqual(pH, 1)
        self.assertEqual(ionic_strength, 0.1)
        self.assertEqual(temperature, 0)


class TestDatabase(unittest.TestCase):

    def setUp(self):
        self.database = Database()
        warnings.filterwarnings('ignore')

    def test_import(self):
        [ion for ion in self.database]

    def test_search(self):
        for ion_name in self.database.keys():
            search_result = self.database.search(ion_name)
            self.assertTrue(ion_name in search_result,
                            '{} not in {}'.format(ion_name, search_result))


class TestSolution(unittest.TestCase):

    def setUp(self):
        self.solutions = [Solution(['tris', 'hydrochloric acid'], [k, 0.1-k])
                          for k in np.linspace(0, 0.1, 5)]

    def test_walk_concentration(self):
        """Increase Acid concentration and observe pH decrease."""
        c_tris = 0.3
        pH_old = 14
        n = 100
        c_hcl_set = [c_tris * i * 2.0 / n for i in range(n)]
        for c_hcl in c_hcl_set:
            buf = Solution(['tris', 'hydrochloric acid'], [c_tris, c_hcl])
            self.assertTrue(buf.pH < pH_old)
            self.solutions.append(buf)

    def test_titrate(self):
        """Test pH titration."""
        base = Solution(['tris'], [0.1])
        for pH in (1, 3, 5, 7):
            self.assertAlmostEqual(base.titrate('hydrochloric acid', pH).pH,
                                   pH)

    def test_solution_properties(self):
        for buf in self.solutions:
            buf.conductivity()
            buf.debye()
            buf.buffering_capacity()

    def test_transference(self):
        buf = self.solutions[-2]
        self.assertNotEqual(buf.transference('hydrochloric acid'), 0,
                            'HCl should have a non-zero '
                            'transference number.')
        self.assertNotEqual(buf.transference(Database()['tris']), 0,
                            'Tris should have a non-zero '
                            'transference number.')
        self.assertEqual(buf.transference(Database()['bis-tris']), 0)

    def test_zone_transfer(self):
        buf = self.solutions[-2]
        self.assertNotEqual(buf.zone_transfer('hydrochloric acid'), 0,
                            'HCl should have a non-zero '
                            'transference number.')
        self.assertNotEqual(buf.zone_transfer(Database()['tris']), 0,
                            'Tris should have a non-zero '
                            'transference number.')
        self.assertNotEqual(buf.zone_transfer(Database()['bis-tris']), 0)

    def test_conservation_functions(self):
        for buf in self.solutions:
            buf.kohlrausch()
            buf.alberty()
            buf.jovin()
            buf.gas()

    def test_water_properties(self):
        water = Solution([], [])
        self.assertAlmostEqual(water.pH, 7.0, 2)
        self.assertAlmostEqual(water.ionic_strength, 0, 4)

    def test_equality(self):
        sol = Solution(['tris', 'hydrochloric acid'],
                       [0.1, 0.05])

        sol2 = Solution(['tris', 'hydrochloric acid'],
                        [0.1, 0.05])

        template = 'Solution equality failed. {} != {}'
        self.assertEqual(sol, sol2,
                         template.format(sol.serialize(),
                                         sol2.serialize())
                         )

    def test_lookup(self):
        sol = Solution(['tris', 'hydrochloric acid'],
                       [0.1, 0.05])

        self.assertEqual(sol.concentration('tris'), 0.1,
                         'Failed to find tris by name.')

        self.assertEqual(sol.concentration(sol.ions[0]), 0.1,
                         'Failed to find tris by ion.')

    def test_getitem(self):
        sol = Solution(['tris', 'hydrochloric acid'],
                       [0.1, 0.05])

        self.assertEqual(sol['tris'], sol[Database()['tris']],
                         'Failed to get equivilent items from solution.')

    def test_repr(self):
        sol = Solution(['tris', 'hydrochloric acid'],
                       [0.1, 0.05])

        self.assertEqual(sol, eval(repr(sol)),
                         'Solution malformed by repr.')

    def test_equilibrate_CO2(self):
        sol = Solution().equilibrate_CO2()
        self.assertAlmostEqual(sol.pH, 5.6, 1)

    def test_displace(self):
        sol = Solution(['tris', 'acetic acid'], [0.01, 0.005])
        cycle = sol.displace('tris', 'bis-tris').displace('bis-tris', 'tris')
        self.assertAlmostEqual(sol.pH, cycle.pH, 1)
        # Check that guesses work for solutions where the initial concentration
        # does not converge to the correct value
        sol = Solution(['tris', 'chloride', 'hepes'], [0.01, 0.004, 0.001])
        cycle = sol.displace('chloride', guess=[0.009, 0.004])
        self.assertAlmostEqual(sol.pH, cycle.pH, 0)

    def test_safe(self):
        """Test safe pH evaluation."""
        sol = Solution('chloride', 0.001)
        self.assertFalse(sol.safe(), 'Unbuffered acid evaluated as safe.')
        self.assertTrue(sol.titrate('tris', 8).safe(), 'Tris buffer evaluated as unsafe.')

    def test_moderate(self):
        """Test moderate pH evaluation."""
        sol = Solution('chloride', 0.001)
        self.assertFalse(sol.moderate(), 'Unbuffered acid evaluated as safe.')
        self.assertTrue(sol.titrate('tris', 8).moderate(), 'Tris buffer evaluated as unsafe.')

    def test_hash(self):
        sol1 = Solution('chloride', 0.001)
        sol2 = Solution('chloride', 0.001)
        self.assertEqual(hash(sol1), hash(sol2))

class TestNucleicAcid(unittest.TestCase):

    def test_mobility(self):
       mu = 0
       for n in [10, 100, 1000, 10000, 100000, None]:
           mup = NucleicAcid(size=n).mobility()
           self.assertLess(mup, 0, "Mobility is not negative.")
           self.assertGreater(mu, mup, "Mobility didn't increase with size.")
           mu = mup

    def test_charge(self):
        ch = 0
        for n in [10, 100, 1000, 10000, 100000, None]:
            chp = NucleicAcid(size=n).charge()
            self.assertLess(chp, 0, "Charge is not negative.")
            self.assertGreater(abs(chp), abs(ch), "Charge doesn't increase with size.")
            ch = chp


class TestPeptide(unittest.TestCase):

    # def test_serialize(self):
    #     Peptide('2AVI').serialize()
    #     Peptide(sequence='DTHKSEIAHRFKDLGEEHFKGL'
    #                      'VLIAFSQYLQQCPFDEHVKLVNE').serialize()

    def test_sequence_input(self):
        pep = Peptide(sequence='DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNE')
        pep.mobility(7, 0.01)

    def test_mobility(self):
        avi = Peptide(sequence='RETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDE'
                               'FKADEKKFWGKYLYEIARRHPYFYAPELLYYANKY')
        mu = avi.mobility(0, 0.01)
        for pH in range(1, 12):
            mup = avi.mobility(pH, 0.01)
            self.assertLess(mup, mu,
                            "Mobility didn't decrease with pH.")
            mu = mup

    def test_physical(self):
        avi = Peptide(sequence='HKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPL'
                               'LEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQE')
        avi.radius()
        avi.volume()
        avi.molecular_weight


class TestProtein(unittest.TestCase):
    # def test_download(self):
    #     for name in ['2AVI', '3V03']:
    #         p = Protein(name)
    #         p.molecular_weight
    #         p.mobility(pH=8)
    #         with p.context(Solution()):
    #             p.mobility()

    # def test_membership(self):
    #     for peptide in Protein('2AVI'):
    #         self.assertTrue(isinstance(peptide, Peptide))
    pass

class TestCLI(unittest.TestCase):
    def test_database_cli(self):
        runner = CliRunner()
        result = runner.invoke(cli, ['database'])
        self.assertEqual(result.exit_code, 0)

    def test_ion_cli(self):
        runner = CliRunner()
        result = runner.invoke(cli, ['ion', 'tris'])
        self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()

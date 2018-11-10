import unittest
import os.path
import sys

from reactions.nubase import NUBASE_PATH, Nuclide, Nuclides, Energy
from reactions.system import System


class NuclideTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(NUBASE_PATH) as fh:
            cls.lines = list(fh)

    def test_basic_fields(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(0, n.atomic_number)
        self.assertEqual('613.9 s', str(n.half_life))

    def test_column_widths(self):
        n = Nuclide.load(line=self.lines[16])
        self.assertEqual({
            'atomicNumber':             '003',
            'atomicNumberExtra':        '8W',
            'decayModesAndIntensities': 'IT=100',
            'ensdfArchiveFileYear':     '02',
            'excitationEnergy':         '3562.88    0.10',
            'halfLife':                 '56',
            'halfLifeUnit':             'as',
            'massNumber':               '006',
            'massExcess':               '17649.76     0.10',
            'nuclide':                  '6Lii',
            'reference':                '81Ro02e',
            'spinAndParity':            '0+      T=1',
            'unknown':                  '14',
            'yearOfDiscovery':          '1981',
        }, n._row)

    def test_mass_number(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual(1, n.mass_number)

    def test_json(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertEqual({
            'atomicNumber': 0,
            'massNumber':   1,
            'halfLife':     613.9,
        }, n.json())

    def test_isotopic_abundance_1(self):
        "There should be an abundance for a stable isotope."
        n = Nuclide.load(line=self.lines[1])
        self.assertEqual(99.9885, n.isotopic_abundance)
        n = Nuclide.load(line=self.lines[2])
        self.assertEqual(0.0115, n.isotopic_abundance)

    def test_isotopic_abundance_2(self):
        "The abundance should be zero for unstable isotopes."
        n = Nuclide.load(line=self.lines[3])
        self.assertEqual(0., n.isotopic_abundance)

    def test_signature(self):
        n = Nuclide.load(line=self.lines[15])
        self.assertEqual(('6Li', '0'), n.signature)
        n = Nuclide.load(line=self.lines[16])
        self.assertEqual(('6Li', 'i'), n.signature)

    def test_60Ni(self):
        n = Nuclide.load(line=self.lines[781])
        self.assertEqual('60Ni', n.label)
        self.assertEqual(-64472.5, n.mass_excess_kev)

    def test_59Ni(self):
        n = Nuclide.load(line=self.lines[763])
        self.assertEqual('59Ni', n.label)
        self.assertEqual(-61156.1, n.mass_excess_kev)

    def test_beta_decay_note(self):
        n = Nuclide.load(line=self.lines[0])
        self.assertIn('→β-', n.notes)

    def test_neutron_decay_note(self):
        n = Nuclide.load(line=self.lines[6])
        self.assertIn('→n', n.notes)

    def test_proton_decay_note(self):
        n = Nuclide.load(line=self.lines[11])
        self.assertIn('→p', n.notes)

    def test_electron_capture_note(self):
        n = Nuclide.load(line=self.lines[23])
        self.assertIn('→ε', n.notes)

    def test_internal_transition_note(self):
        n = Nuclide.load(line=self.lines[62])
        self.assertIn('→IT', n.notes)

    def test_internal_transition_note2(self):
        n = Nuclide.load(line=self.lines[16])
        self.assertIn('→IT', n.notes)

    def test_2n_note(self):
        n = Nuclide.load(line=self.lines[9])
        self.assertIn('→2n', n.notes)

    def test_2p_note(self):
        n = Nuclide.load(line=self.lines[17])
        self.assertIn('→2p', n.notes)

    def test_alpha_decay_note(self):
        n = Nuclide.load(line=self.lines[29])
        self.assertIn('→α', n.notes)

    def test_B_p_note(self):
        n = Nuclide.load(line=self.lines[41])
        self.assertIn('→β+p', n.notes)

    def test_B_A_note(self):
        n = Nuclide.load(line=self.lines[41])
        self.assertIn('→β+α', n.notes)

    def test_B_2n_note(self):
        n = Nuclide.load(line=self.lines[89])
        self.assertIn('→β-2n', n.notes)

    def test_B_d_note(self):
        n = Nuclide.load(line=self.lines[14])
        self.assertIn('→β-d', n.notes)

    def test_B_3n_note(self):
        n = Nuclide.load(line=self.lines[126])
        self.assertIn('→β-3n', n.notes)

    def test_B_n_note2(self):
        n = Nuclide.load(line=self.lines[478])
        self.assertIn('→β-n', n.notes)

    def test_SF_note(self):
        n = Nuclide.load(line=self.lines[4901])
        self.assertIn('→SF', n.notes)

    def test_mass(self):
        nuclides = Nuclides.db()
        n0 = nuclides.get(('208Pb', '0'))
        self.assertEqual(193729.016588, n0.mass.mev)


class NuclidesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nuclides = Nuclides.db()

    def test_nuclide(self):
        self.assertEqual('7Li', self.nuclides.get(('7Li', '0')).label)

    def test_90Sr(self):
        n0 = self.nuclides.get(('90Sr', '0'))
        self.assertEqual(90, n0.mass_number)

    def test_Y(self):
        n0 = self.nuclides.get(('90Y', '0'))
        self.assertEqual(90, n0.mass_number)

    def test_isomers(self):
        n = self.nuclides.get(('7Li', '0'))
        self.assertEqual((7, 3), n.numbers)
        ns = self.nuclides.isomers[n.numbers]
        self.assertEqual([('7Li', '0'), ('7Li', 'i')], [n.signature for n in ns])

    def test_trace_isotopes(self):
        n0 = self.nuclides.get(('90Sr', '0'))
        self.assertFalse(n0.is_stable)
        self.assertTrue(n0.in_nature)
        self.assertTrue(n0.is_trace)
        self.assertIn('trace', n0.notes)

import unittest
import os.path
import sys

from lenrmc.nubase import DB_PATH, Nuclide, Nuclides, Combinations


class NuclideTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(DB_PATH) as fh:
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
            'id':                       '006',
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


class NuclidesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nuclides = Nuclides.db()

    def test_nuclide(self):
        self.assertEqual('7Li', self.nuclides.get('7Li').label)

    def test_isomers(self):
        n = self.nuclides.get('7Li')
        self.assertEqual((7, 3), n.numbers)
        ns = self.nuclides.isomers[n.numbers]
        self.assertEqual(['7Li', '7Lii'], [n.label for n in ns])


class CombinationsTest(unittest.TestCase):

    def test_outcomes(self):
        c = Combinations.load(reactants=[(1, '7Li'), (1, '1H')])
        self.assertEqual([
            [(1, 0), (7, 4)],
            [(1, 1), (7, 3)],
            [(1, 2), (7, 2)],
            [(1, 3), (7, 1)],
            [(2, 0), (6, 4)],
            [(2, 1), (6, 3)],
            [(2, 2), (6, 2)],
            [(2, 3), (6, 1)],
            [(3, 0), (5, 4)],
            [(3, 1), (5, 3)],
            [(3, 2), (5, 2)],
            [(3, 3), (5, 1)],
            [(4, 0), (4, 4)],
            [(4, 1), (4, 3)],
            [(4, 2), (4, 2)],
            [(4, 3), (4, 1)],
            [(5, 0), (3, 4)],
            [(5, 1), (3, 3)],
            [(5, 2), (3, 2)],
            [(5, 3), (3, 1)],
            [(6, 0), (2, 4)],
            [(6, 1), (2, 3)],
            [(6, 2), (2, 2)],
            [(6, 3), (2, 1)],
            [(7, 0), (1, 4)],
            [(7, 1), (1, 3)],
            [(7, 2), (1, 2)],
            [(7, 3), (1, 1)]]
        , list(c._outcomes()))
